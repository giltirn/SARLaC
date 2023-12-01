#include<distribution.h>
#include<random.h>

using namespace SARLaC;

template<typename T>
bool equals(const T &a, const T &b, const double tol){
  if(a.size() != b.size()) return false;
  
  for(int i=0;i<iterate<T>::size(a);i++)
    if(!equals(iterate<T>::at(i, a), iterate<T>::at(i, b),tol) ) return false;

  return true;
}
template<>
bool equals<double>(const double &a, const double &b, const double tol){
  return fabs(a - b) <= tol;
}


int main(void){
  RNG.initialize(1234);

  jackknifeDistribution<double> j1(100);
  gaussianRandom(j1, 1., 2.);

  jackknifeDistribution<double> j2(200);
  gaussianRandom(j2, -30., 16.0);

  bootstrapDistribution<double> b1(bootstrapInitType(1000));
  gaussianRandom(b1, -500, 200.);
  b1.best() = b1.mean();

  bootstrapDistribution<double> b2(bootstrapInitType(2000));
  gaussianRandom(b2, 700., 16.0);
  b2.best() = b2.mean();

  {
    //Test superMulti and superJackknife equivalent

    superJackknifeLayout sjl;
    sjl.addEnsemble("e1",100);
    sjl.addEnsemble("e2",200);

    superJackknifeDistribution<double> sj(sjl, 199.21);
    sj.setEnsembleJackknife("e1",j1);
    sj.setEnsembleJackknife("e2",j2);

    superMultiLayout sml;
    sml.addEnsemble("e1",MultiType::Jackknife,100);
    sml.addEnsemble("e2",MultiType::Jackknife,200);

    superMultiDistribution<double> sm(sml, 199.21);
    sm.setEnsembleDistribution("e1",j1);
    sm.setEnsembleDistribution("e2",j2);

    assert(sm.size() == sj.size());
    assert(iterate<superJackknifeDistribution<double> >::size(sj) == iterate<superMultiDistribution<double> >::size(sm) );
  
    for(int i=0;i<iterate<superJackknifeDistribution<double> >::size(sj);i++)
      assert( iterate<superJackknifeDistribution<double> >::at(i, sj) == iterate<superMultiDistribution<double> >::at(i, sm) );

    assert( equals( sm.standardError(), sj.standardError(), 1e-10 ) );
  }

  { 
    //Test reorder
    
    superMultiLayout sml;
    sml.addEnsemble("e1",MultiType::Jackknife,100);
    sml.addEnsemble("e2",MultiType::Bootstrap,1000);

    superMultiDistribution<double> sm(sml, 199.21);
    sm.setEnsembleDistribution("e2",b1);
    sm.setEnsembleDistribution("e1",j1);

    
    superMultiLayout sml2;
    sml2.addEnsemble("e3",MultiType::Jackknife,300);
    sml2.addEnsemble("e2",MultiType::Bootstrap,1000);
    sml2.addEnsemble("e1",MultiType::Jackknife,100);

    sm.setLayout(sml2);
    
    {
      generalContainer ens = sm.getEnsembleDistribution(0);
      assert(ens.is<jackknifeDistribution<double> >());
      const jackknifeDistribution<double> &ens_d = ens.value<jackknifeDistribution<double> >();
      assert(ens_d.size() == 300);
      for(int i=0;i<300;i++) assert(ens_d.sample(i) == 199.21);
    }      
    {
      generalContainer ens = sm.getEnsembleDistribution(1);
      assert(ens.is<bootstrapDistribution<double> >());
      const bootstrapDistribution<double> &ens_d = ens.value<bootstrapDistribution<double> >();
      assert(ens_d.size() == 1000);
      assert(ens_d.best() == 199.21);
      for(int i=0;i<1000;i++) assert(ens_d.sample(i) == b1.sample(i));
    }
    {
      generalContainer ens = sm.getEnsembleDistribution(2);
      assert(ens.is<jackknifeDistribution<double> >());
      const jackknifeDistribution<double> &ens_d = ens.value<jackknifeDistribution<double> >();
      assert(ens_d.size() == 100);
      for(int i=0;i<100;i++) assert(ens_d.sample(i) == j1.sample(i));
    }

  }

  //Test that error propagation is correct if an ensemble is broken into parts
  {
    rawDataDistribution<double> Araw(600), Braw(600);
    gaussianRandom(Araw, 1., 2.);
    gaussianRandom(Braw, 2., 3.);

    //Introduce some correlation between A and B
    Araw = 0.5* ( Araw + Braw/2. );

    jackknifeDistribution<double> Ajack(Araw), Bjack(Braw);
    
    int N1 = 200;
    int N2 = 400;

    rawDataDistribution<double> Araw1(N1, [&](const int i){ return Araw.sample(i); });
    rawDataDistribution<double> Araw2(N2, [&](const int i){ return Araw.sample(N1+i); });

    //Treat as independent measurements
    rawDataDistribution<double> Braw1(N1, [&](const int i){ return Braw.sample(i); });
    rawDataDistribution<double> Braw2(N2, [&](const int i){ return Braw.sample(N1+i); });

    jackknifeDistribution<double> Ajack1(Araw1), Ajack2(Araw2);
    jackknifeDistribution<double> Bjack1(Braw1), Bjack2(Braw2);

    superMultiLayout sml;
    sml.addEnsemble("e1",MultiType::Jackknife,N1);
    sml.addEnsemble("e2",MultiType::Jackknife,N2);

    superMultiDistribution<double> Amulti1(sml, Ajack1.mean());
    Amulti1.setEnsembleDistribution("e1", Ajack1);

    superMultiDistribution<double> Amulti2(sml, Ajack2.mean());
    Amulti2.setEnsembleDistribution("e2", Ajack2);

    superMultiDistribution<double> Bmulti1(sml, Bjack1.mean());
    Bmulti1.setEnsembleDistribution("e1", Bjack1);

    superMultiDistribution<double> Bmulti2(sml, Bjack2.mean());
    Bmulti2.setEnsembleDistribution("e2", Bjack2);
    
    //Appropriate combination for independent measurements is error-weighted avg
    superMultiDistribution<double> Amulti = weightedAvg(std::vector<superMultiDistribution<double> const*>({&Amulti1, &Amulti2}));
    superMultiDistribution<double> Bmulti = weightedAvg(std::vector<superMultiDistribution<double> const*>({&Bmulti1, &Bmulti2}));
    

    //Look at some different combination
    {
      jackknifeDistribution<double> AplusB_jack = Ajack + Bjack;
      superMultiDistribution<double> AplusB_multi = Amulti + Bmulti;
      std::cout << "A+B  jack: " << AplusB_jack << "  multi: " << AplusB_multi << std::endl;
    }
    {
      jackknifeDistribution<double> AtimesB_jack = Ajack * Bjack;
      superMultiDistribution<double> AtimesB_multi = Amulti * Bmulti;
      std::cout << "A*B  jack: " << AtimesB_jack << "  multi: " << AtimesB_multi << std::endl;
    }

    //What about the covariance? Should just be sum of individual covariances as no cross-covariance and they should combine like errors
    {
      double covAB_jack = jackknifeDistribution<double>::covariance(Ajack, Bjack);

      double covAB_multi = jackknifeDistribution<double>::covariance(Amulti.getEnsembleDistribution(0).value< jackknifeDistribution<double> >(),
								     Bmulti.getEnsembleDistribution(0).value< jackknifeDistribution<double> >())
	+
	jackknifeDistribution<double>::covariance(Amulti.getEnsembleDistribution(1).value< jackknifeDistribution<double> >(),
						  Bmulti.getEnsembleDistribution(1).value< jackknifeDistribution<double> >());
      double covAB_multi_func = superMultiDistribution<double>::covariance(Amulti,Bmulti);

      std::cout << "cov(A,B)  jack: " << covAB_jack << "  multi: " << covAB_multi << "  multi(func): " << covAB_multi_func << std::endl;
    }


  }


  std::cout << "Done" << std::endl;
  return 0;
}
    
