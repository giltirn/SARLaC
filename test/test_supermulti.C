#include<distribution.h>
#include<random.h>

using namespace CPSfit;

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


  std::cout << "Done" << std::endl;
  return 0;
}
    
