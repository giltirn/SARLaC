#include<fit/fit_wrapper/fit_wrapper_freeze.h>
#include<physics/luscherzeta.h>
#include<physics/delta_0_pheno.h>
#include<common.h>
#include<plot.h>

using namespace CPSfit;

typedef std::array<double,3> MomArray;
typedef std::array<int,3> TwistArray;

//Ptot in units of 2pi/L 

#define EPIPI_MEAS_MEMBERS \
  ( std::string, filename )			\
  ( std::vector<int>, idx )			\
  ( MomArray, ptot )				\
  ( int, isospin )

struct EpipiMeas{
  GENERATE_MEMBERS(EPIPI_MEAS_MEMBERS);

  EpipiMeas(): filename("file"), idx(1,0), ptot({0,0,0}), isospin(0){}
};
GENERATE_PARSER(EpipiMeas, EPIPI_MEAS_MEMBERS);

GENERATE_ENUM_AND_PARSER(MpiSource, (EpiDispnReln)(File) );

#define ARGS_MEMBERS \
  ( std::vector<EpipiMeas>, meas )		\
  ( MpiSource, mpi_src)				\
  ( std::string, Epi_file )			\
  ( std::vector<int>, Epi_idx )			\
  ( std::string, mpi_file )			\
  ( std::vector<int>, mpi_idx )			\
  ( std::string, ainv_file )			\
  ( std::vector<int>, ainv_idx)			\
  ( int, L )					\
  ( TwistArray, twists )			\
  ( std::string, out_stub )



struct Args{
  GENERATE_MEMBERS(ARGS_MEMBERS);

  Args(): meas(1), Epi_file("Epi_file"), Epi_idx(1,1), L(32), ainv_file("ainv"), ainv_idx(1,0), twists({1,1,1}), mpi_file(""), mpi_idx(1,0){}
};
GENERATE_PARSER(Args, ARGS_MEMBERS);

//Note mpi_file and mpi_idx only need to be specified if mpi_src == File


correlationFunction<double, jackknifeDistributionD> genColangeloCurve(const int I, const jackknifeDistributionD &mpi, const jackknifeDistributionD &ainv, const double Emax, const int npoint = 100){
  jackknifeDistributionD mpi_phys = ainv*mpi;

  double Emin = 2.1*mpi_phys.mean();
  double delta = (Emax - Emin)/(npoint-1);
  
  correlationFunction<double, jackknifeDistributionD> out;

  for(int i=0;i<npoint;i++){
    double E = Emin + delta*i;
    jackknifeDistributionD delta(mpi.size(), [&](const int s){ return 180./M_PI * PhenoCurveColangelo::compute(E*E, I, mpi_phys.sample(s)); });
    out.push_back(E, delta);
  }
  return out;
}


//deltaE = E_pipi - 2 E_pi  in lab frame
inline double phaseShiftDeltaEZ(const LuscherZeta &zeta, double deltaE, const double m, const int L){
  const GSLvector &d = zeta.getd();
  const GSLvector &l = zeta.getl(); //G-parity twist directions

  double Pcm2 = pow(2*M_PI/L,2) * d.norm2();
  double Ppi2 = pow(M_PI/L,2) * l.norm2();

  double deltaE2 = deltaE*deltaE;

  double k2 = deltaE2/4. + deltaE*sqrt( m*m + Ppi2 ) - Pcm2/4. + Ppi2;
  double q = sqrt(k2) * L/2./M_PI;

  double Epipi = deltaE + 2*sqrt(m*m + Ppi2);
  double gamma = Epipi/sqrt(Epipi*Epipi - Pcm2);

  double delta = -zeta.calcPhi(q,gamma);
  while(delta > M_PI) delta -= M_PI;
  while(delta < -M_PI) delta += M_PI;
  return delta/M_PI * 180;
}

//d= Pcm /(2pi/L) 
inline double phaseShiftDeltaEZ(const double deltaE, const double m, const int L, const std::array<int,3> &twists = {0,0,0}, const std::array<double,3> &d = {0.,0.,0.}){
  LuscherZeta zeta(twists, d);
  return phaseShiftDeltaEZ(zeta,deltaE,m,L);
}



//deltaE2 = E_pipi^2 - 4 E_pi^2  in lab frame
inline double phaseShiftDeltaE2Z(const LuscherZeta &zeta, double deltaE2, const double m, const int L){
  const GSLvector &d = zeta.getd();
  const GSLvector &l = zeta.getl(); //G-parity twist directions

  double Pcm2 = pow(2*M_PI/L,2) * d.norm2();
  double Ppi2 = pow(M_PI/L,2) * l.norm2();

  double k2 = deltaE2/4. - Pcm2/4. + Ppi2;
  double q = sqrt(k2) * L/2./M_PI;
  
  double Epipi = sqrt( deltaE2 + 4*(m*m + Ppi2) );
  double gamma = Epipi/sqrt(Epipi*Epipi - Pcm2);

  double delta = -zeta.calcPhi(q,gamma);
  while(delta > M_PI) delta -= M_PI;
  while(delta < -M_PI) delta += M_PI;
  return delta/M_PI * 180;
}

//d= Pcm /(2pi/L) 
inline double phaseShiftDeltaE2Z(const double deltaE2, const double m, const int L, const std::array<int,3> &twists = {0,0,0}, const std::array<double,3> &d = {0.,0.,0.}){
  LuscherZeta zeta(twists, d);
  return phaseShiftDeltaE2Z(zeta,deltaE2,m,L);
}




int main(const int argc, const char* argv[]){
  Args args;
  if(argc < 2){
    std::ofstream of("template.args");
    (std::cout << "No parameter file provided: writing template to 'template.args' and exiting\n").flush();
    of << args;
    return 1;
  }    
  
  parse(args, argv[1]);

  bool shift_Epi =false;
  double shift_Epi_x;
  DispersionRelation dispn_mpi = DispersionRelation::Continuum;
  DispersionRelation dispn_gamma = DispersionRelation::Continuum;

  bool phaseshift_deltaE = false;
  bool phaseshift_deltaE2 = false;

  int ii=2;
  while(ii < argc){
    if(std::string(argv[ii]) == "-shift_Epi"){ //experimentally shift Epi by some fraction Epi -> Epi*(1+x)  where user provides x
      shift_Epi = true;
      shift_Epi_x = strToAny<double>(argv[ii+1]);
      ii+=2;
    }else if(std::string(argv[ii]) == "-sin_dispn_reln"){
      dispn_mpi = dispn_gamma = DispersionRelation::Sin;
      ii++;
    }else if(std::string(argv[ii]) == "-sinhsin_dispn_reln"){
      dispn_mpi = dispn_gamma = DispersionRelation::SinhSin;
      ii++;
    }else if(std::string(argv[ii]) == "-sin_dispn_reln_mpionly"){
      dispn_mpi = DispersionRelation::Sin;
      ii++;
    }else if(std::string(argv[ii]) == "-phaseshift_deltaE"){
      phaseshift_deltaE = true;
      ii++;
    }else if(std::string(argv[ii]) == "-phaseshift_deltaE2"){
      phaseshift_deltaE2 = true;
      ii++;
    }else{
      error_exit(std::cout << "Unrecognized argument: " << argv[ii] << std::endl);      
    }
  }

  jackknifeDistribution<double> Epi;
  readHDF5file(Epi, args.Epi_file, args.Epi_idx);

  if(shift_Epi){
    std::cout << "Shifting Epi from " << Epi << " to "; 
    Epi = Epi * (1. + shift_Epi_x);
    std::cout << Epi << std::endl;
  }

  int nsample = Epi.size();

  jackknifeDistributionD ainv;
  readHDF5file(ainv, args.ainv_file, args.ainv_idx);

  assert(ainv.size() == nsample);

  int ntwist = args.twists[0] + args.twists[1] + args.twists[2];

  jackknifeDistribution<double> mpi = Epi;

  if(ntwist > 0){
    if(args.mpi_src == MpiSource::EpiDispnReln){
      GSLvector p({(double)args.twists[0],(double)args.twists[1],(double)args.twists[2]});
      for(int s=0;s<nsample;s++) mpi.sample(s) = dispersionRelationGetMass(dispn_mpi, Epi.sample(s), p, M_PI/args.L);
    }else if(args.mpi_src == MpiSource::File){
      std::cout << "Reading m_pi from file " << args.mpi_file << std::endl;
      readHDF5file(mpi, args.mpi_file, args.mpi_idx);
    }
  }

  std::cout << "Epi = " << Epi << std::endl;
  std::cout << "mpi = " << mpi << std::endl;
  std::cout << "a^-1 = " << ainv << std::endl;

  MatPlotLibScriptGenerate plot;

  correlationFunction<jackknifeDistributionD, jackknifeDistributionD> data_I0;
  correlationFunction<jackknifeDistributionD, jackknifeDistributionD> data_I2;
			       
  for(int m=0;m<args.meas.size();m++){
    jackknifeDistribution<double> Epipi;
    readHDF5file(Epipi, args.meas[m].filename, args.meas[m].idx);
    std::cout << "Epipi = " << Epipi << std::endl;
    assert(Epipi.size() == nsample);    
    
    LuscherZeta zeta(args.twists, args.meas[m].ptot);
    jackknifeDistribution<double> delta(nsample);
    jackknifeDistributionD Epipi_CM_phys(nsample);
    if(phaseshift_deltaE){
      jackknifeDistributionD Epi_contm(nsample);
      jackknifeDistributionD Epipi_contm(nsample);

#pragma omp parallel for
      for(int s=0;s<nsample;s++){
	double deltaE = Epipi.sample(s) - 2*Epi.sample(s);
	delta.sample(s) = phaseShiftDeltaEZ(zeta, deltaE,mpi.sample(s),args.L);
	
	double Pcm2 = pow(2*M_PI/args.L,2) * zeta.getd().norm2();
	double Ppi2 = pow(M_PI/args.L,2) * zeta.getl().norm2();

	Epi_contm.sample(s) = sqrt(mpi.sample(s)*mpi.sample(s) + Ppi2);

	Epipi_contm.sample(s) = deltaE + 2*Epi_contm.sample(s);
	Epipi_CM_phys.sample(s) = sqrt( pow(Epipi_contm.sample(s),2) - Pcm2 )*ainv.sample(s);
      }
      
      std::cout << "E_pi(contm) = " << Epi_contm << std::endl;
      std::cout << "E_pipi(contm) = " << Epipi_contm << std::endl;

    }else if(phaseshift_deltaE2){
      jackknifeDistributionD Epi_contm(nsample);
      jackknifeDistributionD Epipi_contm(nsample);

#pragma omp parallel for
      for(int s=0;s<nsample;s++){
	double deltaE2 = pow(Epipi.sample(s),2) - 4*pow(Epi.sample(s),2);
	delta.sample(s) = phaseShiftDeltaE2Z(zeta, deltaE2,mpi.sample(s),args.L);
	
	double Pcm2 = pow(2*M_PI/args.L,2) * zeta.getd().norm2();
	double Ppi2 = pow(M_PI/args.L,2) * zeta.getl().norm2();

	Epi_contm.sample(s) = sqrt(mpi.sample(s)*mpi.sample(s) + Ppi2);

	Epipi_contm.sample(s) = sqrt( deltaE2 + 4*pow(Epi_contm.sample(s),2) );
	Epipi_CM_phys.sample(s) = sqrt( pow(Epipi_contm.sample(s),2) - Pcm2 )*ainv.sample(s);
      }
      
      std::cout << "E_pi(contm) = " << Epi_contm << std::endl;
      std::cout << "E_pipi(contm) = " << Epipi_contm << std::endl;

    }else{
#pragma omp parallel for
      for(int s=0;s<nsample;s++)
	delta.sample(s) = phaseShiftZ(zeta, Epipi.sample(s),mpi.sample(s),args.L,dispn_gamma);
      
      jackknifeDistributionD den(nsample, [&](const int s){ return dispersionRelationGetMass(dispn_gamma, Epipi.sample(s), zeta.getd(), 2*M_PI/args.L); });
      jackknifeDistributionD gamma = Epipi/den;
      
      std::cout << "gamma = " << gamma << std::endl;
      
      //Boost energy to CM frame
      Epipi_CM_phys = Epipi*ainv / gamma;
    }    
  
    std::cout << "Epipi (CM) = " << Epipi_CM_phys << " GeV" << std::endl;
    std::cout << "delta = " << delta << std::endl << std::endl;

    if(args.meas[m].isospin == 0){
      data_I0.push_back(Epipi_CM_phys, delta);
    }else if(args.meas[m].isospin == 2){
      data_I2.push_back(Epipi_CM_phys, delta);
    }else{
      error_exit(std::cout << "Error: Only support I=0/2\n");
    }
  }

  typedef MatPlotLibScriptGenerate::handleType handle;

  struct DataAccessor{
    const correlationFunction<jackknifeDistributionD, jackknifeDistributionD> &data;
    
    DataAccessor(const correlationFunction<jackknifeDistributionD, jackknifeDistributionD> &data): data(data){}

    int ncurves() const{ return data.size(); }
    std::vector<PythonTuple<double> >  curves(const int i) const{
      double cenx = data.coord(i).mean();
      double errx = data.coord(i).standardError();
      double ceny = data.value(i).mean();
      double erry = data.value(i).standardError();

      //I=0 the phase shift is +ve and anticorrelated with the energy hence at x-dx we have y+|dy|   and x+dx we have y-|dy|
      //I=2 the phase shift is -ve and larger energy gives a more negative value hence at x+dx we have y-|dy|  and at x-dx we have y+|dy|
      //Thus the direction of the curve is isospin invariant
	
      std::vector<PythonTuple<double> > out = {
	PythonTuple<double>(cenx-errx,ceny+erry),
	PythonTuple<double>(cenx,ceny),
	PythonTuple<double>(cenx+errx,ceny-erry) };

      return out;
    }
    int nmarkers() const{ return data.size(); }

    PythonTuple<double>  markers(const int i) const{
      return PythonTuple<double>(data.coord(i).mean(), data.value(i).mean());
    } 
  };


  
  typedef DataAccessor Accessor;

  typedef DataSeriesAccessor<correlationFunction<double, jackknifeDistributionD>,
			     ScalarCoordinateAccessor<double>,
			     DistributionPlotAccessor<jackknifeDistributionD> > CurveAccessor;

  if(data_I0.size() > 0){
    Accessor acc(data_I0);
    //handle handle_I0 = plot.plotData(acc, "data_I0");
    handle handle_I0 = plot.errorCurves(acc, "data_I0");


    correlationFunction<double, jackknifeDistributionD> col = genColangeloCurve(0, mpi, ainv, 0.6, 100);
    CurveAccessor cacc(col);
    typename MatPlotLibScriptGenerate::kwargsType kwargs;
    kwargs["alpha"] = 0.6;
    handle curve_I0 = plot.errorBand(cacc,kwargs);
  }
  if(data_I2.size() > 0){
    Accessor acc(data_I2);
    //handle handle_I2 = plot.plotData(acc, "data_I2");
    handle handle_I2 = plot.errorCurves(acc, "data_I2");

    correlationFunction<double, jackknifeDistributionD> col = genColangeloCurve(2, mpi, ainv, 0.6, 100);
    CurveAccessor cacc(col);
    typename MatPlotLibScriptGenerate::kwargsType kwargs;
    kwargs["alpha"] = 0.6;
    handle curve_I2 = plot.errorBand(cacc,kwargs);
  }



  plot.write(args.out_stub + ".py", args.out_stub + ".pdf");

  return 0;
}
