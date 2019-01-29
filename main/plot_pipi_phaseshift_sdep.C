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

#define ARGS_MEMBERS \
  ( std::vector<EpipiMeas>, meas )		\
  ( std::string, Epi_file )			\
  ( std::vector<int>, Epi_idx )			\
  ( std::string, ainv_file )			\
  ( std::vector<int>, ainv_idx)			\
  ( int, L )					\
  ( TwistArray, twists )			\
  ( std::string, out_stub )



struct Args{
  GENERATE_MEMBERS(ARGS_MEMBERS);

  Args(): meas(1), Epi_file("Epi_file"), Epi_idx(1,1), L(32), ainv_file("ainv"), ainv_idx(1,0), twists({1,1,1}){}
};
GENERATE_PARSER(Args, ARGS_MEMBERS);


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
  DispersionRelation dispn = DispersionRelation::Continuum;
  bool stationary_pion = false;

  int ii=2;
  while(ii < argc){
    if(std::string(argv[ii]) == "-shift_Epi"){ //experimentally shift Epi by some fraction Epi -> Epi*(1+x)  where user provides x
      shift_Epi = true;
      shift_Epi_x = strToAny<double>(argv[ii+1]);
      ii+=2;
    }else if(std::string(argv[ii]) == "-sin_dispn_reln"){
      dispn = DispersionRelation::Sin;
      ii++;
    }else if(std::string(argv[ii]) == "-sinhsin_dispn_reln"){
      dispn = DispersionRelation::SinhSin;
      ii++;
    }else if(std::string(argv[ii]) == "-stationary_pion"){ //In a G-parity environment one can obtain the stationary pion mass from the isoscalar singlet or from a periodic ensemble. This option asserts that the input pion mass is that of a stationary pion regardless of the BC twists
      stationary_pion = true;
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

  if(ntwist > 0 && !stationary_pion){
    GSLvector p({(double)args.twists[0],(double)args.twists[1],(double)args.twists[2]});
    for(int s=0;s<nsample;s++) mpi.sample(s) = dispersionRelationGetMass(dispn, Epi.sample(s), p, M_PI/args.L);
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
#pragma omp parallel for
    for(int s=0;s<nsample;s++)
      delta.sample(s) = phaseShiftZ(zeta, Epipi.sample(s),mpi.sample(s),args.L,dispn);
    
    jackknifeDistributionD den(nsample, [&](const int s){ return dispersionRelationGetMass(dispn, Epipi.sample(s), zeta.getd(), 2*M_PI/args.L); });
    jackknifeDistributionD gamma = Epipi/den;
    
    std::cout << "gamma = " << gamma << std::endl;

    //Boost energy to CM frame
    jackknifeDistributionD Epipi_CM_phys = Epipi*ainv / gamma;
    
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

  // struct DataAccessor{
  //   int sample_freq;
  //   const correlationFunction<jackknifeDistributionD, jackknifeDistributionD> &data;
    
  //   DataAccessor(const correlationFunction<jackknifeDistributionD, jackknifeDistributionD> &data, int sample_freq=10): data(data), sample_freq(sample_freq){}

  //   int ncurves() const{ return data.size(); }
  //   std::vector<PythonTuple<double> >  curves(const int i) const{
  //     assert(data.coord(i).size() == data.value(i).size());
  //     int nsample = data.value(i).size();
  //     std::vector<PythonTuple<double> > out;
  //     for(int s=0;s<nsample;s+=sample_freq){
  // 	out.push_back(PythonTuple<double>(data.coord(i).sample(s), data.value(i).sample(s)));
  //     }
  //     return out;
  //   }
  //   int nmarkers() const{ return data.size(); }

  //   PythonTuple<double>  markers(const int i) const{
  //     return PythonTuple<double>(data.coord(i).mean(), data.value(i).mean());
  //   } 
  // };

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
  }



  plot.write(args.out_stub + ".py", args.out_stub + ".pdf");

  return 0;
}
