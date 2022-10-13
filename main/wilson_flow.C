#include<common.h>
#include<fit.h>
#include<parser.h>
#include<plot.h>
#include<random.h>

using namespace CPSfit;

//includes_t2 -   true:  data is t^2 E(t) [Grid default]  false:  data is E(t)
#define ARGS_MEMBERS \
  ( std::string, file_fmt )	       \
  ( bool, includes_t2 )		       \
  ( int, bin_size)    \
  ( int, traj_start ) \
  ( int, traj_inc ) \
  ( int, traj_lessthan )

struct Args{
  GENERATE_MEMBERS(ARGS_MEMBERS);

  Args(): file_fmt("wflow.%d"), traj_start(0), traj_inc(1), traj_lessthan(2), bin_size(1), includes_t2(true){}
};
GENERATE_PARSER(Args, ARGS_MEMBERS);

//includes_t2 -   true:  data is t^2 E(t) [Grid default]  false:  data is E(t)
rawDataCorrelationFunctionD readData(const std::string &file_fmt, int traj_start, int traj_inc, int traj_lessthan, bool includes_t2){
  int ntraj = (traj_lessthan - traj_start)/traj_inc;
  
  std::map<double, std::vector<double> > dd;

  for(int traj=traj_start; traj < traj_lessthan; traj += traj_inc){
    std::string filename = getFilenameFromFmtString(file_fmt, traj);
    std::cout << "Reading traj " << traj << " with filename "  << filename << std::endl;
    assert(fileExists(filename));

    std::ifstream in(filename);
    assert(!in.fail());

    std::set<double> t_traj;

    double t;
    double v;
    while(in >> t >> v){
      std::cout << t << " " << v << std::endl;
      if(!includes_t2) v = t*t*v;

      if(t_traj.count(t) != 0){
	double v_prev = dd.find(t)->second.back();
	std::cout << "Found duplicate value on trajectory for t=" << t << ". Previous value " << v_prev << " this val " << v << std::endl;
	assert(fabs(v_prev-v)<1e-12);
      }else{
	dd[t].push_back(v);
	t_traj.insert(t);
      }
    }
  }

  rawDataCorrelationFunctionD data(dd.size());
  int i=0;
  for(auto it = dd.begin(); it != dd.end(); it++){
    std::cout << "t=" << it->first << " read " << it->second.size() << " values" << std::endl;

    assert(it->second.size() == ntraj);
    data.coord(i) = it->first;
    data.value(i) = rawDataDistributionD(ntraj, [&](const int i){ return it->second[i]; });
    ++i;
  }
  return data;
}

struct dval{
  jackknifeDistributionD const* v;
  int i;

  inline double sep(double target) const{ return fabs( v->best()-target ); } 
  
  dval(const int i, const jackknifeDistributionD &j): i(i), v(&j){}
};

//return t at which data_j(t) = target
jackknifeDistributionD compute(const jackknifeCorrelationFunctionD &data_j, double target = 0.3){
  int N = data_j.value(0).size();

  //sort points by closeness to 0.3, exploiting std::set
  auto cmp = [target](const dval &a, const dval &b){ return a.sep(target) < b.sep(target); };
  std::set<dval, decltype(cmp)> vsorted(cmp);
  for(int i=0;i<data_j.size();i++){
    vsorted.insert(dval(i, data_j.value(i)) );
  }
  //Get the two closest points
  dval a = *vsorted.begin();
  dval b = *std::next(vsorted.begin(),1);
    
  //Order them in time
  dval aa = a.i < b.i ? a : b;
  dval bb = a.i < b.i ? b : a;
    
  double aa_t = data_j.coord(aa.i);
  const jackknifeDistributionD &aa_v = data_j.value(aa.i);

  double bb_t = data_j.coord(bb.i);
  const jackknifeDistributionD &bb_v = data_j.value(bb.i);
    
  std::cout << "Two closest points to " << target <<" :  [" << aa_t << ", " << aa_v << "]   [" << bb_t << ", " << bb_v<<"]" << std::endl;
    
  jackknifeDistributionD dE_by_dt = (bb_v - aa_v)/(bb_t - aa_t);
  std::cout << "dE/dt = " << dE_by_dt << std::endl;
    
  jackknifeDistributionD dt_by_dE = pow(dE_by_dt, -1);
  std::cout << "dt/dE = " << dt_by_dE << std::endl;
    
  double t_closest = data_j.coord(a.i);
  const jackknifeDistributionD &v_closest = data_j.value(a.i);

  jackknifeDistributionD dE = jackknifeDistributionD(N, target) - v_closest;

  std::cout << "dE from " << v_closest << " : " << dE << std::endl;

  jackknifeDistributionD result = jackknifeDistributionD(N, t_closest) + dt_by_dE * dE;
    
  return result;
}


class InterpPoint{
  jackknifeDistributionD xx;
  double yy;
public:
  InterpPoint(const jackknifeDistributionD &x, double y): xx(x), yy(y){}

  inline double x(const int i) const{ return xx.best(); }
  inline double y(const int i) const{ return yy; }
  inline double dxm(const int i) const{ return xx.standardError(); }
  inline double dxp(const int i) const{ return xx.standardError(); }
  inline double dym(const int i) const{ return 0; }
  inline double dyp(const int i) const{ return 0; }
 
  inline int size() const{ return 1; }
};

//c0ca in physical units
void compute_a_corrected(const jackknifeDistributionD &val_lat, const jackknifeDistributionD &val_cont, const jackknifeDistributionD &c0ca){
  // w0(a) = w0^c + c0ca a^2
  // w0(a)/a = w0^c/a + c0ca a  = w0^lat
  // c0ca a^2 + w0^c - a w0^lat = 0

  jackknifeDistributionD a_sol1 = ( val_lat + sqrt( val_lat*val_lat - 4*val_cont*c0ca ) ) / (2. * c0ca);
  jackknifeDistributionD a_sol2 = ( val_lat - sqrt( val_lat*val_lat - 4*val_cont*c0ca ) ) / (2. * c0ca);
  jackknifeDistributionD ainv_sol1 = pow(a_sol1,-1);
  jackknifeDistributionD ainv_sol2 = pow(a_sol2,-1);
  
  std::cout << "a^-1 (+) = " << ainv_sol1 << std::endl;
  std::cout << "a^-1 (-) = " << ainv_sol2 << std::endl;
}


void fakeRandom(jackknifeDistributionD &out, double cen, double err, double N){
  gaussianRandom(out, cen, err * sqrt( double(N)/double(N - 1) ) );
  double corr = cen - out.best();
  for(int i=0;i<N;i++) out.sample(i) += corr;
}


//Basic fitting
int main(const int argc, const char** argv){
  Args args;
  if(argc < 2){
    std::ofstream of("template.args");
    (std::cout << "No parameter file provided: writing template to 'template.args' and exiting\n").flush();
    of << args;
    return 1;
  }    

  parse(args, argv[1]);

  bool compute_a_sqrtt0 = false;
  double sqrtt0_cont_cen, sqrtt0_cont_err;

  bool compute_a_w0 = false;
  double w0_cont_cen, w0_cont_err;

  int rng_seed = 1234;

  int i=2;
  while(i < argc){
    std::string sarg = argv[i];
    if(sarg == "-compute_a_sqrtt0"){
      compute_a_sqrtt0 = true;
      sqrtt0_cont_cen = strToAny<double>(argv[i+1]);
      sqrtt0_cont_err = strToAny<double>(argv[i+2]);
      i+=3;
    }else if(sarg == "-compute_a_w0"){
      compute_a_w0 = true;
      w0_cont_cen = strToAny<double>(argv[i+1]);
      w0_cont_err = strToAny<double>(argv[i+2]);
      i+=3;
    }else if(sarg == "-rng_seed"){ //for generating gaussian distributions for w0, sqrtt0 continuum values for lattice spacing determination
      rng_seed = strToAny<int>(argv[i+1]);
      i+=2;
    }else{
      error_exit(std::cout << "Unrecognized argument: " << sarg);
    }
  }

  RNG.initialize(rng_seed);

  rawDataDistributionOptions::binAllowCropByDefault() = true;
  
  rawDataCorrelationFunctionD data_raw = readData(args.file_fmt, args.traj_start, args.traj_inc, args.traj_lessthan, args.includes_t2);
  jackknifeCorrelationFunctionD data_j(data_raw.size(), 
				       [&](const int i){ 
					 return jackknifeCorrelationFunctionD::ElementType(data_raw.coord(i), 
											   binResample<jackknifeDistributionD>(data_raw.value(i), args.bin_size)
											   );
				       });

  //Prepare physical values of sqrtt0, w0 for lattice spacing computation (if requested)
  int nbinned = data_j.value(0).size();
  jackknifeDistributionD sqrtt0_cont(nbinned);
  if(compute_a_sqrtt0){
    fakeRandom(sqrtt0_cont, sqrtt0_cont_cen, sqrtt0_cont_err, nbinned);
    //gaussianRandom(sqrtt0_cont, sqrtt0_cont_cen, sqrtt0_cont_err * sqrt( double(nbinned)/double(nbinned - 1) ) );
    std::cout << "Physical value of t0^1/2 =" << sqrtt0_cont << std::endl;
  }
  jackknifeDistributionD w0_cont(nbinned);
  if(compute_a_w0){
    //gaussianRandom(w0_cont, w0_cont_cen, w0_cont_err * sqrt( double(nbinned)/double(nbinned - 1) ) );
    fakeRandom(w0_cont, w0_cont_cen, w0_cont_err, nbinned);
    std::cout << "Physical value of w0 =" << w0_cont << std::endl;
  }
  
  //t0 is defined through  t^2 < E(t) > |t=t0 = 0.3
  //w0 is defined from t d/dt( t^2 < E(t) > ) |t=w0^2 = 0.3
  jackknifeDistributionD sqrt_t0 = pow( compute(data_j, 0.3), 0.5);  
  std::cout << "t0^1/2 = " << sqrt_t0 << std::endl;

  //Use forwards derivative,   df/dt (t) = [ f(t+dt) - f(t) ]/dt
  jackknifeCorrelationFunctionD data_j_der(data_j.size()-1, [&](const int i){
      double t = data_j.coord(i);
      jackknifeDistributionD deriv = ( data_j.value(i+1) - data_j.value(i) )/( data_j.coord(i+1) - t );
      return jackknifeCorrelationFunctionD::ElementType(t, t*deriv);
    });

  jackknifeDistributionD w0 = pow( compute(data_j_der, 0.3), 0.5);  
  std::cout << "w0 = " << w0 << std::endl;


  //Use forwards derivative,   df/dt (t) = [ f(t+dt) - f(t-dt) ]/[2dt]
  jackknifeCorrelationFunctionD data_j_der_2(data_j.size()-2, [&](const int i){
      //i=0 -> t=1
      int tt = i+1;
      double t = data_j.coord(tt);
      double dt = data_j.coord(tt+1) - t;

      jackknifeDistributionD deriv = ( data_j.value(tt+1) - data_j.value(tt-1) )/2./dt;
      return jackknifeCorrelationFunctionD::ElementType(t, t*deriv);
    });

  jackknifeDistributionD w0_2 = pow( compute(data_j_der_2, 0.3), 0.5);  
  std::cout << "[2nd order approx] w0 = " << w0_2 << std::endl;

  //Compute lattice spacing if desired
  if(compute_a_sqrtt0){
    jackknifeDistributionD ainv = sqrt_t0 / sqrtt0_cont;
    std::cout << "a^{-1} = " << ainv << " by t0^1/2" << std::endl;

    double c0ca_cen = 0.7307 * 0.042;
    double c0ca_err = 0.7307 * 0.014; //just use err on ca
    jackknifeDistributionD c0ca(nbinned);
    //gaussianRandom(c0ca, c0ca_cen, c0ca_err * sqrt( double(nbinned)/double(nbinned - 1) ) );
    fakeRandom(c0ca, c0ca_cen, c0ca_err, nbinned);

    compute_a_corrected(sqrt_t0, sqrtt0_cont, c0ca);
  }

  if(compute_a_w0){
    jackknifeDistributionD ainv = w0 / w0_cont;
    std::cout << "a^{-1} = " << ainv << " by w0 (1st order)" << std::endl;

    ainv = w0_2 / w0_cont;
    std::cout << "a^{-1} = " << ainv << " by w0 (2nd order)" << std::endl;

    double c0ca_cen = 0.8787 * 0.023;
    double c0ca_err = 0.8787 * 0.013; //just use err on ca
    jackknifeDistributionD c0ca(nbinned);
    //gaussianRandom(c0ca, c0ca_cen, c0ca_err * sqrt( double(nbinned)/double(nbinned - 1) ) );
    fakeRandom(c0ca, c0ca_cen, c0ca_err, nbinned);
		   
    compute_a_corrected(w0, w0_cont, c0ca);
  }

  
  //Plot energy density and t0
  typedef DataSeriesAccessor<jackknifeCorrelationFunctionD, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistributionD> > acc;
  {
    MatPlotLibScriptGenerate plot;
    plot.errorBand(acc(data_j));
    
    InterpPoint p(sqrt_t0*sqrt_t0, 0.3);
    plot.plotData(p);
    plot.write("plot_density.py", "plot_density.pdf"); //t^2 <E(t)>
  }

  //Plot derivative energy density and w0
  {
    MatPlotLibScriptGenerate plot;
    plot.errorBand(acc(data_j_der));
    
    InterpPoint p(w0*w0, 0.3);
    plot.plotData(p);
    plot.write("plot_density_deriv.py", "plot_density_deriv.pdf"); //t d/dt( t^2 <E(t)> )
  }



  std::cout << "Done" << std::endl;
  return 0;
};
