#include<common.h>
#include<fit.h>
#include<parser.h>
#include<plot.h>
#include<random.h>
#include<limits>

using namespace SARLaC;

//includes_t2 -   true:  data is t^2 E(t) [Grid default]  false:  data is E(t)
//adaptive_smearing : The adaptive Wilson flow smearing was used
#define ARGS_MEMBERS \
  ( std::string, file_fmt )	       \
  ( bool, includes_t2 )		       \
  ( bool, adaptive_smearing )	       \
  ( int, bin_size)    \
  ( int, traj_start ) \
  ( int, traj_inc ) \
  ( int, traj_lessthan )

struct Args{
  GENERATE_MEMBERS(ARGS_MEMBERS);

  Args(): file_fmt("wflow.%d"), traj_start(0), traj_inc(1), traj_lessthan(2), bin_size(1), includes_t2(true), adaptive_smearing(false){}
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

inline double linearlyInterpolate(double t_to, double t1, double v1, double t2, double v2){
  //v1 = a*t1+b
  //v2 = a*t2+b
  double a = (v2-v1)/(t2-t1);
  double b = v1 - a*t1;
  return a*t_to + b;
}

double linearlyInterpolate(double t_to, const std::map<double,double> &data){
  assert(data.size() > 2);
  auto u = data.lower_bound(t_to); //returns first element with t>=t_to
  if(u == data.end()){
    //std::cout << "Found no element with t>=" << t_to << std::endl;
    //All t are < t_to;  linearly interpolate between last two points
    auto i1 = data.rbegin();
    auto i2 = std::next(i1);
    return linearlyInterpolate(t_to, i1->first,i1->second,i2->first,i2->second);
  }else if(u->first == t_to){
    //std::cout << "Found element with exact match to t=" << t_to << std::endl;
    //t==t_to
    return u->second; 
  }else{
    //t>t_to
    if(u == data.begin()){ //t is first data point, linearly interpolate between first two points
      //std::cout << "Found element with t=" << u->first << ">" << t_to << " and it is the first element" << std::endl;
      auto i1 = data.begin();
      auto i2 = std::next(i1);
      return linearlyInterpolate(t_to, i1->first,i1->second,i2->first,i2->second);
    }else{ //t is not first data point, linearly interpolate between t and the previous value
      //std::cout << "Found element with t=" << u->first << ">" << t_to << std::endl;
      auto p = std::prev(u);
      return linearlyInterpolate(t_to, p->first,p->second,u->first,u->second);
    }
  }
}

void testLinearlyInterpolate(){
  double a1=2, b1=5;
  double a2=-3, b2=-4;
  
  double t1=1;
  double v1=a1*t1+b1;

  double t2=2;
  double v2=a1*t2+b1;

  double t3=3;
  double v3=a2*t3+b2;

  double t4=4;
  double v4=a2*t4+b2;

  std::map<double,double> data = { {t1,v1},{t2,v2},{t3,v3},{t4,v4} };

  {
    //Check exact match
    assert( linearlyInterpolate(t1,data) == v1 );
    assert( linearlyInterpolate(t2,data) == v2 );
  }
  {
    //Check t below range
    double t = -1;
    double expect = a1*t + b1;
    assert( fabs( linearlyInterpolate(t,data) - expect ) < 1e-7 );
  }
  {
    //Check t above range
    double t = 5;
    double expect = a2*t + b2;
    assert( fabs( linearlyInterpolate(t,data) - expect ) < 1e-7 );
  }
  {
    //Check t within first range
    double t = 1.5;
    double expect = a1*t + b1;
    assert( fabs( linearlyInterpolate(t,data) - expect ) < 1e-7 );
  }
  {
    //Check t within second range
    double t = 3.5;
    double expect = a2*t + b2;
    assert( fabs( linearlyInterpolate(t,data) - expect ) < 1e-7 );
  }

  std::cout << "Test passed" << std::endl;
}



//Read data generated using the adaptive smearing, for which the t values are not the same between samples
//includes_t2 -   true:  data is t^2 E(t) [Grid default]  false:  data is E(t)
rawDataCorrelationFunctionD readDataAdaptive(const std::string &file_fmt, int traj_start, int traj_inc, int traj_lessthan, bool includes_t2){
  int ntraj = (traj_lessthan - traj_start)/traj_inc;
  
  //We determine the minimum dt value and the largest range, then linearly interpolate all samples to integer multiples of that minimum dt
  //across the full range
  std::vector<std::map<double,double> > sdata(ntraj);

  double tmin = std::numeric_limits<double>::max();
  double tmax = std::numeric_limits<double>::lowest();
  double dt_min = std::numeric_limits<double>::max();

  int sample=0;
  for(int traj=traj_start; traj < traj_lessthan; traj += traj_inc){
    std::string filename = getFilenameFromFmtString(file_fmt, traj);
    std::cout << "Reading traj " << traj << " with filename "  << filename << std::endl;
    assert(fileExists(filename));

    std::ifstream in(filename);
    assert(!in.fail());

    double t_prev=-1;
    double t;
    double v;
    while(in >> t >> v){
      std::cout << t << " " << v << std::endl;
      if(!includes_t2) v = t*t*v;

      if(!sdata[sample].count(t)){ //ignore duplicates
	tmin = std::min(t,tmin);
	tmax = std::max(t,tmax);

	if(t_prev!=-1){
	  assert(t_prev < t);
	  double dt = t-t_prev;
	  dt_min = std::min(dt,dt_min);
	}

	sdata[sample][t] = v;
	t_prev = t;
      }

    }
    ++sample;
  }
  
  std::cout << "Determined range " << tmin << " -> " << tmax << " and min step size " << dt_min << std::endl;
  int npoints = int( ceil((tmax - tmin)/dt_min) ) + 1;
  double dt = (tmax - tmin)/(npoints - 1);

  rawDataCorrelationFunctionD data(npoints);
  for(int i=0;i<npoints;i++){
    double t = tmin + i*dt;
    data.coord(i) = t;
    rawDataDistributionD &dist = data.value(i);
    dist.resize(ntraj);
    std::cout << "Interpolating to t=" << t << std::endl;
    for(int s=0;s<ntraj;s++){
      dist.sample(s) = linearlyInterpolate(t, sdata[s]);
    }
    std::cout << "Interpolated to t=" << t << ", got " << dist << std::endl;
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
void compute_a_corrected(const jackknifeDistributionD &val_lat, const jackknifeDistributionD &val_cont, const jackknifeDistributionD &c0ca, const std::string &comment=""){
  // w0(a) = w0^c + c0ca a^2
  // w0(a)/a = w0^c/a + c0ca a  = w0^lat
  // c0ca a^2 + w0^c - a w0^lat = 0

  jackknifeDistributionD a_sol1 = ( val_lat + sqrt( val_lat*val_lat - 4*val_cont*c0ca ) ) / (2. * c0ca);
  jackknifeDistributionD a_sol2 = ( val_lat - sqrt( val_lat*val_lat - 4*val_cont*c0ca ) ) / (2. * c0ca);
  jackknifeDistributionD ainv_sol1 = pow(a_sol1,-1);
  jackknifeDistributionD ainv_sol2 = pow(a_sol2,-1);
  
  //std::cout << "a^-1 (+) = " << ainv_sol1 << std::endl;
  //std::cout << "a^-1 (-) = " << ainv_sol2 << std::endl;  
  std::cout << "a^-1 = " << ainv_sol2 << " using continuum-limit correction " << comment << std::endl;
}

int main(const int argc, const char** argv){
  //testLinearlyInterpolate();

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
  
  rawDataCorrelationFunctionD data_raw;
  if(args.adaptive_smearing) data_raw = readDataAdaptive(args.file_fmt, args.traj_start, args.traj_inc, args.traj_lessthan, args.includes_t2);
  else data_raw = readData(args.file_fmt, args.traj_start, args.traj_inc, args.traj_lessthan, args.includes_t2);

  jackknifeCorrelationFunctionD data_j(data_raw.size(), 
				       [&](const int i){ 
					 return jackknifeCorrelationFunctionD::ElementType(data_raw.coord(i), 
											   binResample<jackknifeDistributionD>(data_raw.value(i), args.bin_size)
											   );
				       });

  //Prepare physical values of sqrtt0, w0 for lattice spacing computation (if requested)
  int nbinned = data_j.value(0).size();
  jackknifeDistributionD sqrtt0_cont(nbinned), sqrtt0_cont_cen_j(nbinned, sqrtt0_cont_cen);
  if(compute_a_sqrtt0){
    sqrtt0_cont = fakeJackknife(sqrtt0_cont_cen, sqrtt0_cont_err, nbinned, RNG, 5e-2);
    std::cout << "Physical value of t0^1/2 =" << sqrtt0_cont << std::endl;
  }
  jackknifeDistributionD w0_cont(nbinned), w0_cont_cen_j(nbinned, w0_cont_cen);
  if(compute_a_w0){
    w0_cont = fakeJackknife(w0_cont_cen, w0_cont_err, nbinned, RNG, 5e-2);
    std::cout << "Physical value of w0 =" << w0_cont << std::endl;
  }
  
  //t0 is defined through  t^2 < E(t) > |t=t0 = 0.3
  //w0 is defined from t d/dt( t^2 < E(t) > ) |t=w0^2 = 0.3
  jackknifeDistributionD sqrt_t0 = pow( compute(data_j, 0.3), 0.5);  
  std::cout << "t0^1/2 = " << sqrt_t0 << std::endl;
  writeParamsStandard(sqrt_t0, "sqrt_t0.hdf5");

  //Use forwards derivative,   df/dt (t) = [ f(t+dt) - f(t) ]/dt
  jackknifeCorrelationFunctionD data_j_der(data_j.size()-1, [&](const int i){
      double t = data_j.coord(i);
      jackknifeDistributionD deriv = ( data_j.value(i+1) - data_j.value(i) )/( data_j.coord(i+1) - t );
      return jackknifeCorrelationFunctionD::ElementType(t, t*deriv);
    });

  jackknifeDistributionD w0 = pow( compute(data_j_der, 0.3), 0.5);  
  std::cout << "w0 = " << w0 << std::endl;
  writeParamsStandard(w0, "w0.hdf5");

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
  writeParamsStandard(w0_2, "w0_2nd_order.hdf5");
    

  //Compute lattice spacing if desired
  if(compute_a_sqrtt0){
    std::cout << "Computing lattice spacing using t0^1/2 with continuum value " << sqrtt0_cont << std::endl;

    jackknifeDistributionD ainv = sqrt_t0 / sqrtt0_cont;
    std::cout << "a^{-1} = " << ainv << std::endl;

    ainv = sqrt_t0 / sqrtt0_cont_cen_j;
    std::cout << "a^{-1} = " << ainv << " (without error)" << std::endl;

    double c0ca_cen = 0.7307 * 0.042;
    double c0ca_err = 0.7307 * 0.014; //just use err on ca
    jackknifeDistributionD c0ca = fakeJackknife(c0ca_cen, c0ca_err, nbinned, RNG, 5e-2);
    compute_a_corrected(sqrt_t0, sqrtt0_cont, c0ca);
    compute_a_corrected(sqrt_t0, sqrtt0_cont_cen_j, c0ca,"(without error)");
  }

  if(compute_a_w0){
    std::cout << "Computing lattice spacing using w0 with continuum value " << w0_cont << std::endl;

    jackknifeDistributionD ainv = w0 / w0_cont;
    std::cout << "a^{-1} = " << ainv << std::endl;

    //ainv = w0_2 / w0_cont;
    //std::cout << "a^{-1} = " << ainv << " by w0 (2nd order) with continuum value " << w0_cont << std::endl;

    ainv = w0 / w0_cont_cen_j;
    std::cout << "a^{-1} = " << ainv << " (without error)" << std::endl;

    double c0ca_cen = 0.8787 * 0.023;
    double c0ca_err = 0.8787 * 0.013; //just use err on ca
    jackknifeDistributionD c0ca = fakeJackknife(c0ca_cen, c0ca_err, nbinned, RNG, 5e-2);
    compute_a_corrected(w0, w0_cont, c0ca);
    compute_a_corrected(w0, w0_cont_cen_j, c0ca, "(without error)");
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

