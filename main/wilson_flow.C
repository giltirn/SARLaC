#include<common.h>
#include<fit.h>
#include<parser.h>
#include<plot.h>
#include<random.h>

using namespace CPSfit;

#define ARGS_MEMBERS \
  ( std::string, file_fmt )	       \
  ( int, bin_size)    \
  ( int, traj_start ) \
  ( int, traj_inc ) \
  ( int, traj_lessthan )

struct Args{
  GENERATE_MEMBERS(ARGS_MEMBERS);

  Args(): file_fmt("wflow.%d"), traj_start(0), traj_inc(1), traj_lessthan(2), bin_size(1){}
};
GENERATE_PARSER(Args, ARGS_MEMBERS);


rawDataCorrelationFunctionD readData(const std::string &file_fmt, int traj_start, int traj_inc, int traj_lessthan){
  int ntraj = (traj_lessthan - traj_start)/traj_inc;
  
  std::map<double, std::vector<double> > dd;

  for(int traj=traj_start; traj < traj_lessthan; traj += traj_inc){
    std::string filename = getFilenameFromFmtString(file_fmt, traj);
    std::cout << "Reading traj " << traj << " with filename "  << filename << std::endl;
    assert(fileExists(filename));

    std::ifstream in(filename);
    assert(!in.fail());

    double t;
    double v;
    while(in >> t >> v){
      std::cout << t << " " << v << std::endl;
      dd[t].push_back(v);
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

  rawDataDistributionOptions::binAllowCropByDefault() = true;
  
  rawDataCorrelationFunctionD data_raw = readData(args.file_fmt, args.traj_start, args.traj_inc, args.traj_lessthan);
  jackknifeCorrelationFunctionD data_j(data_raw.size(), 
				       [&](const int i){ 
					 return jackknifeCorrelationFunctionD::ElementType(data_raw.coord(i), 
											   binResample<jackknifeDistributionD>(data_raw.value(i), args.bin_size)
											   );
				       });

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



  typedef DataSeriesAccessor<jackknifeCorrelationFunctionD, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistributionD> > acc;
  {
    //Plot energy density and t0
    MatPlotLibScriptGenerate plot;
    plot.errorBand(acc(data_j));
    
    InterpPoint p(sqrt_t0*sqrt_t0, 0.3);
    plot.plotData(p);
    plot.write("plot_density.py", "plot_density.pdf"); //t^2 <E(t)>
  }

  {
    //Plot derivative energy density and w0
    MatPlotLibScriptGenerate plot;
    plot.errorBand(acc(data_j_der));
    
    InterpPoint p(w0*w0, 0.3);
    plot.plotData(p);
    plot.write("plot_density_deriv.py", "plot_density_deriv.pdf"); //t d/dt( t^2 <E(t)> )
  }



  std::cout << "Done" << std::endl;
  return 0;
};
