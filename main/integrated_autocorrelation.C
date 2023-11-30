//Compute the integrated autocorrelation time with errors from the bootstrap technique
#include<distribution/utils.h>
#include<common.h>
#include<plot.h>

using namespace CPSfit;

GENERATE_ENUM_AND_PARSER(Format, (Y)(XY) );
//Y: consecutive values    
//XY: a coordinate (e.g. traj number) followed by a value. Coordinate is ignored

#define PARAMS (std::string, file)(Format, format)(int, bin_size)(int, delta_max)

struct Params{
  _GENERATE_MEMBERS(PARAMS);
  Params(): file("values.dat"), bin_size(1), delta_max(12), format(Format::Y){}
};
_GENERATE_PARSER(Params, PARAMS);


int main(const int argc, const char** argv){
  Params args;
  if(argc < 2){
    std::cout << "Require args file. Writing template to params_template.args" << std::endl;
    std::ofstream out("params_template.args");
    out << args;
    return 0;
  }
  RNG.initialize(1234);

  const std::string arg_file = argv[1];
  parse(args, arg_file);

  std::cout << std::setprecision(16) << std::scientific;

  AutoCorrelationOptions opt;
  int arg = 2;
  
  int discard_first = 0;
  int traj_inc = 1; //scale result by measurement frequency if >1
  while(arg < argc){
    std::string sarg = argv[arg];
    if(sarg == "-precomputed_mean"){
      assert(arg < argc-1);
      std::stringstream ss; ss << argv[arg+1]; ss >> opt.precomputed_mean;
      opt.use_precomputed_mean = true;
      std::cout << "Using precomputed mean: " << opt.precomputed_mean << std::endl;
      arg += 2;
    }else if(sarg == "-discard_first"){ //discard the first N samples, e.g. for thermalization
      assert(arg < argc-1);
      std::stringstream ss; ss << argv[arg+1]; ss >> discard_first;
      std::cout << "Discarding first " << discard_first << " samples" << std::endl;
      assert(discard_first >= 0);
      arg += 2;
    }else if(sarg == "-traj_inc"){
      assert(arg < argc-1);
      traj_inc = strToAny<int>(argv[arg+1]);
      arg += 2;
    }else{
      std::cout << "Error: Unknown argument '" << sarg << "'" << std::endl;
      exit(-1);
    }
  }
 
  rawDataDistributionD data;
  {
    std::ifstream fi(args.file);
    double v;
    int i=0; //count of how many numbers encountered (not necessarily just y-values, depends on format)
    int ycount = 0; //count of how many y-values (data) encountered
    while(!fi.eof()){
      fi >> v;
      if(fi.eof()) break;
      assert(!fi.bad());
      if(args.format == Format::Y || (args.format == Format::XY && i%2==1)){
	if(ycount >= discard_first)
	  data.sampleVector().push_back(v);
	++ycount;
      }
      ++i;
    }
  }
  std::cout << "Read " << data.size() << " samples" << std::endl;
  std::cout << "Value " << data << std::endl;

  if(data.size() <= args.delta_max){
    std::cout << "Error: delta_max is larger than the number of samples " << data.size() << std::endl;
    exit(-1);
  }
  
  writeParamsStandard(data, "quantity.hdf5");
  
  //Compute tau_int
  {
    std::vector<bootstrapDistribution<double> > tau_int;

    std::cout << "Computing error bars using bin/resample of products ( v[s] - <v> )( v[s + delta] - <v> ) with bin size " << args.bin_size << std::endl;
    int nbin = (data.size() - args.delta_max)/args.bin_size;
    std::vector<std::vector<int> > rtable = resampleTable(RNG, nbin);
    tau_int = integratedAutocorrelationMulti(args.delta_max, args.bin_size, args.delta_max, data, rtable, opt);
    for(int i=1;i<tau_int.size();i++) tau_int[i] = tau_int[i] * traj_inc;

    for(int d=0;d<=args.delta_max;d++)
      std::cout << d*traj_inc << " " << tau_int[d] << std::endl;

    std::cout << "Generating plot" << std::endl;
    MatPlotLibScriptGenerate plot;
    
    struct Accessor{
      int traj_inc;
      const std::vector<bootstrapDistributionD> &vec;
      Accessor(const std::vector<bootstrapDistributionD> &vec, int traj_inc): vec(vec), traj_inc(traj_inc){}
      
      inline double x(const int i) const{ return i*traj_inc; }
      inline double upper(const int i) const{ return vec[i].confidenceRegion().second; }
      inline double lower(const int i) const{ return vec[i].confidenceRegion().first; }
      inline int size() const{ return vec.size(); }
    };
    
    Accessor acc(tau_int,traj_inc);
    plot.errorBand(acc);
    plot.setXlabel(R"($\Delta_{\rm cut}$)");
    plot.setYlabel(R"($\tau_{\rm int}(\Delta_{\rm cut}$))");

    plot.write("plot_tau_int.py","plot_tau_int.pdf");
  }
  //Also generate the plot of the autocorrelation function
  {
    std::cout << "Plotting autocorrelation function" << std::endl;
    std::vector<double> sep;
    std::vector<double> autocorr;
    for(int i=0;i<data.size()-1;i++){
      sep.push_back(i*traj_inc);
      autocorr.push_back(autocorrelation(i,data));
    }
    MatPlotLibScriptGenerate plot;
  
    struct Acc{
      const std::vector<double>  &vx;  
      const std::vector<double>  &vy;  
      double x(const int i) const{ return vx[i]; }
      double y(const int i) const{ return vy[i]; }
      double dxm(const int i) const{ return 0; }
      double dxp(const int i) const{ return 0; }
      double dym(const int i) const{ return 0; }
      double dyp(const int i) const{ return 0; }
      int size() const{ return vx.size(); }
      Acc(const std::vector<double> &vx, const std::vector<double> &vy): vx(vx),vy(vy){}
    };
    Acc acc(sep,autocorr);
  
    plot.plotData(acc);

    plot.write("plot_autocorr_func.py","plot_autocorr_func.pdf");
  }

  std::cout << "Done" << std::endl;
  return 0;
}

