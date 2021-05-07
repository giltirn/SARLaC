//Compute the integrated autocorrelation time with errors from the bootstrap technique

#include<distribution/utils.h>
#include<common.h>
#include<plot.h>

using namespace CPSfit;

#define PARAMS (std::string, file)(int, bin_size)(int, delta_max)

struct Params{
  _GENERATE_MEMBERS(PARAMS);
  Params(): file("values.dat"), bin_size(1), delta_max(12){}
};
_GENERATE_PARSER(Params, PARAMS);


int main(const int argc, const char** argv){
  Params args;
  if(argc != 2){
    std::cout << "Require args file. Writing template to params_template.args" << std::endl;
    std::ofstream out("params_template.args");
    out << args;
    return 0;
  }
  RNG.initialize(1234);

  const std::string arg_file = argv[1];
  parse(args, arg_file);

  // std::string file = "/home/ckelly/projects/32nt64_MDWF_DSDR_GparityXYZ_fixedRNG_fullanalysis/evo_plots/plaq/values.dat";
  // int bin_size = 1;
  // int delta_max = 12;
  std::cout << std::setprecision(16) << std::scientific;
  
  rawDataDistributionD data;
  {
    std::ifstream fi(args.file);
    double v;
    while(!fi.eof()){
      fi >> v;
      if(fi.eof()) break;
      assert(!fi.bad());
      data.sampleVector().push_back(v);
    }
  }
  std::cout << "Read " << data.size() << " samples" << std::endl;
  std::cout << "Value " << data << std::endl;

  if(data.size() <= args.delta_max){
    std::cout << "Error: delta_max is larger than the number of samples " << data.size() << std::endl;
    exit(-1);
  }
  
  writeParamsStandard(data, "quantity.hdf5");

  std::vector<bootstrapDistribution<double> > tau_int;

  std::cout << "Computing error bars using bin/resample of products ( v[s] - <v> )( v[s + delta] - <v> ) with bin size " << args.bin_size << std::endl;
  int nbin = (data.size() - args.delta_max)/args.bin_size;
  std::vector<std::vector<int> > rtable = resampleTable(RNG, nbin);
  tau_int = integratedAutocorrelationMulti(args.delta_max, args.bin_size, args.delta_max, data, rtable);

  for(int d=0;d<=args.delta_max;d++)
    std::cout << d << " " << tau_int[d] << std::endl;

  std::cout << "Generating plot" << std::endl;
  MatPlotLibScriptGenerate plot;
  
  struct Acc{
    const std::vector<bootstrapDistribution<double> >  &tau;  
    double x(const int i) const{ return i; }
    double y(const int i) const{ return tau[i].best(); }
    double dxm(const int i) const{ return 0; }
    double dxp(const int i) const{ return 0; }
    double dym(const int i) const{ return tau[i].standardError(); }
    double dyp(const int i) const{ return tau[i].standardError(); }
    int size() const{ return tau.size(); }
    Acc(const std::vector<bootstrapDistribution<double> >  &tau): tau(tau){}
  };
  Acc acc(tau_int);
  
  plot.plotData(acc);

  plot.write("plot.py","plot.pdf");

  std::cout << "Done" << std::endl;
  return 0;
}

