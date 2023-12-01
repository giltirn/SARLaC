#include<config.h>

#ifdef HAVE_HDF5
#include<iostream>
#include<cassert>
#include<sstream>
#include<array>

#include<distribution.h>
#include<parser.h>
#include<plot.h>

using namespace CPSfit;

//A program to print arrays of distributions stored in the conventional format (cf distribution_IO.h writeParamsStandard/readParamsStandard)

struct CmdLine{
  bool spec_elem;
  std::vector<std::string> spec_elem_vals;

  bool spec_format;
  std::string spec_format_val;

  bool spec_pub_sfsrc;
  std::string spec_pub_sfsrc_val;

  bool spec_pub_sf;
  int spec_pub_sf_val;

  bool spec_round_pow;
  int spec_round_pow_val;
  
  bool spec_sci_fmt_threshold;
  int spec_sci_fmt_threshold_val;

  bool spec_pub_exp;
  int spec_pub_exp_val;

  std::string sample_plot_file_stub;
  
  bool spec_sample_plot_type;
  std::string spec_sample_plot_type_val;
  
  bool spec_noindex;

  bool factor_out;
  double factor_out_factor;

  bool spec_pub_cen_only;

  int format_sample_index;

  CmdLine(){
    spec_elem = false;
    spec_format = false;
    spec_pub_sfsrc = false;
    spec_pub_sf = false;
    spec_sci_fmt_threshold = false;
    spec_pub_exp = false;
    spec_round_pow = false;
    spec_sample_plot_type = false;
    spec_noindex = false;
    factor_out = false;
    spec_pub_cen_only = false;
  }
  void parse(const int argc, const char* argv[]){
    int i = 2;
    while(i < argc){
      std::string si(argv[i]);      
      if(si == "-elem"){
	spec_elem_vals.resize(0);
	std::istringstream buffer(argv[i+1]);
	std::string line;
	while (std::getline(buffer, line, ',')){
	  spec_elem_vals.push_back(line);
	}
	spec_elem = true;
	i += 2;
      }else if(si == "-format"){ 
	/*Options:
	  "publication"   print value with error in parentheses
	  "basic" (default)  print value +- error
	  "sample_plot" produce a plot of sample. Next argument must be filename. Optional plot formats specified below
	  "breakdown" produce an error breakdown of a superMultiDistribution
	  "sample" print just a single sample. Next argument should be sample index
	*/

	spec_format = true;
	spec_format_val = argv[i+1];
	i += 2;

	if(spec_format_val == "sample_plot")
	  sample_plot_file_stub = argv[i++];
	else if(spec_format_val == "sample")
	  format_sample_index = strToAny<int>(argv[i++]);
	
      }else if(si == "-noindex"){ //For publication, sample or basic printing, suppress printing the index
	spec_noindex = true;
	i++;

	//-------------------- Options specific to publication format ---------------------------------
	//Defaults: Prints 3 sig figs from the largest of the central-value/error
	
      }else if(si == "-pub_sfsrc"){ //Specify the "operational-value" (Largest, Central, Error) upon which the other criteria are applied (eg number of sig figs)
	spec_pub_sfsrc = true;
	spec_pub_sfsrc_val = argv[i+1];
	i+=2;
      }else if(si == "-pub_sf"){ //Specify the number of significant figures of the operational value (defined above)
	spec_pub_sf = true;
	std::stringstream ss(argv[i+1]); ss >> spec_pub_sf_val;
	i+=2;
      }else if(si == "-pub_rounding"){ //Specify the power-of-10 at which the result is rounded. Overrides sig figs
	spec_round_pow = true;
	std::stringstream ss(argv[i+1]); ss >> spec_round_pow_val;
	i+=2;
      }else if(si == "-sci_fmt_threshold"){ //Specify the absolute power-of-10 at which we switch from decimal to scientific format
	spec_sci_fmt_threshold = true;
	std::stringstream ss(argv[i+1]); ss >> spec_sci_fmt_threshold_val;
	i+=2;
      }else if(si == "-pub_exp"){ //Specify the exponent of the output. Overrides sci_fmt_threshold. Set to 0 to force decimal format
	spec_pub_exp = true;
	std::stringstream ss(argv[i+1]); ss >> spec_pub_exp_val;
	i+=2;
      }else if(si == "-pub_cen_only"){ //Only print the central value
	spec_pub_cen_only = true;
	i++;

	//-------------------- Options specific to sample_plot format -------------------------------------
	//Defaults: scatter plot
      }else if(si == "-sample_plot_type"){ //specify the type of plot produced (Scatter, Histogram), default: Scatter
	spec_sample_plot_type = true;
	spec_sample_plot_type_val = argv[i+1];
	i+=2;


	//---------------------Transformations of data-----------------------------------------------------------------------
      }else if(si == "-factor_out"){ //factor out a constant for all printed values
	factor_out = true;
	factor_out_factor = strToAny<double>(argv[i+1]);
	i+=2;
      }else{
	error_exit(std::cout << "Unknown argument: " << si << std::endl);
      }
    }
  }

};


template<typename V, typename Action, int Depth>
struct visitor;

template<typename T, typename Action, int Depth>
struct visitor<std::vector<T>, Action, Depth>{
  inline static void go(const Action &action, std::vector<T> &v){
    std::vector<int> coord(Depth);
    go_recurse(coord, action, v);
  }
  inline static void go_recurse(std::vector<int> coord, const Action &action, std::vector<T> &v){
    for(int i=0;i<v.size();i++){
      coord[coord.size()-Depth] = i;
      visitor<T,Action,Depth-1>::go_recurse(coord,action,v[i]);
    }
  }
};
template<typename D, typename Action>
struct visitor<D,Action,0>{
  inline static void go_recurse(std::vector<int> coord, const Action &action, D &v){
    action(coord, v);
  }
};

template<typename D>
struct formatter{
  virtual void operator()(const std::vector<int> &coord, const D &v) = 0;
  virtual ~formatter(){}
};

template<typename D>
struct actionBasic{
  formatter<D> &fmt;

  actionBasic(formatter<D> &_fmt): fmt(_fmt){}
  
  void operator()(const std::vector<int> &coord, const D &v) const{
    fmt(coord,v);
  }
};
template<typename D>
class actionFilter{
  formatter<D> &fmt;
  std::vector<int> filter;

public:
  actionFilter(const std::vector<std::string> &f, formatter<D> &_fmt): fmt(_fmt){
    filter.resize(f.size());
    for(int i=0;i<f.size();i++){
      if(f[i] == "*"){
	filter[i] = -1;
      }else{
	std::stringstream ss; ss << f[i]; ss >> filter[i];
      }
    }    
  }
  
  void operator()(const std::vector<int> &coord, const D &v) const{
    assert(filter.size() == coord.size());
    for(int i=0;i<coord.size();i++){
      if(filter[i] != -1 && filter[i] != coord[i]){
	return;
      }
    }
    fmt(coord,v);
  }
};

//Divide out a factor
template<typename D>
class actionFactorOut{
  double inv_fac;
public:
  actionFactorOut(const double fac): inv_fac(1./fac){}
   
  void operator()(const std::vector<int> &coord, D &v) const{
    v = v * inv_fac;
  }
};


template<typename VD,typename D, int depth>
void specVDtype(const std::string &filename, const CmdLine &cmdline, formatter<D> &fmt){
  VD data;
  readParamsStandard(data,filename);

  //Transformations
  if(cmdline.factor_out){
    actionFactorOut<D> action(cmdline.factor_out_factor);
    visitor<VD,actionFactorOut<D>,depth>::go(action, data); 
  }

  //Printing
  if(cmdline.spec_elem){
    actionFilter<D> action(cmdline.spec_elem_vals,fmt);
    visitor<VD,actionFilter<D>,depth>::go(action, data);  
  }else{
    actionBasic<D> action(fmt);
    visitor<VD,actionBasic<D>,depth>::go(action, data);
  }    
}


template<typename D>
struct formatPrint: public formatter<D>{
  bool print_index;
  formatPrint(bool _print_index = true): print_index(_print_index){}

  void operator()(const std::vector<int> &coord, const D &v){
    if(print_index) for(int i=0;i<coord.size();i++) std::cout << coord[i] << " ";
    std::cout << v << std::endl; 
  }
};

template<typename D>
class formatSamplePlot: public formatter<D>{
  MatPlotLibScriptGenerate plot;
  const CmdLine &cmdline;
public:
  formatSamplePlot(const CmdLine &_cmdline): cmdline(_cmdline){}

  void operator()(const std::vector<int> &coord, const D &v){
    DistributionSampleAccessor<D> acc(v);
    typename MatPlotLibScriptGenerate::handleType h;

    if(cmdline.spec_sample_plot_type){
      if(cmdline.spec_sample_plot_type_val == "Scatter") h = plot.plotData(acc);
      else if(cmdline.spec_sample_plot_type_val == "Histogram") h = plot.histogram(acc);
      else error_exit(std::cout << "formatSamplePlot: Unknown plot type " << cmdline.spec_sample_plot_type_val << std::endl);
    }else h = plot.plotData(acc);
      
    std::ostringstream os;
    for(int i=0;i<coord.size();i++) os << coord[i] << " ";
    
    plot.setLegend(h,os.str());
  }

  ~formatSamplePlot(){
    plot.createLegend();
    const std::string &file_stub = cmdline.sample_plot_file_stub;
    plot.write(file_stub + ".py", file_stub + ".eps");    
  }
};


template<typename T>
class formatSamplePlot<superJackknifeDistribution<T> >: public formatter<superJackknifeDistribution<T> >{
public:
  formatSamplePlot(const CmdLine &cmdline){
    error_exit(std::cout << "formatSamplePlot doesn't support superJackknifeDistribution\n");
  }
  void operator()(const std::vector<int> &coord, const superJackknifeDistribution<T> &v){
  }
};
template<typename T>
class formatSamplePlot<superMultiDistribution<T> >: public formatter<superMultiDistribution<T> >{
public:
  formatSamplePlot(const CmdLine &cmdline){
    error_exit(std::cout << "formatSamplePlot doesn't support superMultiDistribution\n");
  }
  void operator()(const std::vector<int> &coord, const superMultiDistribution<T> &v){
  }
};


template<typename D>
class formatSuperMultiBreakdown: public formatter<D>{
public:
  formatSuperMultiBreakdown(const CmdLine &cmdline){
    error_exit(std::cout << "formatSuperMultiBreakdown only supports superMultiDistribution\n");
  }
  void operator()(const std::vector<int> &coord, const D &v){
  }
};
template<typename T>
class formatSuperMultiBreakdown< superMultiDistribution<T> >: public formatter< superMultiDistribution<T> >{
  const CmdLine &cmdline;
public:
  formatSuperMultiBreakdown(const CmdLine &cmdline): cmdline(cmdline){
  }
  void operator()(const std::vector<int> &coord, const superMultiDistribution<T> &v){
    std::cout << "Coord (";
    for(int i=0;i<coord.size();i++) std::cout << coord[i] << " ";
    std::cout << ") value " << v << std::endl;
    
    const superMultiLayout & layout = v.getLayout();
    std::cout << "Layout contains " << layout.nEnsembles() << std::endl;
    for(int e=0;e<layout.nEnsembles();e++){
      MultiType type = layout.ensType(e);
      generalContainer data_e = v.getEnsembleDistribution(e);
      
      std::cout << "Ensemble " << e << " with tag " << layout.ensTag(e) << " size " << layout.nSamplesEns(e) << " and value ";
      if(type == MultiType::Jackknife) std::cout << data_e.value<jackknifeDistribution<T> >() << std::endl;
      else if(type == MultiType::Bootstrap) std::cout << data_e.value<bootstrapDistribution<T> >() << std::endl;
      else assert(0);
    }
    std::cout << std::endl;
  }
};

template<typename D>
struct formatPrintSample: public formatter<D>{
  bool print_index;
  int sample;
  formatPrintSample(const int sample, bool _print_index = true): sample(sample), print_index(_print_index){}

  void operator()(const std::vector<int> &coord, const D &v){
    if(print_index) for(int i=0;i<coord.size();i++) std::cout << coord[i] << " ";
    std::cout << (sample == -1 ? v.best() : v.sample(sample)) << std::endl; 
  }
};

template<typename D>
struct setFormat{
  static inline formatter<D>* doit(const std::string &format, const CmdLine &cmdline){
    if(format == "publication"){

      //Print only central value
      if(cmdline.spec_pub_cen_only){
	int nsf = 3;
	if(cmdline.spec_pub_sf)
	  nsf = cmdline.spec_pub_sf_val;	

	publicationCenOrErrDistributionPrinter<D> *printer = new publicationCenOrErrDistributionPrinter<D>(Central, nsf);
	
	if(cmdline.spec_round_pow)
	  printer->setRoundPower(cmdline.spec_round_pow_val);
      
	if(cmdline.spec_sci_fmt_threshold)
	  printer->setSciFormatThreshold(cmdline.spec_sci_fmt_threshold_val);
      
	if(cmdline.spec_pub_exp)
	  printer->setExponent(cmdline.spec_pub_exp_val);
      
	distributionPrint<D>::printer(printer);

	return new formatPrint<D>(!cmdline.spec_noindex);
      }
      else{ //Print central value and error
	int nsf = 3;
	SigFigsSource sfsrc = Largest;
            
	if(cmdline.spec_pub_sfsrc){
	  if(cmdline.spec_pub_sfsrc_val == "Central") sfsrc = Central;
	  else if(cmdline.spec_pub_sfsrc_val == "Error") sfsrc = Error;
	  else if(cmdline.spec_pub_sfsrc_val != "Largest") error_exit(std::cout << "setFormat unknown sig.figs. src " << cmdline.spec_pub_sfsrc_val << std::endl);
	}
	if(cmdline.spec_pub_sf)
	  nsf = cmdline.spec_pub_sf_val;

	publicationDistributionPrinter<D> *printer = new publicationDistributionPrinter<D>(nsf,sfsrc);

	if(cmdline.spec_round_pow)
	  printer->setRoundPower(cmdline.spec_round_pow_val);
      
	if(cmdline.spec_sci_fmt_threshold)
	  printer->setSciFormatThreshold(cmdline.spec_sci_fmt_threshold_val);
      
	if(cmdline.spec_pub_exp)
	  printer->setExponent(cmdline.spec_pub_exp_val);
      
	distributionPrint<D>::printer(printer);

	return new formatPrint<D>(!cmdline.spec_noindex);
      }

    }else if(format == "basic"){
      return new formatPrint<D>(!cmdline.spec_noindex);
    }else if(format == "sample_plot"){
      return new formatSamplePlot<D>(cmdline);
    }else if(format == "breakdown"){
      return new formatSuperMultiBreakdown<D>(cmdline);
    }else if(format == "sample"){
      return new formatPrintSample<D>(cmdline.format_sample_index, !cmdline.spec_noindex);
    }else{
      error_exit(std::cout << "setFormat: Unknown format " << format << std::endl);
    }
  }
};
template<typename T>
struct setFormat<doubleJackknifeDistribution<T> >{
  static inline formatter<doubleJackknifeDistribution<T> >* doit(const std::string &format, const CmdLine &cmdline){
    if(format == "publication"){
      error_exit(std::cout << "setFormat: Double-jackknife does not support format " << format << std::endl);
    }else if(format == "basic"){
      return new formatPrint<doubleJackknifeDistribution<T> >(!cmdline.spec_noindex);      
    }else{
      error_exit(std::cout << "setFormat: Unknown/unsupported format " << format << std::endl);
    }
  }
};


template<typename D>
void specDtype(const std::string &filename, const int vector_depth, const CmdLine &cmdline){
  formatter<D> *fmt;
  if(cmdline.spec_format){
    fmt = setFormat<D>::doit(cmdline.spec_format_val,cmdline);
  }else fmt = new formatPrint<D>(!cmdline.spec_noindex);
  
  switch(vector_depth){
  case 1:
    specVDtype<std::vector<D>,D, 1 >(filename,cmdline,*fmt); break;
  case 2:
    specVDtype<std::vector<std::vector<D> >,D, 2 >(filename,cmdline,*fmt); break;
  default:
    error_exit(std::cout << "specDtype(const std::string &filename, const int vector_depth) Does not support vector depths other than 1 or 2, got " << vector_depth << std::endl);
  }

  delete fmt;
}

void run(const std::string &filename, const DistributionTypeEnum type, const int vector_depth, const CmdLine &cmdline){
  switch(type){
  case DistributionTypeEnum::Raw:
    specDtype<rawDataDistribution<double> >(filename,vector_depth,cmdline);  break;
  case DistributionTypeEnum::Jackknife:
    specDtype<jackknifeDistribution<double> >(filename,vector_depth,cmdline);  break;
  case DistributionTypeEnum::JackknifeC:
    specDtype<jackknifeCdistribution<double> >(filename,vector_depth,cmdline);  break;
  case DistributionTypeEnum::DoubleJackknife:
    specDtype<doubleJackknifeDistribution<double> >(filename,vector_depth,cmdline);  break;
  case DistributionTypeEnum::SuperJackknife:
    specDtype<superJackknifeDistribution<double> >(filename,vector_depth,cmdline);  break;    
  case DistributionTypeEnum::SuperMulti:
    specDtype<superMultiDistribution<double> >(filename,vector_depth,cmdline);  break;    
  case DistributionTypeEnum::Bootstrap:
    specDtype<bootstrapDistribution<double> >(filename,vector_depth,cmdline);  break;
  default:
    error_exit(std::cout << "run(const DistributionTypeEnum type, const int vector_depth) unknown type " << type << std::endl);
  }
}

int main(const int argc, const char* argv[]){
  assert(argc >= 2);
  std::string filename = argv[1];

  CmdLine cmdline;
  cmdline.parse(argc,argv);
  
  DistributionTypeEnum type;
  int vector_depth;
  getTypeInfo(type,vector_depth,filename);
  run(filename, type, vector_depth, cmdline);
  return 0;
}


#else

int main(const int argc, const char* argv[]){
  std::cout << "Require HDF5\n";
  return 1;
}

#endif
