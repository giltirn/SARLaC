#include<config.h>

#ifdef HAVE_HDF5
#include<iostream>
#include<cassert>
#include<sstream>
#include<array>
#include<distribution.h>
#include<superjackknife.h>
#include<parser.h>
#include<plot.h>
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
  
  CmdLine(){
    spec_elem = false;
    spec_format = false;
    spec_pub_sfsrc = false;
    spec_pub_sf = false;
    spec_sci_fmt_threshold = false;
    spec_pub_exp = false;
    spec_round_pow = false;
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
      }else if(si == "-format"){ //"publication","basic" (default), "sample_plot" (plot of samples)
	spec_format = true;
	spec_format_val = argv[i+1];
	i += 2;

	if(spec_format_val == "sample_plot")
	  sample_plot_file_stub = argv[i++];

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

	//--------------------------------------------------------------------------------------------
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
  inline static void go(const Action &action, const std::vector<T> &v){
    std::vector<int> coord(Depth);
    go_recurse(coord, action, v);
  }
  inline static void go_recurse(std::vector<int> coord, const Action &action, const std::vector<T> &v){
    for(int i=0;i<v.size();i++){
      coord[coord.size()-Depth] = i;
      visitor<T,Action,Depth-1>::go_recurse(coord,action,v[i]);
    }
  }
};
template<typename D, typename Action>
struct visitor<D,Action,0>{
  inline static void go_recurse(std::vector<int> coord, const Action &action, const D &v){
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

template<typename VD,typename D, int depth>
void specVDtype(const std::string &filename, const CmdLine &cmdline, formatter<D> &fmt){
  VD data;
  readParamsStandard(data,filename);

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
  void operator()(const std::vector<int> &coord, const D &v){
    for(int i=0;i<coord.size();i++) std::cout << coord[i] << " ";
    std::cout << v << std::endl; 
  }
};

template<typename D>
class formatSamplePlot: public formatter<D>{
  MatPlotLibScriptGenerate plot;
  std::string file_stub;
public:
  formatSamplePlot(const std::string &_file_stub): file_stub(_file_stub){}

  void operator()(const std::vector<int> &coord, const D &v){
    DistributionSampleAccessor<D> acc(v);
    typename MatPlotLibScriptGenerate::handleType h = plot.plotData(acc);

    std::ostringstream os;
    for(int i=0;i<coord.size();i++) os << coord[i] << " ";
    
    plot.setLegend(h,os.str());
  }

  ~formatSamplePlot(){
    plot.createLegend();
    plot.write(file_stub + ".py", file_stub + ".eps");    
  }
};
template<typename T>
class formatSamplePlot<superJackknifeDistribution<T> >: public formatter<superJackknifeDistribution<T> >{
public:
  formatSamplePlot(const std::string &_file_stub){
    error_exit(std::cout << "formatSamplePlot doesn't support superJackknifeDistribution\n");
  }
  void operator()(const std::vector<int> &coord, const superJackknifeDistribution<T> &v){
  }
};


template<typename D>
struct setFormat{
  static inline formatter<D>* doit(const std::string &format, const CmdLine &cmdline){
    if(format == "publication"){
      int nsf = 3;
      SigFigsSource sfsrc = Largest;
            
      if(cmdline.spec_pub_sfsrc){
	if(cmdline.spec_pub_sfsrc_val == "Central") sfsrc = Central;
	else if(cmdline.spec_pub_sfsrc_val == "Error") sfsrc = Error;
	else if(cmdline.spec_pub_sfsrc_val != "Largest") error_exit(std::cout << "setFormat unknown sig.figs. src " << cmdline.spec_pub_sfsrc_val << std::endl);
      }
      if(cmdline.spec_pub_sf){
	nsf = cmdline.spec_pub_sf_val;
      }
      publicationDistributionPrinter<D> *printer = new publicationDistributionPrinter<D>(nsf,sfsrc);

      if(cmdline.spec_round_pow)
	printer->setRoundPower(cmdline.spec_round_pow_val);
      
      if(cmdline.spec_sci_fmt_threshold)
	printer->setSciFormatThreshold(cmdline.spec_sci_fmt_threshold_val);
      
      if(cmdline.spec_pub_exp)
	printer->setExponent(cmdline.spec_pub_exp_val);
      
      distributionPrint<D>::printer(printer);

      return new formatPrint<D>;
    }else if(format == "basic"){
      return new formatPrint<D>;
    }else if(format == "sample_plot"){
      return new formatSamplePlot<D>(cmdline.sample_plot_file_stub);      
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
      return new formatPrint<doubleJackknifeDistribution<T> >;      
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
  }else fmt = new formatPrint<D>;
  
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
  case Raw:
    specDtype<rawDataDistribution<double> >(filename,vector_depth,cmdline);  break;
  case Jackknife:
    specDtype<jackknifeDistribution<double> >(filename,vector_depth,cmdline);  break;
  case JackknifeC:
    specDtype<jackknifeCdistribution<double> >(filename,vector_depth,cmdline);  break;
  case DoubleJackknife:
    specDtype<doubleJackknifeDistribution<double> >(filename,vector_depth,cmdline);  break;
  case SuperJackknife:
    specDtype<superJackknifeDistribution<double> >(filename,vector_depth,cmdline);  break;    
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
