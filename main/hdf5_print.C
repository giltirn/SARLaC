#include<config.h>

#ifdef HAVE_HDF5
#include<iostream>
#include<cassert>
#include<sstream>
#include<array>
#include<distribution.h>
#include<superjackknife.h>
#include<parser.h>
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

  bool spec_sci_fmt_threshold;
  int spec_sci_fmt_threshold_val;
  
  CmdLine(){
    spec_elem = false;
    spec_format = false;
    spec_pub_sfsrc = false;
    spec_pub_sf = false;
    spec_sci_fmt_threshold = false;
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
      }else if(si == "-format"){ //"publication","basic" (default)
	spec_format = true;
	spec_format_val = argv[i+1];
	i += 2;
      }else if(si == "-pub_sfsrc"){
	spec_pub_sfsrc = true;
	spec_pub_sfsrc_val = argv[i+1];
	i+=2;
      }else if(si == "-pub_sf"){
	spec_pub_sf = true;
	std::stringstream ss(argv[i+1]); ss >> spec_pub_sf_val;
	i+=2;
      }else if(si == "-sci_fmt_threshold"){
	spec_sci_fmt_threshold = true;
	std::stringstream ss(argv[i+1]); ss >> spec_sci_fmt_threshold_val;
	i+=2;
      }else{
	error_exit(std::cout << "Unknown argument: " << si << std::endl);
      }
    }
  }

};



GENERATE_ENUM_AND_PARSER( DistributionType, (Jackknife)(JackknifeC)(Raw)(DoubleJackknife)(SuperJackknife)  );

void getTypeInfo(DistributionType &type, int & vector_depth, const std::string &filename){
  HDF5reader rd(filename);
  std::string typestr;
  read(rd, typestr, "distributionType");
  rd.enter("value");
  if(rd.contains("size2")) vector_depth = 2;
  else vector_depth =1;

  if(typestr == "rawDataDistribution<double>"){
    type = Raw;
  }else if(typestr == "jackknifeDistribution<double>"){
    type = Jackknife;
  }else if(typestr == "jackknifeCdistribution<double>"){
    type = JackknifeC;
  }else if(typestr == "doubleJackknifeDistribution<double>"){
    type = DoubleJackknife;
  }else if(typestr == "superJackknifeDistribution<double>"){
    type = SuperJackknife;    
  }else error_exit(std::cout << "getTypeInfo type " << typestr << " unimplemented\n");
}

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
struct ActionPrint{
  void operator()(const std::vector<int> &coord, const D &v) const{
    for(int i=0;i<coord.size();i++) std::cout << coord[i] << " ";
    std::cout << v << std::endl; 
  }
};
template<typename D>
class ActionFilteredPrint{
  std::vector<int> filter;

public:
  ActionFilteredPrint(const std::vector<std::string> &f){
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
	//std::cout << "Coord " << coord << " fail filter " << filter << std::endl;
	return;
      }
    }
    //std::cout << "Coord " << coord << " pass filter " << filter << std::endl;
    
    for(int i=0;i<coord.size();i++) std::cout << coord[i] << " ";
    std::cout << v << std::endl; 
  }
};

template<typename VD,typename D, int depth>
void specVDtype(const std::string &filename, const CmdLine &cmdline){
  //std::cout << "Parsing " << printType<VD>() << std::endl;

  VD data;
  readParamsStandard(data,filename);

  //std::cout << data << std::endl;
  if(cmdline.spec_elem){
    ActionFilteredPrint<D> action(cmdline.spec_elem_vals);
    visitor<VD,ActionFilteredPrint<D>,depth>::go(action, data);  
  }else{
    ActionPrint<D> action;
    visitor<VD,ActionPrint<D>,depth>::go(action, data);
  }    
}

template<typename D>
struct setFormat{
  static inline void doit(const std::string &format, const CmdLine &cmdline){
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
      
      if(cmdline.spec_sci_fmt_threshold){
	printer->setSciFormatThreshold(cmdline.spec_sci_fmt_threshold_val);
      }      
      distributionPrint<D>::printer(printer);
    }else if(format != "basic"){
      error_exit(std::cout << "setFormat: Unknown format " << format << std::endl);
    }
  }
};
template<typename T>
struct setFormat<doubleJackknifeDistribution<T> >{
  static inline void doit(const std::string &format, const CmdLine &cmdline){
    if(format == "publication"){
      error_exit(std::cout << "setFormat: Double-jackknife does not support format " << format << std::endl);
    }else if(format != "basic"){
      error_exit(std::cout << "setFormat: Unknown format " << format << std::endl);
    }
  }
};


template<typename D>
void specDtype(const std::string &filename, const int vector_depth, const CmdLine &cmdline){
  if(cmdline.spec_format){
    setFormat<D>::doit(cmdline.spec_format_val,cmdline);
  }
  
  switch(vector_depth){
  case 1:
    specVDtype<std::vector<D>,D, 1 >(filename,cmdline); break;
  case 2:
    specVDtype<std::vector<std::vector<D> >,D, 2 >(filename,cmdline); break;
  default:
    error_exit(std::cout << "specDtype(const std::string &filename, const int vector_depth) Does not support vector depths other than 1 or 2, got " << vector_depth << std::endl);
  }
}

void run(const std::string &filename, const DistributionType type, const int vector_depth, const CmdLine &cmdline){
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
    error_exit(std::cout << "run(const DistributionType type, const int vector_depth) unknown type " << type << std::endl);
  }
}

int main(const int argc, const char* argv[]){
  assert(argc >= 2);
  std::string filename = argv[1];

  CmdLine cmdline;
  cmdline.parse(argc,argv);
  
  DistributionType type;
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
