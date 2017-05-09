#ifndef _FIT_MPI_GPARITY_AMA_FITFUNC_H
#define _FIT_MPI_GPARITY_AMA_FITFUNC_H

struct FitParams{
  double A;
  double m;
  FitParams(){}
  FitParams(const double _A, const double _m): A(_A), m(_m){}
};
struct FitParamDerivs{
  double dA;
  double dm;

  inline const double &operator()(const int i) const{ return i == 1 ? dm : dA; }
  inline size_t size() const{ return 2;}
};
  

class FitCosh{
  const double Lt;
public:
  typedef double ValueType;
  typedef FitParams ParameterType;
  typedef FitParamDerivs ValueDerivativeType; //derivative wrt parameters
  typedef double GeneralizedCoordinate; //time coord

  FitCosh(const double _Lt): Lt(_Lt){}
  
  //Params are A, m  
  ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const{
    return p.A * ( exp(-p.m*t) + exp(-p.m*(double(Lt)-t)) );
  }
  void parameterDerivatives(ValueDerivativeType &yderivs, const GeneralizedCoordinate &t, const ParameterType &p) const{
    yderivs.dA = exp(-p.m*t) + exp(-p.m*(Lt-t));
    yderivs.dm = p.A * ( -t*exp(-p.m*t) + -(Lt-t)*exp(-p.m*(Lt-t)) );
  }

  inline int Nparams() const{ return 2; }
};

enum DataType { PP_LW_data, AA_LW_data };
std::string toStr(const DataType d){
  const static std::vector<std::string> str = {"PP_LW_data", "AA_LW_data"};
  int dd(d);
  if(dd < 0 || dd >= str.size()){
    std::ostringstream os; os << "Unknown DataType idx " << dd;
    return os.str();
  }else return str[int(d)];
}

struct AllFitParams{
  double m;
  double A_PP_LW;
  double A_AA_LW;

  typedef AllFitParams ET_tag;
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,AllFitParams>::value, int>::type = 0>
  AllFitParams(U&& expr){
    m = expr[0];
    A_PP_LW = expr[1];
    A_AA_LW = expr[2];
  }
  AllFitParams() = default;
  AllFitParams &operator=(const AllFitParams &r) = default;
  AllFitParams &operator=(const double r){ m=r; A_PP_LW=r; A_AA_LW=r; }
  inline double & operator()(const int i){
    switch(i){
    case 0:
      return m;
    case 1:
      return A_PP_LW;
    case 2:
      return A_AA_LW;
    default:
      assert(0);
    };
  }
  inline const double & operator()(const int i) const{ return const_cast<const double &>(const_cast<AllFitParams*>(this)->operator()(i)); }
  
  static inline const int size(){ return 3; }
  inline void zero(){ m=A_PP_LW=A_AA_LW=0.; }
  
  std::string print() const;
};
#define STRUCT_TYPE AllFitParams
#define STRUCT_ARGS ELEM(double, m) ELEM(double, A_PP_LW) ELEM(double, A_AA_LW)
#include<parser_gen.incl>

std::string AllFitParams::print() const{ std::ostringstream os; os << *this; return os.str(); }

template<>
struct getElem<AllFitParams>{
  static inline double elem(const AllFitParams &v, const int i){
    return v(i);
  }    
  inline static int common_properties(const AllFitParams &v){ return 0; }
};



struct AllFitParamDerivs{
  double dm;
  double dA_PP_LW;
  double dA_AA_LW;
  
  inline AllFitParamDerivs(): dm(0.), dA_PP_LW(0.), dA_AA_LW(0.){}

  inline const double &operator()(const int i) const{
    switch(i){
    case 0:
      return dm;
    case 1:
      return dA_PP_LW;
    case 2:
      return dA_AA_LW;      
    default:
      assert(0);
    };
  }
    
  inline static size_t size(){ return 3; }
  
  inline void import(const FitParamDerivs &p, const DataType type){    
    dm = p.dm;
    switch(type){
    case PP_LW_data:
      dA_PP_LW = p.dA; break;
    case AA_LW_data:
      dA_AA_LW = p.dA; break;      
    default:
      std::cout << "Error: attempting to evaluate AllFitParamDerivs::import with data type " << toStr(type) << std::endl;
      exit(-1);
    }
  }
};


//Coordinate, parameters and param derivatives for aggregate fit func
struct Coord{
  double t; //time
  DataType type;
  Coord(){}
  Coord(const double _t, const DataType _type): t(_t), type(_type){}
};
std::ostream & operator<<(std::ostream &os, const Coord &c){
  os << "(" << c.t << "," << toStr(c.type) << ")";
  return os;
}

class FitMpi{
  FitCosh fcosh;  
public:
  typedef double ValueType;
  typedef AllFitParams ParameterType;
  typedef AllFitParamDerivs ValueDerivativeType;
  typedef Coord GeneralizedCoordinate;

  FitMpi(const double _Lt): fcosh(_Lt){}
  
  ValueType value(const GeneralizedCoordinate &c, const ParameterType &p) const{
    switch(c.type){
    case PP_LW_data:
      return fcosh.value(c.t, FitParams(p.A_PP_LW, p.m));
    case AA_LW_data:
      return fcosh.value(c.t, FitParams(p.A_AA_LW, p.m));      
    default:
      std::cout << "Error: attempting to evaluate FitMpi::value with coordinate " << c << std::endl;
      exit(-1);
    };    
  }
  void parameterDerivatives(ValueDerivativeType &yderivs, const GeneralizedCoordinate &c, const ParameterType &p) const{
    FitParamDerivs subderivs;
    switch(c.type){
    case PP_LW_data:
      fcosh.parameterDerivatives(subderivs, c.t, FitParams(p.A_PP_LW, p.m)); break;
    case AA_LW_data:
      fcosh.parameterDerivatives(subderivs, c.t, FitParams(p.A_AA_LW, p.m)); break;      
    default:
      std::cout << "Error: attempting to evaluate FitMpi::parameterDerivatives with coordinate " << c << std::endl;
      exit(-1);
    };    
    yderivs.import(subderivs, c.type);
  }

  inline int Nparams() const{ return AllFitParams::size(); }
};


class FitMpiFrozen{
  FitMpi ff;
  AllFitParams base_params;
  std::vector<bool> do_freeze;
  int n_frozen;
  
public:
  typedef double ValueType;
  typedef NumericVector<double> ParameterType;
  typedef NumericVector<double> ValueDerivativeType;
  typedef Coord GeneralizedCoordinate;

  FitMpiFrozen(const double _Lt): ff(_Lt), do_freeze(AllFitParams::size(),false), n_frozen(0){}

  inline ValueType value(const GeneralizedCoordinate &c, const ParameterType &p) const{
    AllFitParams pp; expand(pp,p);
    return ff.value(c,pp);
  }
  inline void parameterDerivatives(ValueDerivativeType &yderivs, const GeneralizedCoordinate &c, const ParameterType &p) const{
    AllFitParams pp; expand(pp,p);
    AllFitParamDerivs ppd;
    ff.parameterDerivatives(ppd, c, pp);
    reduce(yderivs, ppd);
  }
  inline int Nparams() const{ return AllFitParams::size() - n_frozen; }

  inline void freeze(const int pidx, const double value){
    do_freeze[pidx] = true;
    base_params(pidx) = value;
    ++n_frozen;
  }
  inline void expand(AllFitParams &o, const NumericVector<double> &r) const{
    o = base_params;
    int j=0;
    for(int i=0;i<o.size();i++) if(!do_freeze[i]) o(i) = r(j++);
    assert(j==Nparams());
  }
  inline void reduce(NumericVector<double> &r, const AllFitParams &o) const{
    r.resize(Nparams());
    int j=0;
    for(int i=0;i<o.size();i++) if(!do_freeze[i]) r(j++) = o(i);
    assert(j==Nparams());
  }  
  inline void reduce(NumericVector<double> &r, const AllFitParamDerivs &o) const{
    r.resize(Nparams());
    int j=0;
    for(int i=0;i<o.size();i++) if(!do_freeze[i]) r(j++) = o(i);
    assert(j==Nparams());
  }
  
};


#endif
