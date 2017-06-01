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

struct FitFuncBase{
  typedef double ValueType;
  typedef FitParams ParameterType;
  typedef FitParamDerivs ValueDerivativeType; //derivative wrt parameters
  typedef double GeneralizedCoordinate; //time coord

  virtual ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const = 0;
  virtual ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const = 0;
  virtual int Nparams() const = 0;
};


class FitCosh: public FitFuncBase{
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
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    ValueDerivativeType yderivs;
    yderivs.dA = exp(-p.m*t) + exp(-p.m*(Lt-t));
    yderivs.dm = p.A * ( -t*exp(-p.m*t) + -(Lt-t)*exp(-p.m*(Lt-t)) );
    return yderivs;
  }

  inline int Nparams() const{ return 2; }
};
class FitSinh: public FitFuncBase{
  const double Lt;
public:
  typedef double ValueType;
  typedef FitParams ParameterType;
  typedef FitParamDerivs ValueDerivativeType; //derivative wrt parameters
  typedef double GeneralizedCoordinate; //time coord

  FitSinh(const double _Lt): Lt(_Lt){}
  
  //Params are A, m  
  ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const{
    return p.A * ( exp(-p.m*t) - exp(-p.m*(double(Lt)-t)) );
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    ValueDerivativeType yderivs;
    yderivs.dA = exp(-p.m*t) - exp(-p.m*(Lt-t));
    yderivs.dm = p.A * ( -t*exp(-p.m*t) - -(Lt-t)*exp(-p.m*(Lt-t)) );
    return yderivs;
  }

  inline int Nparams() const{ return 2; }
};



//Define the amplitudes, their enumerations (starting from 1) and their associated fit forms
GENERATE_ENUM_AND_PARSER(Params, (Mass)(A_PP_LW)(A_AA_LW)(A_AP_LW)(A_PP_WW)(A_AP_WW)(Params_size) )

struct FitControls{

  //Get the amplitude enum associated with a particular data type
  static inline Params amplitude(const DataType type){
    switch(type){
    case PP_LW_data:
      return A_PP_LW;
    case AA_LW_data:
      return A_AA_LW;
    case AP_LW_data:
      return A_AP_LW;
    case PP_WW_data:
      return A_PP_WW;
    case AP_WW_data:
      return A_AP_WW;
    default:
      error_exit(std::cout << "FitControls::amplitude invalid type " << type << std::endl);
    };
  }
  static inline FitFuncBase const* baseFitFunc(const DataType type, FitFuncBase const* fcosh, FitFuncBase const* fsinh){
    switch(type){
    case PP_LW_data:
    case AA_LW_data:
    case PP_WW_data:
      return fcosh;
    case AP_LW_data:
    case AP_WW_data:
      return fsinh;
    default:
      error_exit(std::cout << "FitControls::baseFitFunc invalid type\n");
    }
  }

};





class AllParamMap{
public:

  typedef Params tagType;
private:
  std::vector<int> pmap; //(int)Params -> index
  std::vector<Params> upmap; //index -> Params

public:
  AllParamMap(const Args &args): pmap((int)Params_size), upmap(1){
    pmap[(int)Mass] = 0;
    upmap[0] = Mass;

    int idx = 1;
    for(int i=0;i<args.data.size();i++){
      if(args.data[i].FF_data.include_data || args.data[i].BB_data.include_data){      
	Params amp = FitControls::amplitude(args.data[i].type);
	pmap[(int)amp] = idx++;
	upmap.push_back(amp);
      }
    }
  }
  inline int map(const Params param) const{
    return pmap[(int)param];
  }
  inline Params unmap(const int idx) const{
    return upmap[idx];
  }
  inline int size() const{ return upmap.size(); }
  
  inline std::string paramName(const int idx) const{
    return toString(upmap[idx]);
  }
  inline int nParams() const{ return upmap.size(); }
};



//Coordinate, parameters and param derivatives for aggregate fit func
struct Coord{
  double t; //time
  DataType type;
  Coord(){}
  Coord(const double _t, const DataType _type): t(_t), type(_type){}
};
std::ostream & operator<<(std::ostream &os, const Coord &c){
  os << "(" << c.t << "," << c.type << ")";
  return os;
}

class FitMpi{
public:
  typedef double ValueType;
  typedef mappedVector<double, AllParamMap> ParameterType;
  typedef mappedVector<double, AllParamMap> ValueDerivativeType;
  typedef Coord GeneralizedCoordinate;
private:
  FitCosh fcosh;
  FitSinh fsinh;
    
  inline FitParams reduce(const ParameterType &p, const DataType type) const{
    const double &A = p(FitControls::amplitude(type));
    const double &m = p(Mass);
    return FitParams(A,m);
  }
  int nparams;
public:

  FitMpi(const double _Lt, const AllParamMap &param_map): fcosh(_Lt), fsinh(_Lt), nparams(param_map.nParams()){  }
  
  inline ValueType value(const GeneralizedCoordinate &c, const ParameterType &p) const{
    return FitControls::baseFitFunc(c.type, &fcosh, &fsinh)->value(c.t,reduce(p,c.type));
  }

  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &c, const ParameterType &p) const{
    ValueDerivativeType yderivs(p.getMapping());
    yderivs.zero();
    FitParamDerivs subderivs = FitControls::baseFitFunc(c.type, &fcosh, &fsinh)->parameterDerivatives(c.t,reduce(p,c.type));
    yderivs(Mass) = subderivs.dm;
    yderivs(FitControls::amplitude(c.type)) = subderivs.dA;
    return yderivs;
  }

  inline int Nparams() const{ return nparams; }
};



template<typename T>
class MLwrapper{
  T t;
public:
  ENABLE_GENERIC_ET(MLwrapper, MLwrapper<T>);
  //typedef MLwrapper<T> ET_tag;
  MLwrapper(): t(0.){}
  explicit MLwrapper(const T _t): t(_t){}
  inline T & operator()(const int i){ assert(i==0); return t; }
  inline const T &operator()(const int i) const{ assert(i==0); return t; }
  inline int size() const{ return 1; }
  inline std::string print() const{ return anyToStr<T>(t); }
  inline void resize(const int i){ assert(i==1); }
  inline void zero(){ t=0.; }
  inline T& operator*(){ return t; }
  inline const T& operator*() const{ return t; }
};
template<typename T>
struct getElem<MLwrapper<T> >{
  static inline T& elem(MLwrapper<T> &v, const int i){ return *v; }
  static inline const T& elem(const MLwrapper<T> &v, const int i){ return *v; }    
  inline static int common_properties(const MLwrapper<T> &v){ return 0; }
};




class FitEffectiveMass{
  FitCosh fcosh;
  FitSinh fsinh;
  DataType type;
  FitFuncBase const* fitfunc;
public:
  typedef double ValueType;
  typedef MLwrapper<double> ParameterType;
  typedef MLwrapper<double> ValueDerivativeType;
  typedef double GeneralizedCoordinate;

  FitEffectiveMass(const double _Lt, const DataType _type): fcosh(_Lt), fsinh(_Lt), type(_type), fitfunc(FitControls::baseFitFunc(_type, &fcosh, &fsinh)){}

  inline double value(const double t, const ParameterType &params) const{    
    FitParams p(1000, *params);
    return fitfunc->value(t,p)/fitfunc->value(t+1,p);
  }
  
  ValueDerivativeType parameterDerivatives(const double t, const ParameterType &params) const{
    ValueDerivativeType yderivs;
    FitParams p(1000, *params);
    
    double value_t = fitfunc->value(t,p);
    FitParamDerivs subderivs_t = fitfunc->parameterDerivatives(t,p);

    double value_tp1 = fitfunc->value(t+1,p);
    FitParamDerivs subderivs_tp1 = fitfunc->parameterDerivatives(t+1,p);
 
    *yderivs = subderivs_t.dm/value_tp1 - value_t/value_tp1/value_tp1 * subderivs_tp1.dm;
    return yderivs;
  }

  inline int Nparams() const{ return 1; }
  
  
};





#endif
