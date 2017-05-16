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
class FitSinh{
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
  void parameterDerivatives(ValueDerivativeType &yderivs, const GeneralizedCoordinate &t, const ParameterType &p) const{
    yderivs.dA = exp(-p.m*t) - exp(-p.m*(Lt-t));
    yderivs.dm = p.A * ( -t*exp(-p.m*t) - -(Lt-t)*exp(-p.m*(Lt-t)) );
  }

  inline int Nparams() const{ return 2; }
};




//Define the amplitudes, their enumerations (starting from 1) and their associated fit forms
#define ALL_FIT_PARAMS_AMPL_TYPES \
  (1,PP_LW,fcosh)\
  (2,AA_LW,fcosh)\
  (3,AP_LW,fsinh)\
  (4,PP_WW,fcosh)\
  (5,AP_WW,fsinh)

#define _AFP_IDX(ENUM_IDX,AMP,FORM) ENUM_IDX
#define _AFP_AMP(ENUM_IDX,AMP,FORM) AMP
#define _AFP_FORM(ENUM_IDX,AMP,FORM) FORM


//Define the enum datatypes
#define _SETUP_DATATYPE_ELEM_E(r,data,i,elem) BOOST_PP_COMMA_IF(i) BOOST_PP_CAT(_AFP_AMP elem,_data)
#define _SETUP_DATATYPE_ELEMSTR_E(r,data,i,elem) BOOST_PP_COMMA_IF(i) BOOST_PP_STRINGIZE(BOOST_PP_CAT(_AFP_AMP elem,_data))
#define _SETUP_DATATYPE_ELEM_DERIVMAP_E(r,data,elem) (BOOST_PP_CAT(_AFP_AMP elem,_data), BOOST_PP_CAT(dA_,_AFP_AMP elem) )
#define _SETUP_DATATYPE_ELEM_AMPFORMMAP_E(r,data,elem) (BOOST_PP_CAT(_AFP_AMP elem,_data), BOOST_PP_CAT(A_,_AFP_AMP elem), _AFP_FORM elem)

#define ALL_FIT_PARAMS_DATATYPES TUPLE_SEQUENCE_FOR_EACH_I(_SETUP_DATATYPE_ELEM_E, , ALL_FIT_PARAMS_AMPL_TYPES)
#define ALL_FIT_PARAMS_DATATYPES_STR TUPLE_SEQUENCE_FOR_EACH_I(_SETUP_DATATYPE_ELEMSTR_E, , ALL_FIT_PARAMS_AMPL_TYPES)
#define ALL_FIT_PARAMS_DATATYPE_DERIV_MAP TUPLE_SEQUENCE_FOR_EACH(_SETUP_DATATYPE_ELEM_DERIVMAP_E, , ALL_FIT_PARAMS_AMPL_TYPES)
#define ALL_FIT_PARAMS_DATATYPE_AMP_FORM_MAP TUPLE_SEQUENCE_FOR_EACH(_SETUP_DATATYPE_ELEM_AMPFORMMAP_E, , ALL_FIT_PARAMS_AMPL_TYPES)

enum DataType { ALL_FIT_PARAMS_DATATYPES, NDataType };
std::string toStr(const DataType d){
  const static std::vector<std::string> str = {ALL_FIT_PARAMS_DATATYPES_STR, "NDataType"};
  int dd(d);
  if(dd < 0 || dd >= str.size()){
    std::ostringstream os; os << "Unknown DataType idx " << dd;
    return os.str();
  }else return str[int(d)];
}

//Create the preprocessor tuple sequences containing the information needed to generate the parameter and derivative classes
#define _SETUP_MEMBER_ELEM(I,AMP,FORM) (I,double,BOOST_PP_CAT(A_,AMP) )
#define _SETUP_MEMBER_ELEM_E(r,data,elem) _SETUP_MEMBER_ELEM elem

//    List of (enum, type, name) for all data members of AllFitParams
#define ALL_FIT_PARAMS_ENUM_MEMBERS \
  (0,double,m) \
  TUPLE_SEQUENCE_FOR_EACH(_SETUP_MEMBER_ELEM_E, , ALL_FIT_PARAMS_AMPL_TYPES)

#define _EXTRACT23(A,B,C) (B,C)
#define _EXTRACT23_E(r,data,elem) _EXTRACT23 elem
#define _EXTRACT13(A,B,C) (A,C)
#define _EXTRACT13_E(r,data,elem) _EXTRACT13 elem

//    List of (enum, name) and (type, name) for all data members of AllFitParams
#define ALL_FIT_PARAMS_MEMBERS TUPLE_SEQUENCE_FOR_EACH(_EXTRACT23_E, , ALL_FIT_PARAMS_ENUM_MEMBERS)
#define ALL_FIT_PARAMS_ENUMS TUPLE_SEQUENCE_FOR_EACH(_EXTRACT13_E, , ALL_FIT_PARAMS_ENUM_MEMBERS)

//    Macros for generating function bodies
#define _ALLFITPARAMS_MEM_EXPR(r,data,elem) BOOST_PP_TUPLE_ELEM(1,elem) = expr[BOOST_PP_TUPLE_ELEM(0,elem)];
#define _ALLFITPARAMS_MEM_ACCESSOR(r,data,elem) case BOOST_PP_TUPLE_ELEM(0,elem): return BOOST_PP_TUPLE_ELEM(1,elem);
#define _ALLFITPARAMS_MEM_SET_R(rpt,data,elem) BOOST_PP_TUPLE_ELEM(1,elem) = r;
#define _ALLFITPARAMS_GET_INDEX(rpt,str,elem) if(str == BOOST_PP_STRINGIZE(BOOST_PP_TUPLE_ELEM(1,elem))) return BOOST_PP_TUPLE_ELEM(0,elem);
  
//    The parameter class
struct AllFitParams{
  GENERATE_MEMBERS(ALL_FIT_PARAMS_MEMBERS)

  typedef AllFitParams ET_tag;
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,AllFitParams>::value, int>::type = 0>
  AllFitParams(U&& expr){
    TUPLE_SEQUENCE_FOR_EACH(_ALLFITPARAMS_MEM_EXPR,,ALL_FIT_PARAMS_ENUMS)
  }
  AllFitParams() = default;
  AllFitParams &operator=(const AllFitParams &r) = default;
  AllFitParams &operator=(const double r){
    TUPLE_SEQUENCE_FOR_EACH(_ALLFITPARAMS_MEM_SET_R,,ALL_FIT_PARAMS_MEMBERS)
  }
  inline double & operator()(const int i){
    switch(i){
      TUPLE_SEQUENCE_FOR_EACH(_ALLFITPARAMS_MEM_ACCESSOR,,ALL_FIT_PARAMS_ENUMS)
    default:
      assert(0);
    };
  }
  inline const double & operator()(const int i) const{ return const_cast<const double &>(const_cast<AllFitParams*>(this)->operator()(i)); }
  
  static inline const int size(){ return TUPLE_SEQUENCE_SIZE(ALL_FIT_PARAMS_ENUM_MEMBERS); }
  inline void zero(){ (*this) = 0.; }
  
  std::string print() const;

  static inline int getIndex(const std::string param){
    TUPLE_SEQUENCE_FOR_EACH(_ALLFITPARAMS_GET_INDEX,param,ALL_FIT_PARAMS_ENUMS)
    error_exit(std::cout << "Error AllFitParams::getIndex  unknown parameter string \"" << param << "\"\n");
  }
};

GENERATE_PARSER(AllFitParams,ALL_FIT_PARAMS_MEMBERS);

std::string AllFitParams::print() const{ std::ostringstream os; os << *this; return os.str(); }

template<>
struct getElem<AllFitParams>{
  static inline double elem(const AllFitParams &v, const int i){
    return v(i);
  }    
  inline static int common_properties(const AllFitParams &v){ return 0; }
};

//   List of (enum, name) and (type, name) for all data members of AllFitParamsDerivs
#define _ADD_D(A,B) (A, BOOST_PP_CAT(d,B))
#define _ADD_D_ELEM(r,data,elem) _ADD_D elem
#define ALL_FIT_PARAMS_DMEMBERS TUPLE_SEQUENCE_FOR_EACH(_ADD_D_ELEM, , ALL_FIT_PARAMS_MEMBERS)
#define ALL_FIT_PARAMS_DENUMS TUPLE_SEQUENCE_FOR_EACH(_ADD_D_ELEM, , ALL_FIT_PARAMS_ENUMS)

#define ALL_FIT_PARAMS_ELEM_DCASE_MAP(TYPE, AMP) case TYPE: AMP = p.dA; break;
#define ALL_FIT_PARAMS_ELEM_DCASE_MAP_E(r,data,elem) ALL_FIT_PARAMS_ELEM_DCASE_MAP elem

//   The derivative class
struct AllFitParamDerivs{
  GENERATE_MEMBERS(ALL_FIT_PARAMS_DMEMBERS)
  
  inline AllFitParamDerivs() = default;

  inline const double &operator()(const int i) const{
    switch(i){
    TUPLE_SEQUENCE_FOR_EACH(_ALLFITPARAMS_MEM_ACCESSOR,,ALL_FIT_PARAMS_DENUMS)   
    default:
      assert(0);
    };
  }
    
  inline static size_t size(){ return TUPLE_SEQUENCE_SIZE(ALL_FIT_PARAMS_ENUM_MEMBERS); }
  
  inline void import(const FitParamDerivs &p, const DataType type){    
    dm = p.dm;
    switch(type){
    TUPLE_SEQUENCE_FOR_EACH(ALL_FIT_PARAMS_ELEM_DCASE_MAP_E,,ALL_FIT_PARAMS_DATATYPE_DERIV_MAP)
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

#define ALL_FIT_PARAMS_ELEM_VALUE_CASE_MAP(TYPE, AMP, FORM) case TYPE: return FORM.value(c.t, FitParams(p. AMP, p.m));
#define ALL_FIT_PARAMS_ELEM_VALUE_CASE_MAP_E(r,data,elem) ALL_FIT_PARAMS_ELEM_VALUE_CASE_MAP elem

#define ALL_FIT_PARAMS_ELEM_DERIVS_CASE_MAP(TYPE, AMP, FORM) case TYPE: FORM.parameterDerivatives(subderivs, c.t, FitParams(p. AMP, p.m)); break;
#define ALL_FIT_PARAMS_ELEM_DERIVS_CASE_MAP_E(r,data,elem) ALL_FIT_PARAMS_ELEM_DERIVS_CASE_MAP elem

class FitMpi{
  FitCosh fcosh;
  FitSinh fsinh;
public:
  typedef double ValueType;
  typedef AllFitParams ParameterType;
  typedef AllFitParamDerivs ValueDerivativeType;
  typedef Coord GeneralizedCoordinate;

  FitMpi(const double _Lt): fcosh(_Lt), fsinh(_Lt){}
  
  ValueType value(const GeneralizedCoordinate &c, const ParameterType &p) const{
    switch(c.type){
    TUPLE_SEQUENCE_FOR_EACH(ALL_FIT_PARAMS_ELEM_VALUE_CASE_MAP_E,,ALL_FIT_PARAMS_DATATYPE_AMP_FORM_MAP)  
    default:
      std::cout << "Error: attempting to evaluate FitMpi::value with coordinate " << c << std::endl;
      exit(-1);
    };    
  }
  void parameterDerivatives(ValueDerivativeType &yderivs, const GeneralizedCoordinate &c, const ParameterType &p) const{
    FitParamDerivs subderivs;
    switch(c.type){
    TUPLE_SEQUENCE_FOR_EACH(ALL_FIT_PARAMS_ELEM_DERIVS_CASE_MAP_E,,ALL_FIT_PARAMS_DATATYPE_AMP_FORM_MAP)   
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
  inline void freeze(const std::string &paramname, const double value){
    freeze(AllFitParams::getIndex(paramname),value);
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


#define FIT_EFFMASS_ELEM_VALUE_CASE_MAP(TYPE, AMP, FORM) case TYPE: return FORM.value(t,p)/FORM.value(t+1,p);
#define FIT_EFFMASS_ELEM_VALUE_CASE_MAP_E(r,data,elem) FIT_EFFMASS_ELEM_VALUE_CASE_MAP elem

#define FIT_EFFMASS_ELEM_DERIVS_CASE_MAP(TYPE, AMP, FORM) \
  case TYPE: \
  FORM.parameterDerivatives(subderivs_t, t,p); \
  value_t = FORM.value(t,p); \
  FORM.parameterDerivatives(subderivs_tp1, t+1,p); \
  value_tp1 = FORM.value(t+1,p); \
  break;
#define FIT_EFFMASS_ELEM_DERIVS_CASE_MAP_E(r,data,elem) FIT_EFFMASS_ELEM_DERIVS_CASE_MAP elem

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
public:
  typedef double ValueType;
  typedef MLwrapper<double> ParameterType;
  typedef MLwrapper<double> ValueDerivativeType;
  typedef double GeneralizedCoordinate;

  FitEffectiveMass(const double _Lt, const DataType _type): fcosh(_Lt), fsinh(_Lt), type(_type){}

  double value(const double t, const ParameterType &params) const{    
    FitParams p(1000, *params);
    switch(type){
      TUPLE_SEQUENCE_FOR_EACH(FIT_EFFMASS_ELEM_VALUE_CASE_MAP_E,,ALL_FIT_PARAMS_DATATYPE_AMP_FORM_MAP)  
    default:
      error_exit(std::cout << "FitEffectiveMass::value unknown type\n");
    }
  }
  
  void parameterDerivatives(ValueDerivativeType &yderivs, const double t, const ParameterType &params) const{
    FitParams p(1000, *params);
    FitParamDerivs subderivs_t;
    FitParamDerivs subderivs_tp1;
    double value_t, value_tp1;
    switch(type){
    TUPLE_SEQUENCE_FOR_EACH(FIT_EFFMASS_ELEM_DERIVS_CASE_MAP_E,,ALL_FIT_PARAMS_DATATYPE_AMP_FORM_MAP)   
    default:
      error_exit(std::cout << "FitEffectiveMass::parameterDerivatives unknown type\n");
    };    
    *yderivs = subderivs_t.dm/value_tp1 - value_t/value_tp1/value_tp1 * subderivs_tp1.dm;
  }

  inline int Nparams() const{ return 1; }
  
  
};





#endif
