#ifndef FIT_SIMULTANEOUS_FITFUNC_H_
#define FIT_SIMULTANEOUS_FITFUNC_H_

struct FitSimCoord{
  int idx;
  double t;
  FitSimCoord(const int idx, const double t): idx(idx),t(t){}
};

typedef correlationFunction<FitSimCoord,jackknifeDistributionD> jackknifeSimFitCorrelationFunction;
typedef correlationFunction<FitSimCoord,doubleJackknifeDistributionD> doubleJackknifeSimFitCorrelationFunction;


struct FitSimParams{
  double m;
  std::vector<double> A;

  explicit FitSimParams(const int n): A(n-1,0.), m(0.){}
  FitSimParams() = default;
  FitSimParams(const FitSimParams &r) = default;
  FitSimParams(FitSimParams &&r) = default;
    
  FitSimParams & operator=(const FitSimParams &r) = default;
  FitSimParams & operator=(FitSimParams &&r) = default;

  ENABLE_GENERIC_ET(FitSimParams, FitSimParams, FitSimParams);

  inline StandardFitParams getParams(const int i) const{
    return StandardFitParams(A[i],m);
  }
  inline void setParams(const int i, const StandardFitParams &in){
    m = in.m;
    A[i] = in.A;
  }
  
  void resize(const int n){ A.resize(n-1); }
  inline int size() const{ return 1 + A.size(); }
  
  inline double & operator()(const int i){ return i==0 ? m : A[i-1]; }
  inline double operator()(const int i) const{ return i==0 ? m : A[i-1]; }
  inline double & operator[](const int i){ return i==0 ? m : A[i-1]; }
  inline double operator[](const int i) const{ return i==0 ? m : A[i-1]; }

  std::string print() const{
    std::ostringstream os; 
    os << "(m=" << m;
    for(int i=0;i<A.size();i++) os << " A" << i << "=" << A[i];
    os << ")";
    return os.str();
  }

  void zero(){ m=0.; for(int i=0;i<A.size();i++) A[i] =0.; }

  GENERATE_HDF5_SERIALIZE_METHOD( (m)(A) );
};

inline std::ostream & operator<<(std::ostream &os, const FitSimParams &p){
  os << p.print(); return os;
}

GENERATE_HDF5_SERIALIZE_FUNC(FitSimParams);

class FitFuncSimultaneous{
public:
  typedef double ValueType;
  typedef FitSimParams ParameterType;
  typedef FitSimParams ValueDerivativeType; //derivative wrt parameters
  typedef FitSimCoord GeneralizedCoordinate; //time coord

private:
  std::vector<StandardFitFuncBase const*> fitfuncs;

public:
  void setNcorrelators(const int n){ fitfuncs.resize(n,NULL); }
  void setCorrelatorFitFunc(const int idx, StandardFitFuncBase const* to){
    assert(fitfuncs.size() > idx);
    fitfuncs[idx] = to;
  }

  ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const{
    StandardFitParams sfp = p.getParams(t.idx);
    return fitfuncs[t.idx]->value(t.t, sfp);
  }

  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    StandardFitParams sfp = p.getParams(t.idx);
    StandardFitParams derivs = fitfuncs[t.idx]->parameterDerivatives(t.t, sfp);
    ValueDerivativeType out(1+fitfuncs.size());
    out.setParams(t.idx, derivs);
    return out;
  }

  int Nparams() const{ return fitfuncs.size() + 1; }

  FitSimParams guess() const{ 
    FitSimParams out(1+fitfuncs.size());
    out.m=0.5;
    for(int i=0;i<fitfuncs.size();i++) out.A[i] = 1e6;
    return out;
  }
};

#endif
