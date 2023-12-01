#ifndef _GLOBALFIT_PARAMS_H_
#define _GLOBALFIT_PARAMS_H_

//BaseParams has the set of named parameters, this class adds the scaling parameters and wraps the bindings
template<typename BaseParams>
class GlobalFitParams{
public:
  std::vector<double> Zl;
  std::vector<double> Zh;
  std::vector<double> Ra;
  BaseParams params;
private:

  std::vector<double*> pmap;
  void gen_pmap(){
    int nlat = Zl.size();
    pmap.resize(3*nlat + params.size());
    for(int i=0;i<nlat;i++){
      pmap[i] = &Zl[i];
      pmap[i+nlat] = &Zh[i];
      pmap[i+2*nlat] = &Ra[i];
    }
    for(int i=0;i<params.size();i++)
      pmap[3*nlat + i] = &params(i);
  }
public:

  struct _cprops: public OstreamHook{ 
    int nlat; 
    _cprops(int nlat): nlat(nlat){} 
    inline bool operator!=(const _cprops &r) const{ return nlat != r.nlat; }
    void write(std::ostream &os) const{ os << "(nlat = " << nlat << ")"; }
  };
  inline _cprops cprops() const{ return _cprops(Zl.size()); }  

  GlobalFitParams() = default;

  GlobalFitParams(const int nlat): Zl(nlat,1.), Zh(nlat,1.), Ra(nlat,1.){
    gen_pmap();
  }
  GlobalFitParams(const _cprops &props): GlobalFitParams(props.nlat){}  

  GlobalFitParams(const GlobalFitParams &r): Zl(r.Zl), Zh(r.Zh), Ra(r.Ra), params(r.params){
    gen_pmap();
  }
  GlobalFitParams(GlobalFitParams &&r): Zl(std::move(r.Zl)), Zh(std::move(r.Zh)), Ra(std::move(r.Ra)), params(std::move(r.params)){
    gen_pmap();
  }

  GlobalFitParams &operator=(const GlobalFitParams &r){
    Zl = r.Zl; Zh = r.Zh; Ra = r.Ra; params = r.params; gen_pmap(); return *this;
  }
  GlobalFitParams &operator=(GlobalFitParams &&r){
    Zl = std::move(r.Zl); Zh = std::move(r.Zh); Ra = std::move(r.Ra); params = std::move(r.params); gen_pmap(); return *this;
  }


  template<typename Binding>
  static void supersetToSubset(typename Binding::subsetType &psub, const GlobalFitParams<BaseParams> &psup, const DataParams &t){
    psub.Zl = psup.Zl[t.lattice];
    psub.Zh = psup.Zh[t.lattice];
    psub.Ra = psup.Ra[t.lattice];
    Binding::supersetToSubset(psub,psup.params);
  }

  template<typename Binding>
  static void subsetToSuperset(GlobalFitParams<BaseParams> &psup, const typename Binding::subsetType &psub, const DataParams &t){
    psup.Zl[t.lattice] = psub.Zl;
    psup.Zh[t.lattice] = psub.Zh;
    psup.Ra[t.lattice] = psub.Ra;
    Binding::subsetToSuperset(psup.params,psub);
  }

  inline int size() const{ return 3*Zl.size() + params.size(); }

  inline void resize(const _cprops &props){ Zl.resize(props.nlat); Zh.resize(props.nlat); Ra.resize(props.nlat); }

  inline double &operator()(const int i){
    return *pmap[i];
  }
  inline const double &operator()(const int i) const{
    return *pmap[i];
  }

  inline void zero(){ for(int i=0;i<size();i++) this->operator()(i) = 0.; }

  //O(N) search for member with pointer 'mem'. Useful to lookup indices of particular members
  inline int paramIdx(double const* mem){
    for(int i=0;i<pmap.size();i++) if(pmap[i] == mem) return i;
    return -1;
  }

  ENABLE_GENERIC_ET(GlobalFitParams, GlobalFitParams<BaseParams>, GlobalFitParams<BaseParams> );
};

template<typename T>
struct SARLaC::getElem< GlobalFitParams<T> >{
  static inline double elem(const GlobalFitParams<T> &v, const int i){ return v(i); }
  static inline double& elem(GlobalFitParams<T> &v, const int i){ return v(i); }
  static inline typename GlobalFitParams<T>::_cprops common_properties(const GlobalFitParams<T> &v){ return v.cprops(); }
};


#endif
