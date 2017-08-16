#ifndef _FIT_KTOPIPI_MAIN_H
#define _FIT_KTOPIPI_MAIN_H


//Functor for tensor reduction by average
template<typename T, int Rank>
struct averageDimensionFunctor{
  const int dim;
  std::vector<int> const* use;
  
  averageDimensionFunctor(const int _dim, std::vector<int> const*_use = NULL): dim(_dim), use(_use){}
  
  void operator()(T &o, int const *coord, const NumericTensor<T,Rank> &from) const{
    int full_coord[Rank];
    int i=0; for(int ii=0;ii<Rank;ii++) if(ii!=dim) full_coord[ii] = coord[i++];    
    zeroit(o);
    if(use != NULL){
      assert(use->size()> 0);
      full_coord[dim] = use->at(0);
      o = from(full_coord);      
      for(int i=1;i<use->size();i++){
	full_coord[dim] = use->at(i);
	o = o + from(full_coord);
      }
      o = o/double(use->size());
    }else{
      assert(from.size(dim)>0);
      full_coord[dim] = 0;
      o = from(full_coord);
      for(int i=1;i<from.size(dim);i++){
	full_coord[dim] = i;
	o = o + from(full_coord);
      }
      o = o/double(from.size(dim));
    }
  }
};
    
template<typename Resampled, typename Raw>
struct resampleFunctor{
  inline Resampled operator()(int const* coord,const Raw &from) const{
    Resampled o(from.size());
    o.resample(from);
    return o;
  }
};

template<typename distributionType>
struct iterate;

template<typename T>
struct iterate<doubleJackknifeDistribution<T> >{
  static inline int size(const doubleJackknifeDistribution<T> &from){ return from.size() * (from.size()-1); } //j + (from.size()-1)*i
  static inline const T& at(const int i, const doubleJackknifeDistribution<T> &from){
    const int nn = from.size()-1;
    return from.sample(i/nn).sample(i%nn);
  }
  static inline T & at(const int i, doubleJackknifeDistribution<T> &from){
    const int nn = from.size()-1;
    return from.sample(i/nn).sample(i%nn);
  }    
};
template<typename T>
struct iterate<rawDataDistribution<T> >{
  static inline int size(const rawDataDistribution<T> &from){ return from.size(); } 
  static inline const T& at(const int i, const rawDataDistribution<T> &from){
    return from.sample(i);
  }
  static inline T & at(const int i, rawDataDistribution<T> &from){
    return from.sample(i);
  }    
};

  

template<typename distributionType>
struct multiplyBubbleFunctor{
  const NumericTensor<distributionType,1> &bubble;
  int tsep_k_pi;
  int tsep_pipi;
  int Lt;
  int tK_dim;
  multiplyBubbleFunctor(const NumericTensor<distributionType,1> &_bubble, const int _tsep_k_pi, const int _tsep_pipi, const int _Lt, const int _tK_dim): bubble(_bubble), tsep_k_pi(_tsep_k_pi), tsep_pipi(_tsep_pipi), Lt(_Lt), tK_dim(_tK_dim){}
  inline distributionType operator()(int const* coord, const distributionType &from) const{
    typedef iterate<distributionType> iter;
    int tK = coord[tK_dim]; int tB = (tK + tsep_k_pi + tsep_pipi) % Lt;
    distributionType out(from);
    for(int i=0;i<iter::size(out);i++)
      iter::at(i,out) = iter::at(i,out) * iter::at(i,bubble({tB}));
    return out;
  }
};



#endif
