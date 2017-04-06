#ifndef _SUPERJACKKNIFE_H_
#define _SUPERJACKKNIFE_H_

#include <distribution.h>
#include <generic_ET.h>

class superJackknifeLayout{
  std::vector<std::string > ens_tags;
  std::vector<int> ens_sizes;
  std::vector< std::pair<int,int> > ens_sample_map;
  int total_size;
public:
  superJackknifeLayout(): total_size(0){}
  
  //returns -1 for non-existent
  inline int ensIdx(const std::string &e) const{
    for(int i=0;i<ens_tags.size();i++) if(e == ens_tags[i]) return i;
    return -1;
  }

  void addEnsemble(const std::string &tag, const int size){
    int ens_idx = ens_tags.size();
    ens_tags.push_back(tag);
    ens_sizes.push_back(size);
    for(int i=0;i<size;i++) ens_sample_map.push_back(std::pair<int,int>(ens_idx,i) );
    total_size += size;
  }

  inline int nSamplesTotal() const{ return total_size; }

  inline int nEnsembles() const{ return ens_tags.size(); }

  inline int nSamplesEns(const int ens_idx) const{ return ens_sizes[ens_idx]; }
  inline int nSamplesEns(const std::string &ens_tag) const{ return nSamplesEns(ensIdx(ens_tag)); }
  
  inline const std::pair<int,int> & sampleMap(int gsample) const{
    return ens_sample_map[gsample];
  }
};

template<typename _DataType>
class superJackknifeDistribution{
public:
  typedef _DataType DataType;
  typedef superJackknifeDistribution<DataType> Generic_ET_base_type;
protected:
  std::vector< jackknifeDistribution<DataType> > ens_jacks;
  DataType avg;
  superJackknifeLayout const* layout;

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & ens_jacks;
    ar & avg;
  }

public:
  superJackknifeDistribution(const superJackknifeLayout &_layout, const DataType &mean): layout(&_layout), ens_jacks(_layout.nEnsembles()){
    for(int i=0;i<ens_jacks.size();i++) ens_jacks[i].resize( layout->nSamplesEns(i), mean ); 
    avg = mean;
  }

  superJackknifeDistribution(const superJackknifeLayout &_layout, const int first_ens_idx, const jackknifeDistribution<DataType> &first): superJackknifeDistribution(_layout, first.mean()) {
    assert(first.size() == layout->nSamplesEns(first_ens_idx));
    ens_jacks[first_ens_idx] = first;
  }
  
  superJackknifeDistribution(const superJackknifeDistribution &r): layout(r.layout), ens_jacks(r.ens_jacks), avg(r.avg){}
  superJackknifeDistribution(superJackknifeDistribution&& o) noexcept : layout(o.layout), ens_jacks(std::move(o.ens_jacks)), avg(o.avg){ }

  template<typename Expr, IS_EXPRESSION_WITH_GENERIC_BASE_TYPE(Expr, superJackknifeDistribution<DataType>, superJackknifeDistribution<DataType>) > \
  superJackknifeDistribution(const Expr &e){    
    layout = &e.getLayout();
    ens_jacks.resize(layout->nEnsembles());
    for(int i=0;i<ens_jacks.size();i++) ens_jacks[i].resize( layout->nSamplesEns(i) ); 
    for(int i=0;i<this->size();i++) this->sample(i) = e.sample(i);
    avg = e.mean();
  }
  
  superJackknifeDistribution & operator=(const superJackknifeDistribution &r){ layout = r.layout; ens_jacks = r.ens_jacks; avg = r.avg; return *this; }
  
  int size() const{ return layout->nSamplesTotal(); }
  
  const DataType & sample(const int idx) const{
    return ens_jacks[ layout->sampleMap(idx).first ].sample( layout->sampleMap(idx).second );
  }
  DataType & sample(const int idx){
    return ens_jacks[ layout->sampleMap(idx).first ].sample( layout->sampleMap(idx).second );
  }

  DataType mean() const{ return avg; }

  DataType standardError() const{
    DataType v = 0;
    for(int i=0;i<ens_jacks.size();i++){
      DataType vi = ens_jacks[i].standardError();
      v = v + vi*vi;
    }
    return sqrt(v);
  }

  void setEnsembleJackknife(const int ens_idx, const jackknifeDistribution<DataType> &j){
    assert(j.mean() == avg);
    assert(j.size() == layout->nSamplesEns(ens_idx));
    ens_jacks[ens_idx] = j;    
  }
  inline void setEnsembleJackknife(const std::string & ens_tag, const jackknifeDistribution<DataType> &j){
    setEnsembleJackknife(layout->ensIdx(ens_tag),j);
  }

  inline const jackknifeDistribution<DataType> & getEnsembleJackknife(const int ens_idx) const{
    return ens_jacks[ens_idx];
  }
  inline const jackknifeDistribution<DataType> & getEnsembleJackknife(const std::string &ens_tag) const{
    return ens_jacks[layout->ensIdx(ens_tag)];
  }
  
  const superJackknifeLayout & getLayout() const{ return *layout; }
};






#define DEF_SJACK_BINOP(CLASS, OP, OP_ACTION_AB, OP_ACTION_MEAN)			\
  template<typename T, typename U, BASE_TYPES_EQUAL_TO(T,U,superJackknifeDistribution<typename T::DataType>) > \
class CLASS{ \
public: \
  typedef typename T::Generic_ET_base_type Generic_ET_base_type; \
  typedef typename T::DataType DataType; \
  T const& a; \
  U const& b; \
  CLASS(const T &aa, const U &bb): a(aa), b(bb){ assert(&a.getLayout() == &b.getLayout()); } \
  \
  inline auto sample(const int i) const ->decltype(OP_ACTION_AB){ return OP_ACTION_AB; } \
  inline auto mean() const ->decltype(OP_ACTION_MEAN){ return OP_ACTION_MEAN; } \
  inline int size() const{ return a.size(); } \
  inline const superJackknifeLayout & getLayout() const{ return a.getLayout(); } \
}; \
\
template<typename T, typename U, BASE_TYPES_EQUAL_TO(T,U,superJackknifeDistribution<typename T::DataType>) >  \
CLASS<T,U> OP(const T &a, const U &b){		\
  return CLASS<T,U>(a,b); \
}  

DEF_SJACK_BINOP(ETsjackPlus, operator+, a.sample(i) + b.sample(i), a.mean() + b.mean() );
DEF_SJACK_BINOP(ETsjackMinus, operator-, a.sample(i) - b.sample(i), a.mean() - b.mean() );
DEF_SJACK_BINOP(ETsjackTimes, operator*, a.sample(i)*b.sample(i), a.mean()*b.mean() );
DEF_SJACK_BINOP(ETsjackDivide, operator/, a.sample(i)/b.sample(i), a.mean()/b.mean() );





//The below are broken into a base and 2 ops because it doesn't always make sense to perform a vector-scalar binop in both orders  eg  vector/scalar but not scalar/vector

//vector-scalar type
#define DEF_SJACK_BINOP_VS_BASE(CLASS, OP, OP_ACTION_AB, OP_ACTION_MEAN) \
template<typename T, typename U, BASE_TYPE_EQUAL_TO_VECTORSCALAR(T,U,superJackknifeDistribution<typename T::DataType>) > \
class CLASS{ \
public: \
  typedef typename T::Generic_ET_base_type Generic_ET_base_type; \
  typedef typename T::DataType DataType; \
  T const& a; \
  U const& b; \
  CLASS(const T &aa, const U &bb): a(aa), b(bb){ } \
  \
  inline auto sample(const int i) const ->decltype(OP_ACTION_AB){ return OP_ACTION_AB; } \
  inline auto mean() const ->decltype(OP_ACTION_MEAN){ return OP_ACTION_MEAN; } \
  inline int size() const{ return a.size(); } \
  inline const superJackknifeLayout & getLayout() const{ return a.getLayout(); } \
};

//T is vector-type, U is scalar type
#define DEF_SJACK_BINOP_VS(CLASS,OP)		     \
template<typename T, typename U, BASE_TYPE_EQUAL_TO_VECTORSCALAR(T,U,superJackknifeDistribution<typename T::DataType>)>  \
CLASS<T,U> OP(const T &a, const U &b){		\
  return CLASS<T,U>(a,b); \
}

//T is scalar-type, U is vector type
#define DEF_SJACK_BINOP_SV(CLASS,OP)		     \
template<typename T, typename U, BASE_TYPE_EQUAL_TO_VECTORSCALAR(U,T,superJackknifeDistribution<typename U::DataType>) >  \
CLASS<U,T> OP(const T &a, const U &b){		\
  return CLASS<U,T>(b,a); \
}

DEF_SJACK_BINOP_VS_BASE( ETsjackScalarMult, operator*, a.sample(i)*b, a.mean()*b);
DEF_SJACK_BINOP_VS(ETsjackScalarMult, operator*);
DEF_SJACK_BINOP_SV(ETsjackScalarMult, operator*);
DEF_SJACK_BINOP_VS_BASE( ETsjackScalarDiv, operator/, a.sample(i)/b, a.mean()/b);
DEF_SJACK_BINOP_VS( ETsjackScalarDiv, operator/); 


#define DEF_SJACK_UNOP(CLASS, OP, OP_ACTION, OP_ACTION_MEAN)				\
template<typename T, BASE_TYPE_EQUAL_TO_UNOP(T,superJackknifeDistribution<typename T::DataType>) > \
class CLASS{ \
public: \
  typedef typename T::Generic_ET_base_type Generic_ET_base_type; \
  typedef typename T::DataType DataType;			 \
  T const& a; \
  CLASS(const T &aa): a(aa){ } \
  \
  inline auto sample(const int i) const ->decltype(OP_ACTION){ return OP_ACTION; } \
  inline auto mean() const ->decltype(OP_ACTION_MEAN){ return OP_ACTION_MEAN; } \
  inline int size() const{ return a.size(); } \
  inline const superJackknifeLayout & getLayout() const{ return a.getLayout(); } \
}; \
template<typename T, BASE_TYPE_EQUAL_TO_UNOP(T,superJackknifeDistribution<typename T::DataType>) >  \
CLASS<T> OP(const T &a){		\
  return CLASS<T>(a); \
}

DEF_SJACK_UNOP( ETsjackSqrt, sqrt, sqrt(a.sample(i)), sqrt(a.mean()) );


#endif
