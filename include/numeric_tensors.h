#ifndef _NUMERIC_TENSORS_H
#define _NUMERIC_TENSORS_H

#include<map>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include<template_wizardry.h>
#include<generic_ET.h>
#include<hdf5_serialize.h>
#include<distribution.h>

template<typename Numeric>
class NumericVector{
  std::vector<Numeric> v;

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & v;
  }
public:
  NumericVector():v(){}
  explicit NumericVector(const int n):v(n){}
  NumericVector(const int n, const Numeric &def):v(n,def){}
  NumericVector(const NumericVector &r) = default;
  NumericVector(NumericVector &&r) = default;
  
  template<typename Initializer> //Initializer is a lambda-type with operator()(const int)
  inline NumericVector(const int n, const Initializer &initializer): v(n){
    for(int i=0;i<n;i++) v[i] = initializer(i);
  }
  
  NumericVector & operator=(const NumericVector &r) = default;
  NumericVector & operator=(NumericVector &&r) = default;

  ENABLE_GENERIC_ET(NumericVector, NumericVector<Numeric>, NumericVector<Numeric>);
  
  int size() const{ return v.size(); }
  
  void resize(const int n){
    v.resize(n);
  }
  void resize(const int n, const Numeric &init){
    v.resize(n,init);
  }
  template<typename Initializer>
  inline void resize(const int n, const Initializer &initializer){
    v.resize(n);
    for(int i=0;i<n;i++) v[i] = initializer(i);
  }
  
  void zero(){ for(int i=0;i<v.size();i++) v[i] = 0.; }

  std::string print() const{
    std::ostringstream os;
    os << v[0];
    for(int i=1;i<v.size();i++)
      os << " " << v[i];
    return os.str();
  }

  Numeric & operator()(const int i){ return v[i]; }
  const Numeric & operator()(const int i) const { return v[i]; }

  Numeric & operator[](const int i){ return v[i]; }
  const Numeric & operator[](const int i) const { return v[i]; }
  
  NumericVector &operator+=(const NumericVector &r){
    for(int i=0;i<v.size();i++) v[i] += r.v[i];
  }

  inline void push_back(const Numeric &b){ v.push_back(b); }

  GENERATE_HDF5_SERIALIZE_METHOD((v));
};
#ifdef HAVE_HDF5
template<typename D>
inline void write(HDF5writer &writer, const NumericVector<D> &d, const std::string &tag){ d.write(writer,tag); }
template<typename D>
inline void read(HDF5reader &reader, NumericVector<D> &d, const std::string &tag){ d.read(reader,tag); }
#endif


template<typename Numeric> 
std::ostream & operator<<(std::ostream & stream, const NumericVector<Numeric> &vec){
  stream << "(";
  for(int i=0;i<vec.size();i++)
    stream << vec[i] << (i != vec.size()-1 ? " " : ")");
  return stream;
}

//Vector dot product
template<typename T>
T dot(const NumericVector<T> &a, const NumericVector<T> &b){
  assert(a.size() == b.size());
  T out(a(0)); zeroit(out);
  for(int i=0;i<a.size();i++) out = out + a(i) * b(i);
  return out;
}
template<typename T>
T mod2(const NumericVector<T> &m){
  return dot(m,m);
}

template<typename Numeric>
class NumericSquareMatrix{ //square matrix
  std::vector<std::vector<Numeric> > m;

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & m;
  }
public:
  NumericSquareMatrix():m(){}
  explicit NumericSquareMatrix(const int n): m(n, std::vector<Numeric>(n)){}
  NumericSquareMatrix(const int n, const Numeric &init): m(n, std::vector<Numeric>(n,init)){}
  NumericSquareMatrix(const NumericSquareMatrix &r) = default;
  NumericSquareMatrix(NumericSquareMatrix &&r) = default;

  template<typename Initializer> //Initializer is a lambda-type with operator()(const int)
  inline NumericSquareMatrix(const int n, const Initializer &initializer): m(n, std::vector<Numeric>(n)){
    for(int i=0;i<n;i++)
      for(int j=0;j<n;j++)
	m[i][j] = initializer(i,j);
  }
  
  NumericSquareMatrix & operator=(const NumericSquareMatrix &r) = default;
  NumericSquareMatrix & operator=(NumericSquareMatrix &&r) = default;
  
  typedef NumericSquareMatrix<Numeric> ET_tag;
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,NumericSquareMatrix<Numeric> >::value, int>::type = 0>
  NumericSquareMatrix(U&& expr): NumericSquareMatrix(expr.common_properties()){
#pragma omp parallel for
    for(int i=0;i<this->size()*this->size();i++)
      getElem<NumericSquareMatrix<Numeric> >::elem(*this, i) = expr[i];
  }
  
  inline int size() const{ return m.size(); }
  
  void resize(const int n){
    if(m.size() == n) return;
    m.resize(n);
    for(int i=0;i<n;i++)
      m[i].resize(n);      
  }
  void resize(const int n, const Numeric &init){
    if(m.size() == n) return;
    m.resize(n);
    for(int i=0;i<n;i++)
      m[i].resize(n,init);     
  }
  template<typename Initializer>
  inline void resize(const int n, const Initializer &initializer){
    m.resize(n);
    for(int i=0;i<n;i++){
      m[i].resize(n);
      for(int j=0;j<n;j++) m[i][j] = initializer(i,j);
    }
  }
  
  void zero(){
    const int n = m.size();
    for(int i=0;i<n;i++)
      for(int j=0;j<n;j++)
	m[i][j] = 0.;
  }

  std::string print() const{
    std::ostringstream os;
    for(int i=0;i<m.size();i++){
      os << m[i][0];
      for(int j=1;j<m[i].size();j++)
  	os << " " << m[i][j];
      os << std::endl;
    }
    return os.str();
  }
  
  Numeric & operator()(const int i, const int j){ return m[i][j]; }
  const Numeric & operator()(const int i, const int j) const { return m[i][j]; }

  GENERATE_HDF5_SERIALIZE_METHOD((m));
};
#ifdef HAVE_HDF5
template<typename D>
inline void write(HDF5writer &writer, const NumericSquareMatrix<D> &d, const std::string &tag){ d.write(writer,tag); }
template<typename D>
inline void read(HDF5reader &reader, NumericSquareMatrix<D> &d, const std::string &tag){ d.read(reader,tag); }
#endif

//Modulus operator for matrices
template<typename T>
T mod2(const NumericSquareMatrix<T> &m){
  assert(m.size() > 0);
  T out(m(0,0));
  zeroit(out);
  for(int i=0;i<m.size();i++)
    for(int j=0;j<m.size();j++)
      out = out + m(i,j)*m(i,j); //only correct for real matrices!
  return out;
}


template<typename Numeric, typename StreamType, typename std::enable_if< isStreamType<StreamType>::value, int>::type = 0> 
StreamType & operator<<(StreamType & stream, const NumericSquareMatrix<Numeric> &mat){
  for(int i=0;i<mat.size();i++){
    for(int j=0;j<mat.size();j++){
      stream << mat(i,j) << " ";
    }
    stream << '\n';
  }
  return stream;
}
		     



template<typename NumericSquareMatrixType>
class NumericSquareMatrixSampleView{
  typedef typename _get_elem_type<NumericSquareMatrixType>::type DistributionType;
  typedef typename std::remove_const<typename std::remove_reference<decltype( ((DistributionType*)(NULL))->sample(0) )>::type>::type SampleType;

  NumericSquareMatrixType &M;
  int sample;
public:
  NumericSquareMatrixSampleView(NumericSquareMatrixType &_M, const int _sample): M(_M), sample(_sample){}
  
  inline int size() const{ return M.size(); }

  inline const SampleType& operator()(const int i, const int j) const{ return M(i,j).sample(sample); }
  
  template<typename U = NumericSquareMatrixType>
  inline typename std::enable_if< !std::is_const<U>::value, SampleType& >::type operator()(const int i, const int j){ return M(i,j).sample(sample); }
};


//A vector type with an additional global mapping between some generic tag and it's elements
template<typename T, typename MapType>
class mappedVector: public NumericVector<T>{
  //MapType must conform to the following form
  // struct MapTypeExample{
  //   typedef std::string tagType;    //choose a tag type - here a string but can be anything
  //   inline int map(const tagType &tag) const; //mapping between tag and index
  //   inline tagType unmap(const int idx) const; //unmapping between index and tag
  //   inline int size() const; //number of elements
  // };

  MapType const* mapping;

public:
  typedef typename MapType::tagType tagType;

  mappedVector(const MapType & _mapping): mapping(&_mapping), NumericVector<T>(_mapping.size()){}
  
  typedef mappedVector<T,MapType> ET_tag;
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,mappedVector<T,MapType> >::value, int>::type = 0>
  mappedVector(U&& expr){
    mapping = expr.common_properties();
    this->resize(mapping->size());
    for(int i=0;i<mapping->size();i++) this->operator()(i) = expr[i];
  }
  
  inline T & operator()(const tagType &tag){ return this->operator()(mapping->map(tag)); }
  inline const T & operator()(const tagType &tag) const{ return this->operator()(mapping->map(tag)); }
  inline T & operator()(const int idx){ return NumericVector<T>::operator()(idx); }
  inline const T & operator()(const int idx) const{ return NumericVector<T>::operator()(idx); }
  
  const MapType & getMapping() const{ return *mapping; }
};
template<typename T, typename MapType>
struct getElem<mappedVector<T,MapType> >{
  static inline auto elem(const mappedVector<T,MapType> &v, const int i)->decltype(v(i)){ return v(i); }
  // static inline auto elem(mappedVector<T,MapType> &v, const int i)->decltype(v(i)){ return v(i); }
  static inline MapType const* common_properties(const mappedVector<T,MapType> &v){ return &v.getMapping(); }
};

template<typename T, typename MapType>
inline void debug_print(const mappedVector<T,MapType> &v){ std::cout << &v.getMapping() << std::endl; std::cout.flush(); }


class stringTagMap{
  std::map<std::string,int> mp;
  std::map<int,std::string> ump;
public:
  typedef std::string tagType; 
  inline int map(const tagType &tag) const{
    auto it = mp.find(tag);
    assert(it != mp.end());
    return it->second;
  }
  inline tagType unmap(const int idx) const{
    auto it = ump.find(idx);
    assert(it != ump.end());
    return it->second;
  }
  void add(const std::string &tag, const int idx){
    mp[tag] = idx;
    ump[idx] = tag;
  }
  inline int size() const{ return mp.size(); }
};





//A true general tensor type
template<typename T,int Rank,int D>
struct _tensor_helper{
  static inline size_t vol(const int* sizes){
    return sizes[0]*_tensor_helper<T,Rank,D-1>::vol(sizes+1);
  }
  static inline size_t map(const size_t accum, const int *elem, const int *sizes){
    return _tensor_helper<T,Rank,D-1>::map( elem[0] + accum*sizes[0], elem+1,sizes+1);
  }
  static inline void unmap(int *into, const size_t off, const int *sizes, const size_t vol){
    size_t subvol = vol/sizes[0];
    *into = off / subvol;
    _tensor_helper<T,Rank,D-1>::unmap(into+1, off % subvol, sizes+1, subvol);
  }
};
template<typename T,int Rank>
struct _tensor_helper<T,Rank,-1>{
  static inline size_t vol(const int* lst){
    return 1;
  }
  static inline size_t map(const size_t accum, const int *elem, const int *sizes){
    return accum;
  }
  static inline void unmap(int *into, const size_t off, const int *sizes, const size_t vol){};
};

template<typename T, int R>
class NumericTensor;

template<typename T, int R>
struct _NTprinter{
  static void print(std::ostream &os, const NumericTensor<T,R> &t);
};




template<typename DataType, int Rank>
class NumericTensor{
  typedef _tensor_helper<DataType,Rank,Rank-1> helper;
  std::vector<int> dsizes;
  size_t vol;
  std::vector<DataType> data;

  inline size_t map(const int *elem) const{ return helper::map(0,elem,dsizes.data()); }
  inline void unmap(int* into, const size_t off) const{ return helper::unmap(into,off,dsizes.data(),vol); }

  friend struct getElem<NumericTensor<DataType,Rank> >;
  friend struct _NTprinter<DataType,Rank>;

  template<typename T, int R>
  friend class NumericTensor;
public:
  inline NumericTensor() = default;
  inline explicit NumericTensor(int const* _size): dsizes(_size,_size+Rank), vol(helper::vol(_size)), data(vol){ }
  inline explicit NumericTensor(const std::vector<int> &_size): NumericTensor(_size.data()){}
  inline explicit NumericTensor(std::initializer_list<int> _size): NumericTensor(_size.begin()){}

  
  inline NumericTensor(int const* _size, const DataType &init): dsizes(_size,_size+Rank), vol(helper::vol(_size)), data(vol,init){ }
  inline NumericTensor(const std::vector<int> &_size, const DataType &init): NumericTensor(_size.data(),init){}
  inline NumericTensor(std::initializer_list<int> _size, const DataType &init): NumericTensor(_size.begin(),init){}
  
  
  inline NumericTensor(const NumericTensor<DataType,Rank> &r) = default;
  inline NumericTensor(NumericTensor<DataType,Rank> &&r) = default;

  
  template<typename Initializer>
  inline NumericTensor(int const*_size, const Initializer &init): NumericTensor(_size){
    for(size_t i=0;i<vol;i++){
      int oe[Rank]; this->unmap(oe,i);
      data[i] = init(oe); //Initializer must have DataType operator()(int const* coord)
    }
  }
  
  template<typename Initializer>
  inline NumericTensor(const std::vector<int> &_size, const Initializer &init): NumericTensor(_size.data(),init){}

  template<typename Initializer>
  inline NumericTensor(std::initializer_list<int> _size, const Initializer &init): NumericTensor(_size.begin(),init){}

      
  typedef NumericTensor<DataType,Rank> ET_tag;
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,ET_tag>::value, int>::type = 0>
  NumericTensor(U&& expr){
    this->resize(expr.common_properties().sizes);
    for(size_t i=0;i<this->vol;i++)
      data[i] = expr[i];
  }

  inline NumericTensor& operator=(const NumericTensor<DataType,Rank> &r) = default;
  inline NumericTensor& operator=(NumericTensor<DataType,Rank> &&r) = default;
  
   //will invalidate existing values
  inline void resize(int const* _size){
    dsizes = std::vector<int>(_size,_size+Rank);
    vol = helper::vol(_size);
    data.resize(vol);
  }
  inline void resize(const std::vector<int> & _size){ return resize(_size.data()); } 
  inline void resize(std::initializer_list<int> _size){ return resize(_size.begin()); } 

  inline void resize(int const* _size, const DataType &init){
    dsizes = std::vector<int>(_size,_size+Rank);
    vol = helper::vol(_size);
    data.resize(vol,init);
  }
  inline void resize(const std::vector<int> & _size, const DataType &init){ return resize(_size.data(),init); } 
  inline void resize(std::initializer_list<int> _size, const DataType &init){ return resize(_size.begin(),init); } 

  inline const std::vector<DataType> & internalVector() const{ return data; }
  
  inline const DataType &operator()(std::initializer_list<int> elem) const{
    return data[map(elem.begin())];
  }
  inline DataType &operator()(std::initializer_list<int> elem){
    return data[map(elem.begin())];
  }
  inline const DataType &operator()(int const* elem) const{
    return data[map(elem)];
  }
  inline DataType &operator()(int const* elem){
    return data[map(elem)];
  }  
  
  inline const std::vector<int> &sizes() const{ return dsizes; }
  inline int size(const int dir) const{ return dsizes[dir]; }

  template<typename Rule>
  auto transform(const Rule &rule = Rule()) const -> NumericTensor<typename std::decay<decltype(rule( (int const*)(NULL), data[0]))>::type,Rank>{ 
    NumericTensor<typename std::decay<decltype(rule((int const*)(NULL), data[0]))>::type,Rank> out(dsizes.data());
    for(size_t i=0;i<vol;i++){
      int oe[Rank]; out.unmap(oe,i);
      out.data[i] = rule(oe,data[i]); //Rule must have NewDataType operator()(int const* coord, const DataType &from)
    }
    return out;
  }

  template<typename Rule>
  NumericTensor<DataType,Rank-1> reduce(const int dim, const Rule &rule = Rule()) const{
    int new_size[Rank-1];
    int i=0;
    for(int ii=0;ii<Rank;ii++) if(ii != dim) new_size[i++] = dsizes[ii];    
    NumericTensor<DataType,Rank-1> out(new_size);
    for(size_t i=0;i<out.vol;i++){
      int oe[Rank-1]; out.unmap(oe,i);
      rule(out.data[i], oe, *this); //Rule must have operator()(DataType &into, int const* coord, const NumericTensor<DataType,Rank> &from)
    }
    return out;
  }
  
};

//Contract two tensors over a single index. Preserves ordering of remaining indices
//eg A_ijkl B_mnio -> C_jklmno

template<typename T, int R1, int R2>
NumericTensor<T,R1+R2-2> contract(const NumericTensor<T,R1> &A, const NumericTensor<T,R2> &B, const int rankA, const int rankB){
  if(A.size(rankA) != B.size(rankB)) error_exit(std::cout << "NumericTensor contract(...) contracting over ranks with different sizes!\n");
  const int conRankSize = A.size(rankA);
  int dimC[R1+R2-2];
  int c=0;
  for(int i=0;i<R1;i++)
    if(i!=rankA) dimC[c++] = A.size(i);
  for(int i=0;i<R2;i++)
    if(i!=rankB) dimC[c++] = B.size(i);
  
  return NumericTensor<T,R1+R2-2>(dimC,
				  [&](const int* cC){
				    int cA[R1], cB[R2];
				    int c=0;
				    for(int i=0;i<R1;i++)
				      if(i!=rankA) cA[i]=cC[c++];
				    for(int i=0;i<R2;i++)
				      if(i!=rankB) cB[i]=cC[c++];

				    T out; zeroit(out);
				    for(int i=0;i<conRankSize;i++){
				      cA[rankA] = cB[rankB] = i;				      
				      out = out + A(cA)*B(cB);
				    }
				    return out;
				  });
}


template<typename T, int R>
void _NTprinter<T,R>::print(std::ostream &os, const NumericTensor<T,R> &t){
  int e[R];
  for(size_t i=0;i<t.vol;i++){
    t.unmap(e,i);
    os << '(';
    for(int d=0;d<R-1;d++) os << e[d] << ',';
    os << e[R-1] << ") : " << t.data[i] << std::endl;
  }
};

//Matrix specialization
template<typename T>
struct _NTprinter<T,2>{
  static void print(std::ostream &os, const NumericTensor<T,2> &t){
    for(int i=0;i<t.size(0);i++){
      for(int j=0;j<t.size(1)-1;j++){
	os << t({i,j}) << ", ";
      }
      os << t({i,t.size(1)-1}) << "\n";
    }
  }
};
//Vector specialization
template<typename T>
struct _NTprinter<T,1>{
  static void print(std::ostream &os, const NumericTensor<T,1> &t){
    for(int j=0;j<t.size(0)-1;j++){
      os << t({j}) << ", ";
    }
    os << t({t.size(0)-1});
  }
};

template<typename T, int R>
inline std::ostream & operator<<(std::ostream &os, const NumericTensor<T,R> &t){
  _NTprinter<T,R>::print(os,t);
  return os;
}


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


#include<numeric_tensors_ET.h>

#endif
