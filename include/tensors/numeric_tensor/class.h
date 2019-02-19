#ifndef _CPSFIT_NUMERIC_TENSOR_CLASS_H_
#define _CPSFIT_NUMERIC_TENSOR_CLASS_H_

#include<iostream>
#include<vector>

#include<config.h>
#include<utils/macros.h>
#include<ET/generic_ET.h>
#include<tensors/numeric_tensor/helper_structs.h>

CPSFIT_START_NAMESPACE

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

  template<int R = Rank, typename std::enable_if<R == 2, int>::type = 0>
  inline DataType &operator()(const int i, const int j){ int ij[2] = {i,j}; return (*this)(ij); }
  
  template<int R = Rank, typename std::enable_if<R == 2, int>::type = 0>
  inline const DataType &operator()(const int i, const int j) const{ int ij[2] = {i,j}; return (*this)(ij); }


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

CPSFIT_END_NAMESPACE
#endif
