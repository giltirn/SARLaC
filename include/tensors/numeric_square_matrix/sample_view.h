#ifndef _SARLAC_NUMERIC_SQUARE_MATRIX_SAMPLE_VIEW_H_
#define _SARLAC_NUMERIC_SQUARE_MATRIX_SAMPLE_VIEW_H_

//Create a "view" of a square matrix of distribution type, picking out a single sample

#include<config.h>
#include<utils/macros.h>
#include<tensors/numeric_square_matrix/class.h>

SARLAC_START_NAMESPACE

template<typename NumericSquareMatrixType>
class NumericSquareMatrixSampleView{
  typedef typename _get_elem_type<NumericSquareMatrixType>::type DistributionType;
  typedef typename std::remove_const<typename std::remove_reference<decltype( ((DistributionType*)(NULL))->sample(0) )>::type>::type SampleType;

  NumericSquareMatrixType &M;
  int sample;
public:
  NumericSquareMatrixSampleView(NumericSquareMatrixType &_M, const int _sample): M(_M), sample(_sample){}
  
  inline int size() const{ return M.size(); }

  inline const SampleType& operator()(const int i, const int j) const{ return iterate<DistributionType>::at(sample, M(i,j)); }
  
  template<typename U = NumericSquareMatrixType>
  inline typename std::enable_if< !std::is_const<U>::value, SampleType& >::type operator()(const int i, const int j){ return iterate<DistributionType>::at(sample, M(i,j)); }
};

SARLAC_END_NAMESPACE
#endif
