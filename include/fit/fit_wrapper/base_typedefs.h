#ifndef _SARLAC_BASE_TYPEDEFS_H_
#define _SARLAC_BASE_TYPEDEFS_H_

//There are a number of basic types and compound types needed for the framework. These are controlled by a baseFitTypedefs, which takes a policy some user-defined input types

#include<config.h>
#include<utils/macros.h>
#include<tensors/numeric_vector.h>
#include<tensors/numeric_square_matrix.h>
#include<distribution/jackknife.h>
#include<data_series/correlationfunction.h>
#include<data_series/sample_series.h>

SARLAC_START_NAMESPACE

#define INHERIT_TYPEDEF(FROM,DEF) typedef typename FROM::DEF DEF
#define INHERIT_USING(FROM, DEF, TEMPL) template<typename TEMPL> using DEF = typename FROM::template DEF<TEMPL>

//BaseTypes must typedef the following:
#define INHERIT_INPUT_FIT_TYPEDEFS(FROM)			       \
  INHERIT_TYPEDEF(FROM,DistributionType); \
  INHERIT_TYPEDEF(FROM,CorrelationFunctionDistribution)

template<typename BaseTypes>
struct baseFitTypedefs{
  INHERIT_INPUT_FIT_TYPEDEFS(BaseTypes);
  
  typedef NumericSquareMatrix<DistributionType> MatrixDistribution;
  typedef NumericVector<DistributionType> VectorDistribution;

  typedef sampleSeries<const CorrelationFunctionDistribution> sampleSeriesType;
  typedef NumericSquareMatrixSampleView<const MatrixDistribution> sampleInvCorrType;
};

#define INHERIT_BASE_FIT_TYPEDEFS(FROM)			       \
  INHERIT_INPUT_FIT_TYPEDEFS(FROM);			       \
  INHERIT_TYPEDEF(FROM,MatrixDistribution); \
  INHERIT_TYPEDEF(FROM,VectorDistribution); \
  INHERIT_TYPEDEF(FROM,sampleSeriesType); \
  INHERIT_TYPEDEF(FROM,sampleInvCorrType)



//This is the standard implementation of BaseTypes for some generate coordinate type, with the data contained in a correlationFunction object comprising some kind of distributions (usually jackknife)
//Minimizer is assumed to take the cost function as its sole template parameter
template<typename GeneralizedCoordinate, typename _DistributionType = jackknifeDistribution<double> >
struct standardInputFitTypes{
  typedef _DistributionType DistributionType;
  typedef correlationFunction<GeneralizedCoordinate,DistributionType> CorrelationFunctionDistribution;
};


SARLAC_END_NAMESPACE
#endif
