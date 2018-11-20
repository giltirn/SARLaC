#ifndef _COMMON_DEFS_H
#define _COMMON_DEFS_H

//A set of useful typedefs - you don't have to use these but it makes things a bit cleaner
#include<config.h>
#include<utils/macros.h>

#include<distribution.h>
#include<data_series.h>
#include<tensors.h>

CPSFIT_START_NAMESPACE

//Distribution types
typedef distribution<double> distributionD;
typedef distribution<float> distributionF;

typedef rawDataDistribution<double> rawDataDistributionD;
typedef rawDataDistribution<float> rawDataDistributionF;

typedef jackknifeDistribution<double> jackknifeDistributionD;
typedef jackknifeDistribution<float> jackknifeDistributionF;

typedef jackknifeCdistribution<double> jackknifeCdistributionD;
typedef jackknifeCdistribution<float> jackknifeCdistributionF;

typedef doubleJackknifeDistribution<double> doubleJackknifeDistributionD;
typedef doubleJackknifeDistribution<float> doubleJackknifeDistributionF;

typedef bootstrapDistribution<double> bootstrapDistributionD;
typedef bootstrapDistribution<float> bootstrapDistributionF;


//Series types
typedef dataSeries<int, rawDataDistribution<double> > rawTimeSeriesD; //time, dist(value)
typedef dataSeries<int, rawDataDistribution<float> > rawTimeSeriesF;

typedef dataSeries<int, jackknifeDistribution<double> > jackknifeTimeSeriesD; 
typedef dataSeries<int, jackknifeDistribution<float> > jackknifeTimeSeriesF;

typedef dataSeries<int, jackknifeCdistribution<double> > jackknifeCtimeSeriesD; 
typedef dataSeries<int, jackknifeCdistribution<float> > jackknifeCtimeSeriesF;

typedef dataSeries<int, doubleJackknifeDistribution<double> > doubleJackknifeTimeSeriesD; 
typedef dataSeries<int, doubleJackknifeDistribution<float> > doubleJackknifeTimeSeriesF;

typedef filteredDataSeries<rawTimeSeriesD> filteredRawTimeSeriesD;
typedef filteredDataSeries<rawTimeSeriesF> filteredRawTimeSeriesF;

typedef filteredDataSeries<jackknifeTimeSeriesD> filteredJackknifeTimeSeriesD;
typedef filteredDataSeries<jackknifeTimeSeriesF> filteredJackknifeTimeSeriesF;

typedef filteredDataSeries<doubleJackknifeTimeSeriesD> filteredDoubleJackknifeTimeSeriesD;
typedef filteredDataSeries<doubleJackknifeTimeSeriesF> filteredDoubleJackknifeTimeSeriesF;

typedef correlationFunction<double,rawDataDistributionD> rawDataCorrelationFunctionD;
typedef correlationFunction<double,jackknifeDistributionD> jackknifeCorrelationFunctionD;
typedef correlationFunction<double,doubleJackknifeDistributionD> doubleJackknifeCorrelationFunctionD;

//Tensor types
typedef NumericSquareMatrix<jackknifeDistribution<double> >  jackknifeSquareMatrixD;
typedef NumericSquareMatrix<jackknifeDistribution<float> >  jackknifeSquareMatrixF;

CPSFIT_END_NAMESPACE

#endif
