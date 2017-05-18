#ifndef _COMMON_DEFS_H
#define _COMMON_DEFS_H

//A set of useful typedefs - you don't have to use these but it makes things a bit cleaner
#include<distribution.h>
#include<data_series.h>
#include<numeric_tensors.h>

//Distribution types
typedef distribution<double> distributionD;
typedef distribution<float> distributionF;

typedef jackknifeDistribution<double> jackknifeDistributionD;
typedef jackknifeDistribution<float> jackknifeDistributionF;

typedef jackknifeCdistribution<double> jackknifeCdistributionD;
typedef jackknifeCdistribution<float> jackknifeCdistributionF;

typedef doubleJackknifeDistribution<double> doubleJackknifeDistributionD;
typedef doubleJackknifeDistribution<float> doubleJackknifeDistributionF;

//Series types
typedef dataSeries<int, distribution<double> > rawTimeSeriesD; //time, dist(value)
typedef dataSeries<int, distribution<float> > rawTimeSeriesF;

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

//Tensor types
typedef NumericMatrix<jackknifeDistribution<double> >  jackknifeMatrixD;
typedef NumericMatrix<jackknifeDistribution<float> >  jackknifeMatrixF;



#endif
