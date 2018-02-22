#ifndef _SAMPLE_SERIES_H_
#define _SAMPLE_SERIES_H_

//An accessor class that allows access to individual samples of a series of distributions with different coordinates in some general array
//This is used because the fitter expects single values for each coordinate
//Prepend const to seriesType if const access required

#include<utils/macros.h>
#include<utils/template_wizardry.h>
#include<distribution/distribution_iterate.h>

CPSFIT_START_NAMESPACE

template<typename SeriesType>
class sampleSeries{
public:
  typedef typename SeriesType::DataType ElementType; 
  typedef typename SeriesType::GeneralizedCoordinate GeneralizedCoordinate;
  typedef typename ElementType::DataType DataType;
private:
  const int sample;
  SeriesType &series;
public:
  sampleSeries(SeriesType &_series, const int _sample):series(_series), sample(_sample){ }
  inline int size() const{ return series.size(); }
  inline typename add_const_if<DataType, SeriesType>::type & value(const int i) const{ return iterate<ElementType>::at(sample, series.value(i)); }
  inline typename add_const_if<GeneralizedCoordinate, SeriesType>::type & coord(const int i) const{ return series.coord(i); }
};

CPSFIT_END_NAMESPACE
#endif
