#ifndef _FILTERED_DATA_SERIES_CLASS_H
#define _FILTERED_DATA_SERIES_CLASS_H

//Generic filter interface for a data series, allowing for example restriction of data for fitting
//Prepend const to SeriesType if const access required

#include<vector>

#include<config.h>
#include<utils/template_wizardry.h>
#include<utils/macros.h>

SARLAC_START_NAMESPACE

template<typename SeriesType>
class filteredDataSeries{
public:
  typedef typename SeriesType::GeneralizedCoordinate GeneralizedCoordinate;
  typedef typename SeriesType::DataType DataType;
  std::vector<int> subset_map;
  SeriesType &series;
public:
  template<typename Filter> 
  filteredDataSeries(SeriesType &_series,const Filter &filter):series(_series){
    for(int i=0;i<series.size();i++){
      if(filter.accept(series.coord(i), series.value(i)))
	subset_map.push_back(i);      
    }
  }
  inline int size() const{ return subset_map.size(); }
  inline typename add_const_if<DataType, SeriesType>::type & value(const int i) const{ return series.value(subset_map[i]); }
  inline typename add_const_if<GeneralizedCoordinate, SeriesType>::type & coord(const int i) const{ return series.coord(subset_map[i]); }
};

SARLAC_END_NAMESPACE

#endif
