#ifndef _FILTERED_DATA_SERIES_FILTERS_H
#define _FILTERED_DATA_SERIES_FILTERS_H

#include<config.h>
#include<data_series/filtered_data_series/class.h>

SARLAC_START_NAMESPACE

//Only include coordinates within (outside if invert active) some range.
template<typename GeneralizedCoordinate>
class filterXrange{
  //accept if   x_min<=x<=x_max  
  GeneralizedCoordinate x_min;
  GeneralizedCoordinate x_max;
  bool _invert;
public:
  filterXrange(const GeneralizedCoordinate &_x_min, const GeneralizedCoordinate &_x_max): x_min(_x_min), x_max(_x_max), _invert(false){}

  void invert(){ _invert = !_invert; } //switch accept inside/outside of range
  
  template<typename T>
  bool accept(const GeneralizedCoordinate &x, const T &y) const{
    bool cond = x >= x_min && x <= x_max;
    return _invert ? !cond : cond;    
  }
};

SARLAC_END_NAMESPACE

#endif
