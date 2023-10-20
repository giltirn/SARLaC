#ifndef _CPSFIT_PLOT_ACCESSORS_2D_H_
#define _CPSFIT_PLOT_ACCESSORS_2D_H_

//Some common accessor classes for feeding data into the 2D plotters
#include<cassert>
#include<vector>
#include<cmath>

#include<config.h>
#include<utils/macros.h>
#include<distribution/distribution_iterate.h>

CPSFIT_START_NAMESPACE


//Example accessor for x,y data in std::vectors with symmetric errors
//Can also be used to define error bands
class DataVectorAccessor{
  const std::vector<double> &_x;
  const std::vector<double> &_y;
  const std::vector<double> &_dx;
  const std::vector<double> &_dy;
  int sz;
public:
  DataVectorAccessor(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &dx, const std::vector<double> &dy): _x(x), _y(y), _dx(dx), _dy(dy), sz(y.size()){
    assert(x.size() == sz && dx.size() == sz && dy.size() == sz);
  }

  inline double x(const int i) const{ return _x[i]; }
  inline double y(const int i) const{ return _y[i]; }
  inline double dxm(const int i) const{ return _dx[i]; }
  inline double dxp(const int i) const{ return _dx[i]; }
  inline double dym(const int i) const{ return _dy[i]; }
  inline double dyp(const int i) const{ return _dy[i]; }

  inline double upper(const int i) const{ return _y[i]+_dy[i]; }
  inline double lower(const int i) const{ return _y[i]-_dy[i]; }
  
  inline int size() const{ return sz; }
};

class DataVectorInternalAccessor{
  std::vector<double> _x;
  std::vector<double> _y;
  std::vector<double> _dx;
  std::vector<double> _dy;
  int sz;
public:
  DataVectorInternalAccessor(): sz(0){}
  DataVectorInternalAccessor(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &dx, const std::vector<double> &dy): _x(x), _y(y), _dx(dx), _dy(dy), sz(y.size()){
    assert(x.size() == sz && dx.size() == sz && dy.size() == sz);
  }

  void push_back(const double x, const double y, const double dx, const double dy){ 
    _x.push_back(x);
    _y.push_back(y);
    _dx.push_back(dx);
    _dy.push_back(dy);
    ++sz;
  }

  inline double x(const int i) const{ return _x[i]; }
  inline double y(const int i) const{ return _y[i]; }
  inline double dxm(const int i) const{ return _dx[i]; }
  inline double dxp(const int i) const{ return _dx[i]; }
  inline double dym(const int i) const{ return _dy[i]; }
  inline double dyp(const int i) const{ return _dy[i]; }

  inline double upper(const int i) const{ return _y[i]+_dy[i]; }
  inline double lower(const int i) const{ return _y[i]-_dy[i]; }
  
  inline int size() const{ return sz; }
};





//Example accessors for error band with x, upper, lower separate std::vectors
class BandVectorAccessor{
  const std::vector<double> &_x;
  const std::vector<double> &_upper;
  const std::vector<double> &_lower;
  int sz;
public:
  BandVectorAccessor(const std::vector<double> &x, const std::vector<double> &upper, const std::vector<double> &lower): _x(x), _upper(upper), _lower(lower), sz(upper.size()){
    assert(lower.size() == sz && x.size() == sz);
  }
  inline double x(const int i) const{ return _x[i]; }
  inline double upper(const int i) const{ return _upper[i]; }
  inline double lower(const int i) const{ return _lower[i]; }
  inline int size() const{ return sz; }
};

//Give a constant error band across some range
template<typename DistributionType, typename ValuePolicy>
class BandRangeConstantDistributionValue: public ValuePolicy{
  const int min;
  const int max;
  double hi;
  double lo;
public:
  BandRangeConstantDistributionValue(const int _min, const int _max, const DistributionType &d): min(_min), max(_max){
    double cen = this->ValuePolicy::value(d);
    hi = cen + this->ValuePolicy::errplus(d); 
    lo = cen - this->ValuePolicy::errminus(d); 
  }
  inline double x(const int i) const{ return min + i; }
  inline double upper(const int i) const{ return hi; }
  inline double lower(const int i) const{ return lo; }
  inline int size() const{ return max-min+1; }
};
  
  
//Use to plot all elements of a distribution type
template<typename DistributionType>
class DistributionSampleAccessor{
  const DistributionType &d;
public:
  DistributionSampleAccessor(const DistributionType &_d): d(_d){}

  inline double x(const int i) const{ return i; }
  inline double y(const int i) const{ return iterate<DistributionType>::at(i,d); }
  inline double dxm(const int i) const{ return 0; }
  inline double dxp(const int i) const{ return 0; }
  inline double dym(const int i) const{ return 0; }
  inline double dyp(const int i) const{ return 0; }

  inline double upper(const int i) const{ return y(i); }
  inline double lower(const int i) const{ return y(i); }
  
  inline int size() const{ return iterate<DistributionType>::size(d); }
};

//A generic accessor for a data series type (any container with a size() and operator() ), with the element accesses controlled by policies:
//CoordinatePolicy is converts the underlying coordinate type into double central values and errors
//ValuePolicy does the same for the y data
template<typename DataSeriesType, typename CoordinatePolicy, typename ValuePolicy>
class DataSeriesAccessor: public ValuePolicy, public CoordinatePolicy{
  const DataSeriesType &series;
public:
  DataSeriesAccessor(const DataSeriesType &_series):series(_series){}
  
  int size() const{ return series.size(); }

  inline double x(const int i) const{ return this->CoordinatePolicy::value(series.coord(i)); }
  inline double dxp(const int i) const{ return this->CoordinatePolicy::errplus(series.coord(i)); }
  inline double dxm(const int i) const{ return this->CoordinatePolicy::errminus(series.coord(i)); }  

  inline double y(const int i) const{ return this->ValuePolicy::value(series.value(i)); }
  inline double dyp(const int i) const{ return this->ValuePolicy::errplus(series.value(i)); }
  inline double dym(const int i) const{ return this->ValuePolicy::errminus(series.value(i)); }  

  inline double upper(const int i) const{ return y(i)+dyp(i); }
  inline double lower(const int i) const{ return y(i)-dym(i); }
};

//Accessor for a std::vector<Distribution> or other array type that has an operator[] and size(). Coordinates are just the array indices 
template<typename DistributionArrayType, typename ValuePolicy>
class DistributionArrayAccessor: public ValuePolicy{
  const DistributionArrayType &v;
public:
  DistributionArrayAccessor(const DistributionArrayType &_v):v(_v){}
  
  int size() const{ return v.size(); }

  inline double x(const int i) const{ return i; }
  inline double dxp(const int i) const{ return 0; }
  inline double dxm(const int i) const{ return 0; }  

  inline double y(const int i) const{ return this->ValuePolicy::value(v[i]); }
  inline double dyp(const int i) const{ return this->ValuePolicy::errplus(v[i]); }
  inline double dym(const int i) const{ return this->ValuePolicy::errminus(v[i]); }  

  inline double upper(const int i) const{ return y(i)+dyp(i); }
  inline double lower(const int i) const{ return y(i)-dym(i); }
};

//Example ValuePolicy types
template<typename DistributionType> //specialize if necessary
class DistributionPlotAccessor{  
public:
  static inline double value(const DistributionType &d){ return d.mean(); }
  static inline double errplus(const DistributionType &d){ return d.standardError(); }
  static inline double errminus(const DistributionType &d){ return d.standardError(); }  
};

template<typename DistributionType>
class DistributionPlotAccessorNanRepl{
  bool replace_nan;
  double replace_nan_with;

  inline double repl(const double v) const { return (replace_nan && std::isnan(v)) ? replace_nan_with : v; }

public:
  DistributionPlotAccessorNanRepl(): replace_nan(false){}
  
  void replaceNan(const double with = 0.){ replace_nan = true; replace_nan_with = with; }
  
  double value(const DistributionType &d) const{ return repl(d.mean()); }
  double errplus(const DistributionType &d) const{ return repl(d.standardError()); }
  double errminus(const DistributionType &d) const{ return repl(d.standardError()); }  
};

template<typename T>
class ScalarValueAccessor{
public:
  static inline double value(const T &d){ return d; }
  static inline double errplus(const T &d){ return 0; } //zero error
  static inline double errminus(const T &d){ return 0; }  
};

//Example CoordinatePolicy types
//Coordinate accessor that relies on implicit conversion of type to double
template<typename T>
class ScalarCoordinateAccessor{
public:
  static inline double value(const T &d){ return d; }
  static inline double errplus(const T &d){ return 0; } //zero error
  static inline double errminus(const T &d){ return 0; }  
};

template<typename DistributionType> //specialize if necessary
class DistributionCoordinateAccessor{
public:
  static inline double value(const DistributionType &d){ return d.mean(); }
  static inline double errplus(const DistributionType &d){ return d.standardError(); }
  static inline double errminus(const DistributionType &d){ return d.standardError(); }  
};


//Virtual base class accessors for reducing the number of function definitionst
template<typename T>
class CurveDataAccessorBase{
public:

  virtual double x(const int i) const = 0;
  virtual double y(const int i) const = 0;
  double dxm(const int i) const{ return 0; }
  double dxp(const int i) const{ return 0; }
  double dym(const int i) const{ return 0; }
  double dyp(const int i) const{ return 0; }

  double upper(const int i) const{ return this->y(i); }
  double lower(const int i) const{ return this->y(i); }
  
  virtual int size() const = 0;
};



CPSFIT_END_NAMESPACE
#endif
