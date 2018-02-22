#ifndef _CPSFIT_MPL_DATA_SETS_2D_H_
#define _CPSFIT_MPL_DATA_SETS_2D_H_

#include<iostream>
#include<vector>

#include<config.h>
#include<utils/macros.h>
#include<plot/plot/write_python.h>

CPSFIT_START_NAMESPACE

enum SetType{ DataSetType, ErrorBandType, HistogramType };

//Contain and write data for a series of data points with errors
class PythonDataContainer{
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> dxm;
  std::vector<double> dxp;
  std::vector<double> dym;
  std::vector<double> dyp;
  std::string set_tag;

public:

  template<typename Data>
  void import(const Data &data, const std::string &tag){
    set_tag = tag;
    
    int sz = data.size();
    for(int i=0;i<sz;i++){
      x.push_back(data.x(i));
      y.push_back(data.y(i));
      
      dxm.push_back(data.dxm(i));
      dxp.push_back(data.dxp(i));
      
      dym.push_back(data.dym(i));
      dyp.push_back(data.dyp(i));
    }
  }

  void write(std::ostream &os) const{
    os << set_tag << "= pyplot.DataSet()\n";
    os << set_tag << ".x = " << ListPrint<double>(x) << '\n';
    os << set_tag << ".y = " << ListPrint<double>(y) << '\n';
    os << set_tag << ".dxm = " << ListPrint<double>(dxm) << '\n';
    os << set_tag << ".dxp = " << ListPrint<double>(dxp) << '\n';
    os << set_tag << ".dym = " << ListPrint<double>(dym) << '\n';
    os << set_tag << ".dyp = " << ListPrint<double>(dyp) << "\n\n";
  }

  const std::string &tag() const{ return set_tag; }
  
};

//Contain and write data for an error band
class PythonErrorBandContainer{
  std::vector<double> x;
  std::vector<double> upper;
  std::vector<double> lower;
  std::string set_tag;

public:
  
  template<typename Band>
  void import(const Band &band, const std::string &tag){
    set_tag = tag;
    
    int sz = band.size();
    for(int i=0;i<sz;i++){
      x.push_back(band.x(i));
      upper.push_back(band.upper(i));
      lower.push_back(band.lower(i));
    }
  }

  void write(std::ostream &os) const{
    os << set_tag << "= pyplot.ErrorBand()\n";
    os << set_tag << ".x = " << ListPrint<double>(x) << '\n';
    os << set_tag << ".upper = " << ListPrint<double>(upper) << '\n';
    os << set_tag << ".lower = " << ListPrint<double>(lower) << '\n';
  }

  const std::string &tag() const{ return set_tag; }  
};

//Contain and write data for a histogram
class PythonHistogramContainer{
  std::vector<double> y;
  std::string set_tag;

public:
  
  template<typename Data>
  void import(const Data &data, const std::string &tag){
    set_tag = tag;
    
    int sz = data.size();
    for(int i=0;i<sz;i++){
      y.push_back(data.y(i));
    }
  }

  void write(std::ostream &os) const{
    os << set_tag << "= pyplot.HistogramData()\n";
    os << set_tag << ".y = " << ListPrint<double>(y) << '\n';
  }

  const std::string &tag() const{ return set_tag; }  
};

CPSFIT_END_NAMESPACE
#endif
