#ifndef _SARLAC_MPL_DATA_SETS_2D_H_
#define _SARLAC_MPL_DATA_SETS_2D_H_

#include<iostream>
#include<vector>

#include<config.h>
#include<utils/macros.h>
#include<plot/plot/write_python.h>

SARLAC_START_NAMESPACE

enum SetType{ DataSetType, ErrorBandType, HistogramType, ErrorLineType, ErrorCurveType };

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


class PythonErrorLineContainer{
  typedef PythonTuple<double> TupleD;

  std::vector<TupleD> start; //array of tuples giving the start coordinates of each line
  std::vector<TupleD> end; //array of tuples giving the end coordinates of each line
  std::string set_tag;

public:
  
  template<typename Data>
  void import(const Data &data, const std::string &tag){
    set_tag = tag;
    
    int sz = data.size();
    for(int i=0;i<sz;i++){
      start.push_back(data.start(i));
      end.push_back(data.end(i));
    }
  }

  void write(std::ostream &os) const{
    os << set_tag << "= pyplot.ErrorLine()\n";
    os << set_tag << ".start = " << ListPrint<TupleD>(start) << '\n';
    os << set_tag << ".end = " << ListPrint<TupleD>(end) << '\n';
  }

  const std::string &tag() const{ return set_tag; }  
};

class PythonErrorCurveContainer{
  typedef PythonTuple<double> TupleD;

  std::vector<std::vector<TupleD> > curves; //an array of tuples for each curve describing several points along the curve
  std::vector<TupleD> markers; //array of tuples giving the coordinates of markers (optional)
  std::string set_tag;

public:
  
  template<typename Data>
  void import(const Data &data, const std::string &tag){
    set_tag = tag;
   
    for(int i=0;i<data.ncurves();i++)
      curves.push_back(data.curves(i));
    
    for(int i=0;i<data.nmarkers();i++)
      markers.push_back(data.markers(i));
  }

  void write(std::ostream &os) const{
    os << set_tag << "= pyplot.ErrorCurve()\n";
    os << set_tag << ".curves = [";
    for(int c=0;c<curves.size();c++){
      os << ListPrint<TupleD>(curves[c]);
      if(c != curves.size()-1) os << ", ";
    }
    os << "]\n";
    if(markers.size() != 0){   
      os << set_tag << ".markers = " << ListPrint<TupleD>(markers) << '\n';
    }else{
      os << set_tag << ".markers = None\n";
    }

  }

  const std::string &tag() const{ return set_tag; }  
};


SARLAC_END_NAMESPACE
#endif
