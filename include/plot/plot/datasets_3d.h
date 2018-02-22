#ifndef _CPSFIT_MPL_DATA_SETS_3D_H_
#define _CPSFIT_MPL_DATA_SETS_3D_H_

#include<iostream>
#include<vector>

#include<config.h>
#include<utils/macros.h>
#include<plot/plot/write_python.h>

CPSFIT_START_NAMESPACE

enum SetType3D{ WireframeType, ScatterType };

class Python3DdataContainer{
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::string set_tag;

public:
  
  template<typename Accessor>
  void import(const Accessor &acc, const std::string &tag){
    set_tag = tag;
    
    int sz = acc.size();
    for(int i=0;i<sz;i++){
      x.push_back(acc.x(i));
      y.push_back(acc.y(i));
      z.push_back(acc.z(i));
    }
  }

  void write(std::ostream &os) const{
    os << set_tag << "= pyplot.DataSet3D()\n";
    os << set_tag << ".x = " << ListPrint<double>(x) << '\n';
    os << set_tag << ".y = " << ListPrint<double>(y) << '\n';
    os << set_tag << ".z = " << ListPrint<double>(z) << '\n';
  }

  const std::string &tag() const{ return set_tag; }  
};

CPSFIT_END_NAMESPACE
#endif
