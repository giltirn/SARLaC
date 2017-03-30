#ifndef _CPSFIT_DISTRIBUTION_PRINT_H
#define _CPSFIT_DISTRIBUTION_PRINT_H

#include<cmath>
#include<cfenv>
#include<iostream>
#include<cassert>
#include<sstream>


//Base type for printers that use the mean and standard error of a distribution
template<typename Derived>
class printerBase{
  std::ostream &os;
 public:
  printerBase(std::ostream &_os = std::cout): os(_os){}

  template<typename DistributionType, typename std::enable_if< hasSampleMethod<DistributionType>::value && hasMeanMethod<DistributionType>::value, int>::type = 0>
  Derived &operator<<(const DistributionType &d){
    Derived* td = (Derived*)(this);
    td->print(os,d);
    return *td;
  }
  template<typename OtherType, typename std::enable_if< !hasSampleMethod<OtherType>::value || !hasMeanMethod<OtherType>::value, int>::type = 0>
  Derived &operator<<(const OtherType &d){
    Derived* td = (Derived*)(this);
    os << d;
    return *td;
  }
  Derived & operator<<( std::ostream& (*fp)(std::ostream&) ){
    Derived* td = (Derived*)(this);
    os << fp;
    return *td;
  }

};
class basicPrint: public printerBase<basicPrint>{
private:
  friend class printerBase<basicPrint>;
  template<typename Dist>
  void print(std::ostream &os, const Dist &d){ os << "(" << d.mean() << " +- " << d.standardError() << ")"; }
public:
  basicPrint(std::ostream &_os = std::cout): printerBase<basicPrint>(_os){}
};

class publicationPrint: public printerBase<publicationPrint>{
public:
  enum SigFigsSource { Central, Error, Largest };
  
private:
  int nsf; //number of sig figs
  SigFigsSource sfsrc; //whether the sig.figs specified is based on the error or the central value
  int min_width; //pad with trailing spaces if width < min_width. Use 0 for no padding
  
  friend class printerBase<publicationPrint>;
  template<typename Dist>
  void print(std::ostream &os_out, const Dist &d){
    // std::ios oldState(nullptr);
    // oldState.copyfmt(os);
    const auto oldRound = std::fegetround();

    std::stringstream os;
    
    std::ios::streampos init_pos = os.tellp();
    
    typedef decltype(d.mean()) valueType;
    valueType mu = d.mean();
    valueType err = d.standardError();

    SigFigsSource src = sfsrc;
    if(src == Largest) src = fabs(mu) > fabs(err) ? Central : Error;
    //os << "(" << mu << " +- " << err << ") -> ";

    os << std::fixed;

    valueType &srcv = src == Central ? mu : err;
    std::fesetround(FE_DOWNWARD); 

    int pow10 = int( rint(log10(fabs(srcv))) );
    
    std::fesetround(FE_TONEAREST);
    
    valueType coeffpow10mu = mu / pow(10.,pow10);
    valueType coeffpow10err = err / pow(10.,pow10);

    os.precision(nsf-1);
    if(coeffpow10mu >= 0.) os << ' '; //equal spacing if +/-
    os << coeffpow10mu << "(" << coeffpow10err << ")";    
    
    if(pow10 < 0) os << "*10^{" << pow10 << "}";
    else if(pow10 > 0 && pow10 < 10) os << "*10^" << pow10;
    else if(pow10 > 10) os << "*10^{" << pow10 << "}";

    int width = os.tellp() - init_pos;

    for(int i=0;i<min_width - width;i++) os << ' ';
    
    std::fesetround(oldRound);

    os_out << os.rdbuf();
    
    //os.copyfmt(oldState);
  }
public:
  publicationPrint(const int _nsf = 3, const SigFigsSource _sfsrc = Largest, std::ostream &_os = std::cout): printerBase<publicationPrint>(_os), nsf(_nsf), sfsrc(_sfsrc), min_width(0){}

  inline void setSigFigs(const int _nsf, const SigFigsSource _sfsrc = Largest){ nsf = _nsf; sfsrc  = _sfsrc; }
  inline void setMinWidth(const int _min_width){ min_width = _min_width; }
};




#endif
