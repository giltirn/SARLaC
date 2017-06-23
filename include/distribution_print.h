#ifndef _CPSFIT_DISTRIBUTION_PRINT_H
#define _CPSFIT_DISTRIBUTION_PRINT_H

#include<cmath>
#include<cfenv>
#include<iostream>
#include<cassert>
#include<sstream>

/* template<typename DistributionType> */
/* struct distributionPrinter{ */
/*   virtual void print(std::ostream &os, const DistributionType &dist) const = 0; */
/* }; */



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
    ::operator<<(os,d);
    //os << d;
    return *td;
  }
  Derived & operator<<( std::ostream& (*fp)(std::ostream&) ){
    Derived* td = (Derived*)(this);
    os << fp;
    return *td;
  }

};

//Allow manual override of source of central value and error by type
template<typename DistributionType>
struct printStats{
  inline static auto centralValue(const DistributionType &d)->decltype(d.best()){ return d.best(); }
  inline static auto error(const DistributionType &d)->decltype(d.standardError()){ return d.standardError(); }
};
template<typename T>
struct printStats< doubleJackknifeDistribution<T> >{
  inline static std::string centralValue(const doubleJackknifeDistribution<T> &d){
    std::ostringstream os; os << "[";
    for(int s=0;s<d.size()-1;s++) os << d.sample(s).best() << ", ";
    os << d.sample(d.size()-1).best() << "]";
    return os.str();
  }
  inline static std::string error(const doubleJackknifeDistribution<T> &d){ 
    std::ostringstream os; os << "[";
    for(int s=0;s<d.size()-1;s++) os << d.sample(s).standardError() << ", ";
    os << d.sample(d.size()-1).standardError() << "]";
    return os.str();
  }

};

template< template<typename> class ValuePolicy = printStats>
class basicPrint: public printerBase<basicPrint<ValuePolicy> >{
private:
  friend class printerBase<basicPrint<ValuePolicy> >;
  template<typename Dist>
  void print(std::ostream &os, const Dist &d){ os << "(" << printStats<Dist>::centralValue(d) << " +- " << printStats<Dist>::error(d) << ")"; }
public:
  basicPrint(std::ostream &_os = std::cout): printerBase<basicPrint<ValuePolicy> >(_os){}
};

template< template<typename> class ValuePolicy = printStats>
class publicationPrint: public printerBase<publicationPrint<ValuePolicy> >{
public:
  enum SigFigsSource { Central, Error, Largest };
  
private:
  int nsf; //number of sig figs
  SigFigsSource sfsrc; //whether the sig.figs specified is based on the error or the central value
  int min_width; //pad with trailing spaces if width < min_width. Use 0 for no padding
  
  friend class printerBase<publicationPrint<ValuePolicy> >;
  template<typename Dist>
  void print(std::ostream &os_out, const Dist &d){
    // std::ios oldState(nullptr);
    // oldState.copyfmt(os);
    const auto oldRound = std::fegetround();

    std::stringstream os;
    
    std::ios::streampos init_pos = os.tellp();
    
    typedef decltype(ValuePolicy<Dist>::centralValue(d)) valueType;
    valueType mu = ValuePolicy<Dist>::centralValue(d);
    valueType err = ValuePolicy<Dist>::error(d);

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
  publicationPrint(const int _nsf = 3, const SigFigsSource _sfsrc = Largest, std::ostream &_os = std::cout): printerBase<publicationPrint<ValuePolicy> >(_os), nsf(_nsf), sfsrc(_sfsrc), min_width(0){}

  inline void setSigFigs(const int _nsf, const SigFigsSource _sfsrc = Largest){ nsf = _nsf; sfsrc  = _sfsrc; }
  inline void setMinWidth(const int _min_width){ min_width = _min_width; }
};




#endif
