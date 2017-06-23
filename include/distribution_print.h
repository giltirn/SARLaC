#ifndef _CPSFIT_DISTRIBUTION_PRINT_H
#define _CPSFIT_DISTRIBUTION_PRINT_H

#include<cmath>
#include<cfenv>
#include<iostream>
#include<cassert>
#include<sstream>

//Allow manual compile-time override of source of central value and error by type
template<typename DistributionType>
struct printStats{
  inline static auto centralValue(const DistributionType &d)->decltype(d.best()){ return d.best(); }
  inline static auto error(const DistributionType &d)->decltype(d.standardError()){ return d.standardError(); }
};

//Objects that control how distributions are printed
template<typename DistributionType>
struct distributionPrinter{
  virtual void print(std::ostream &os, const DistributionType &dist) const = 0;
};

//( Central value +- error )
template<typename DistributionType, typename ValuePolicy = printStats<DistributionType> >
struct basicDistributionPrinter: public distributionPrinter<DistributionType>{
  void print(std::ostream &os, const DistributionType &dist) const{
    os << "(" << ValuePolicy::centralValue(dist) << " +- " << ValuePolicy::error(dist) << ")";
  }
};

//Central value
template<typename DistributionType, typename ValuePolicy = printStats<DistributionType> >
struct centralValueDistributionPrinter: public distributionPrinter<DistributionType>{
  void print(std::ostream &os, const DistributionType &dist) const{
    os << ValuePolicy::centralValue(dist);
  }
};

//Central value (error)
template<typename DistributionType, typename ValuePolicy = printStats<DistributionType> >
struct publicationDistributionPrinter: public distributionPrinter<DistributionType>{
public:
  enum SigFigsSource { Central, Error, Largest };
  
private:
  int nsf; //number of sig figs
  SigFigsSource sfsrc; //whether the sig.figs specified is based on the error or the central value
  int min_width; //pad with trailing spaces if width < min_width. Use 0 for no padding

public:
  void print(std::ostream &os_out, const DistributionType &d) const{
    // std::ios oldState(nullptr);
    // oldState.copyfmt(os);
    const auto oldRound = std::fegetround();

    std::stringstream os;
    
    std::ios::streampos init_pos = os.tellp();
    
    typedef decltype(ValuePolicy::centralValue(d)) valueType;
    valueType mu = ValuePolicy::centralValue(d);
    valueType err = ValuePolicy::error(d);

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

  publicationDistributionPrinter(const int _nsf = 3, const SigFigsSource _sfsrc = Largest): nsf(_nsf), sfsrc(_sfsrc), min_width(0){}

  inline void setSigFigs(const int _nsf, const SigFigsSource _sfsrc = Largest){ nsf = _nsf; sfsrc  = _sfsrc; }
  inline void setMinWidth(const int _min_width){ min_width = _min_width; }
};

//A class that stores a singleton copy of the printer for a given type. The current printer is used in the stream operators for the distributions. The printer can be overridden at arbitrary time
template<typename DistributionType>
struct distributionPrint{
  static distributionPrinter<DistributionType>* printer(distributionPrinter<DistributionType>* change = NULL, bool delete_old = true){
    static distributionPrinter<DistributionType>* p = NULL;
    static bool initialized = false;
    if(!initialized){ p = new basicDistributionPrinter<DistributionType>; initialized = true; }

    if(change != NULL){ if(p!=NULL && delete_old) delete p; p = change; }
    return p;
  }
};


#endif
