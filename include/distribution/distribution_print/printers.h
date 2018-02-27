#ifndef _CPSFIT_PRINTERS_H_
#define _CPSFIT_PRINTERS_H_

//Objects that control how distributions are printed. All derive from a common base class, distributionPrinter

#include<utils/macros.h>
#include<output_format/cenerr.h>
#include<distribution/distribution_print/print_policy.h>

CPSFIT_START_NAMESPACE

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

//Error
template<typename DistributionType, typename ValuePolicy = printStats<DistributionType> >
struct errorDistributionPrinter: public distributionPrinter<DistributionType>{
  void print(std::ostream &os, const DistributionType &dist) const{
    os << ValuePolicy::error(dist);
  }
};

//Central value or error with control over rounding, sig figs, etc
template<typename DistributionType, typename ValuePolicy = printStats<DistributionType> >
struct publicationCenOrErrDistributionPrinter: public distributionPrinter<DistributionType>{
public:
  SigFigsSource src;

  int min_width; //pad with trailing spaces if width < min_width. Use 0 for no padding
  int sci_fmt_threshold; //when the operative value's exponent is > this value we will use scientific format, otherwise we will use decimal format

  //Optional: set the number of significant figures and which of the two quantities this refers to
  bool set_sig_figs;
  int nsf; //number of sig figs

  //Optional: set the power-of-10 at which the output is rounded. Overrides sig-figs option
  bool set_round_pow;
  int round_pow;
  
  //Optional: set the exponents to a fixed value. Disables sci_fmt_threshold
  bool set_exponent;  
  int set_exponent_to;
public:
  void print(std::ostream &os_out, const DistributionType &d) const{
    std::stringstream os;    
    std::ios::streampos init_pos = os.tellp();
    typedef decltype(ValuePolicy::centralValue(d)) valueType;
    valueType val = src == Central ? ValuePolicy::centralValue(d) : ValuePolicy::error(d);
    decimal dval(val);

    if(set_round_pow){
      dval = dval.roundAtPow(round_pow);
      dval = dval.truncateAtPow(round_pow);
    }else if(set_sig_figs){
      int sig_fig_pow = dval.firstSigFigPow();
      int rp = sig_fig_pow - nsf + 1;
      dval = dval.roundAtPow(rp);
      dval = dval.truncateAtPow(rp);
    }
    
    if(set_exponent){
      dval = dval.setExp(set_exponent_to);
    }else{ //check if we should convert to decimal format
      int abs_pow10_op = abs(decimal::pow10(val));
      if(abs_pow10_op <= sci_fmt_threshold)
	dval = dval.setExp(0);
    }
    
    dval.print(os, false,true);
    if(dval.exponent() != 0) os << "\\times 10^{" << dval.exponent() << "}";

    int width = os.tellp() - init_pos;

    for(int i=0;i<min_width - width;i++) os << ' ';

    os_out << os.rdbuf();
  }

  publicationCenOrErrDistributionPrinter(const SigFigsSource _src = Central,const int _nsf = 3): src(_src), set_sig_figs(true), nsf(_nsf),  min_width(0), sci_fmt_threshold(3), set_exponent(false), set_round_pow(false){
    assert(src == Central || src == Error);    
  }

  inline void setSigFigs(const int _nsf){ set_sig_figs = true; nsf = _nsf; }
  inline void setMinWidth(const int _min_width){ min_width = _min_width; }
  inline void setSciFormatThreshold(const int p){ sci_fmt_threshold = p; }
  inline void setExponent(const int p){ set_exponent = true; set_exponent_to = p; }
  inline void setRoundPower(const int p){ set_round_pow = true; round_pow = p; }
};


//Publication format : Central value (error)
template<typename DistributionType, typename ValuePolicy = printStats<DistributionType> >
struct publicationDistributionPrinter: public distributionPrinter<DistributionType>{
public:
  int min_width; //pad with trailing spaces if width < min_width. Use 0 for no padding
  int sci_fmt_threshold; //when the operative value's exponent is > this value we will use scientific format, otherwise we will use decimal format

  //Optional: set the number of significant figures and which of the two quantities this refers to
  bool set_sig_figs;
  int nsf; //number of sig figs
  SigFigsSource sfsrc; //whether the sig.figs specified is based on the error or the central value

  //Optional: set the power-of-10 at which the output is rounded. Overrides sig-figs option
  bool set_round_pow;
  int round_pow;
  
  //Optional: set the exponents to a fixed value. Disables sci_fmt_threshold
  bool set_exponent;  
  int set_exponent_to;
public:
  void print(std::ostream &os_out, const DistributionType &d) const{
    std::stringstream os;    
    std::ios::streampos init_pos = os.tellp();
    typedef decltype(ValuePolicy::centralValue(d)) valueType;
    valueType mu = ValuePolicy::centralValue(d);
    valueType err = ValuePolicy::error(d);

    cenErr ce(mu,err);

    if(set_round_pow)
      ce.setRoundPower(round_pow);
    else if(set_sig_figs)
      ce.setSigFigs(nsf,sfsrc);
    
    if(set_exponent){
      ce.setExp(set_exponent_to);
    }else{ //check if we should convert to decimal format
      int abs_pow10_op;
      {
	SigFigsSource sff = sfsrc;
	if(sff == Largest) sff = (mu >= err ? Central : Error); 
	abs_pow10_op = abs(decimal::pow10(sff == Central ? mu : err));
      }
      if(abs_pow10_op <= sci_fmt_threshold)
	ce.convertDecimal();
    }
    
    ce.printPublication(os);

    int width = os.tellp() - init_pos;

    for(int i=0;i<min_width - width;i++) os << ' ';

    os_out << os.rdbuf();
  }

  publicationDistributionPrinter(const int _nsf = 3, const SigFigsSource _sfsrc = Largest): set_sig_figs(true), nsf(_nsf), sfsrc(_sfsrc), min_width(0), sci_fmt_threshold(3), set_exponent(false), set_round_pow(false){}

  inline void setSigFigs(const int _nsf, const SigFigsSource _sfsrc = Largest){ set_sig_figs = true; nsf = _nsf; sfsrc  = _sfsrc; }
  inline void setMinWidth(const int _min_width){ min_width = _min_width; }
  inline void setSciFormatThreshold(const int p){ sci_fmt_threshold = p; }
  inline void setExponent(const int p){ set_exponent = true; set_exponent_to = p; }
  inline void setRoundPower(const int p){ set_round_pow = true; round_pow = p; }
};

CPSFIT_END_NAMESPACE
#endif
