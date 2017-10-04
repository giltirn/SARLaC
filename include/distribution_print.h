#ifndef _CPSFIT_DISTRIBUTION_PRINT_H
#define _CPSFIT_DISTRIBUTION_PRINT_H

#include<cmath>
#include<cfenv>
#include<iostream>
#include<cassert>
#include<sstream>
#include<deque>

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


//A class that stores a number in a decimal representation
class decimal{
  std::deque<int> v;
  int dp;
  int base_pow;
  int sgn;
    
  void setup(double val, int init_sz = 16){
    v.resize(0);
    base_pow = pow10(val);
    val/=pow(10.,base_pow); //n.m e exp
    if(val<0.){
      val = -val; sgn = -1;
    }else sgn = 1;
    
    char tmp[init_sz+2]; //leave room for trailing \0 and decimal pt
    snprintf(tmp,init_sz+2,"%.50f",val);

    for(int i=0;i<init_sz+1;i++){
      int d = (int)tmp[i] - 48;
      if(d>=0){
	if(d>9){ std::cout << "Unexpected digit '" << tmp[i] << "' in " << tmp << std::endl; exit(-1); }
	
	assert(d<=9);
	v.push_back(d);
	if(v.size() == init_sz) break;
      }
    }
    while(v.size() < init_sz) v.push_back(0);

    dp = 0;
  }
public:

  //Get the base 10 exponent of the value v
  static inline int pow10(const double v){ return v == 0. ? 0 : int(floor(log10(fabs(v)))); }
  
  decimal(const double d, int init_sz = 16){
    setup(d,init_sz);
  }

  void print(std::ostream &os, bool incl_exponent = true, bool leading_space_if_positive = false) const{
    if(leading_space_if_positive) os << ( sgn == -1 ? '-' : ' ' );
    else if(sgn==-1) os << '-';
    
    for(int i=0;i<v.size();i++){
      os << v[i];
      if(i==dp && i!=v.size()-1) os << '.';
    }
    if(incl_exponent && base_pow !=0 ) os << 'e' << base_pow;
  }

  double value() const{
    std::stringstream ss;
    print(ss);
    double o;
    ss >> o;
    return o;
  }

  //Vector access of individual digits
  inline const int nDigits() const{ return v.size(); }
  inline const int operator[](const int i) const{ return v[i]; }
  
  int exponent() const{ return base_pow; }

  //Right-shift = divide by 10  (keep base power the same)
  decimal operator>>(const unsigned int n) const{
    decimal out(*this);
    for(int i=0;i<n;i++){
      --out.dp;
      if(out.dp<0){
	out.v.push_front(0);
	++out.dp;
      }	
    }
    return out;
  }
  //Left-shift = multiply by 10  (keep base power the same)
  decimal operator<<(const unsigned int n) const{
    decimal out(*this);
    for(unsigned int i=0;i<n;i++){
      ++out.dp;
      if(out.dp == out.v.size()-1)
	out.v.push_back(0);
    }
    return out;
  }

  //Set the base power
  decimal setExp(int p) const{
    decimal out(*this);
    int n = out.base_pow - p;
    if(n>0) out = out << n;
    else if(n<0) out = out >> -n;

    out.base_pow = p;
    return out;
  }
  
  //Round the numerical representation of the vector of integers on index i
  decimal roundAtPow(const int p) const{
    decimal out(*this);
    int pow_first_digit = out.base_pow + out.dp;
    int round_idx = pow_first_digit - p;
    if(round_idx < 0){
      out = out.setExp(p);
      round_idx = 0;
    }
    while(round_idx >= out.v.size()-1){
      out.v.push_back(0);
    }
    if(out.v[round_idx+1]>=5){
      int b=round_idx;
      while(out.v[b]==9){
  	out.v[b] = 0;
  	b--;
  	if(b<0) break;
      };
      if(b>=0) ++out.v[b];
      else{
	out.v.push_front(1);
	++out.dp;
      }
    }
    for(int b=round_idx+1;b<out.v.size();b++) out.v[b]=0;
    return out;
  }

  //Remove all digits after power p
  decimal truncateAtPow(const int p) const{
    decimal out(*this);
    int pow_first_digit = out.base_pow + out.dp;
    int last_idx = pow_first_digit - p;
    if(last_idx < 0){
      out = out.setExp(p);
      last_idx = 0;
    }
    while(last_idx >= out.v.size()-1){
      out.v.push_back(0);
    }
    out.v.resize(last_idx + 1);
    return out;
  }
    
  //Power of 10 at which the first sig fig resides
  int firstSigFigPow() const{
    int pow_first_digit = base_pow + dp;
    
    for(int i=0;i<v.size();i++){
      if(v[i]!=0){
	return pow_first_digit - i;
      }
    }
    return pow_first_digit;
  }
    
  std::vector<int> sigFigs(const int n, int *dp_idx = NULL) const{
    assert(n>0);
    std::vector<int> out(n,0);
    bool start = false;
    int o=0;
    for(int i=0;i<v.size();i++){
      if(!start && v[i]!=0){
	start = true;
	if(dp_idx != NULL) *dp_idx = dp - i; //location of decimal point relative to first sig fig
      }
      if(start){
	out[o++] = v[i];
	if(o==n) break;
      }
    }
    return out; 
  }
};

inline std::ostream & operator<<(std::ostream &os, const decimal &dec){
  dec.print(os); return os;
}

enum SigFigsSource { Central, Error, Largest };

//A class that contains a central value and error in decimal representation for printing purposes
class cenErr{
  decimal cen;
  decimal err;
  int nsf;
public:
  
  std::ostream & printBasic(std::ostream &os) const{
    os << cen << " +/- " << err; return os;
  }
  std::ostream & printPublication(std::ostream &os) const{
    cen.print(os, false,true);
    os << "(";
    if(err.value()!=0. && err[0] == 0){
      assert(err.nDigits() == cen.nDigits());
      //Remove leading zeroes
      bool start = false;
      for(int i=0;i<err.nDigits();i++){
	if(!start && err[i]!=0) start=true;	
	if(start) os << err[i];
      }
    }else err.print(os, false,false);
    os << ")";
    if(cen.exponent() != 0) os << "\\times 10^{" << cen.exponent() << "}";
    return os;
  }
  
  cenErr(double c, double e, SigFigsSource sf = Largest, const int _nsf = 3): cen(c), err(e), nsf(_nsf){
    assert(e>=0.);
    if(sf == Largest)
      sf = fabs(c) >= fabs(e) ? Central : Error;
    
    //Put them under a common exponent
    int exp = cen.exponent() < err.exponent() ? cen.exponent() : err.exponent();
    cen = cen.setExp(exp);
    err = err.setExp(exp);
 
    int sig_fig_pow_cen = cen.firstSigFigPow();
    int sig_fig_pow_err = err.firstSigFigPow();

    int round_pow = (sf == Central ? sig_fig_pow_cen : sig_fig_pow_err) - nsf + 1;

    //Round and truncate
    cen = cen.roundAtPow(round_pow);
    err = err.roundAtPow(round_pow);

    cen = cen.truncateAtPow(round_pow);
    err = err.truncateAtPow(round_pow);
  }

  void convertDecimal(){
    cen = cen.setExp(0);
    err = err.setExp(0);
  }

};





//Central value (error)
template<typename DistributionType, typename ValuePolicy = printStats<DistributionType> >
struct publicationDistributionPrinter: public distributionPrinter<DistributionType>{
public:
  int nsf; //number of sig figs
  SigFigsSource sfsrc; //whether the sig.figs specified is based on the error or the central value
  int min_width; //pad with trailing spaces if width < min_width. Use 0 for no padding
  int sci_fmt_threshold; //when the operative value's exponent is > this value we will use scientific format, otherwise we will use decimal format
public:
  void print(std::ostream &os_out, const DistributionType &d) const{
    std::stringstream os;    
    std::ios::streampos init_pos = os.tellp();
    typedef decltype(ValuePolicy::centralValue(d)) valueType;
    valueType mu = ValuePolicy::centralValue(d);
    valueType err = ValuePolicy::error(d);

    cenErr ce(mu,err,sfsrc,nsf);

    int abs_pow10_op;
    {
      SigFigsSource sff = sfsrc;
      if(sff == Largest) sff = (mu >= err ? Central : Error); 
      abs_pow10_op = abs(decimal::pow10(sff == Central ? mu : err));
    }
    if(abs_pow10_op <= sci_fmt_threshold)
      ce.convertDecimal();    
    
    ce.printPublication(os);

    int width = os.tellp() - init_pos;

    for(int i=0;i<min_width - width;i++) os << ' ';

    os_out << os.rdbuf();

  }

  publicationDistributionPrinter(const int _nsf = 3, const SigFigsSource _sfsrc = Largest): nsf(_nsf), sfsrc(_sfsrc), min_width(0), sci_fmt_threshold(3){}

  inline void setSigFigs(const int _nsf, const SigFigsSource _sfsrc = Largest){ nsf = _nsf; sfsrc  = _sfsrc; }
  inline void setMinWidth(const int _min_width){ min_width = _min_width; }
  inline void setSciFormatThreshold(const int p){ sci_fmt_threshold = p; }
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
