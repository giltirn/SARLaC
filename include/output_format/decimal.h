#ifndef _CPSFIT_DECIMAL_H_
#define _CPSFIT_DECIMAL_H_

#include<cstdlib>
#include<cmath>
#include<cassert>
#include<vector>
#include<sstream>
#include<iostream>
#include<deque>

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

//A class that stores a number in a decimal representation for output formatting and truncation
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

  std::ostream & print(std::ostream &os, bool incl_exponent = true, bool leading_space_if_positive = false) const{
    if(leading_space_if_positive) os << ( sgn == -1 ? '-' : ' ' );
    else if(sgn==-1) os << '-';
    
    for(int i=0;i<v.size();i++){
      os << v[i];
      if(i==dp && i!=v.size()-1) os << '.';
    }
    if(incl_exponent && base_pow !=0 ) os << 'e' << base_pow;
    return os;
  }

  //For debugging, print the full state
  std::ostream & report(std::ostream &os){
    os << "dp="<<dp << " base_pow=" << base_pow << " sgn="<< sgn << " v={";
    for(int i=0;i<v.size();i++){ os << v[i]; if(i==dp) os << '|'; }
    os << "}";
  }

  double value() const{
    std::stringstream ss;
    print(ss);
    double o;
    ss >> o;
    return o;
  }

  //Vector access of individual digits
  inline int nDigits() const{ return v.size(); }
  inline int operator[](const int i) const{ return v[i]; }
  
  int exponent() const{ return base_pow; }

  //Right-shift = divide by 10  (keep base power the same)
  decimal operator>>(const unsigned int n) const{
    decimal out(*this);
    for(unsigned int i=0;i<n;i++){
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
      if(out.dp == out.v.size())
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

  //Remove all digits after power p)
  decimal truncateAtPow(const int p) const{
    int orig_pow = exponent();
    decimal out = this->setExp(p);
    out.v.resize(out.dp+1); //just keep all ints up to the decimal point
    out = out.setExp(orig_pow);
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

CPSFIT_END_NAMESPACE

#endif
