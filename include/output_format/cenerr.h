#ifndef _SARLAC_CENERR_H_
#define _SARLAC_CENERR_H_

//A class that contains a central value and error in decimal representation for printing purposes

#include<config.h>
#include<utils/macros.h>
#include<output_format/decimal.h>

SARLAC_START_NAMESPACE

enum SigFigsSource { Central, Error, Largest };

class cenErr{
  decimal cen;
  decimal err;
public:
  
  const decimal & getCen() const{ return cen; }
  const decimal & getErr() const{ return err; }

  std::ostream & printBasic(std::ostream &os) const{
    os << cen << " +/- " << err; return os;
  }
  std::ostream & printPublication(std::ostream &os) const{
    assert(err.exponent() == cen.exponent());

    //In situations where we have numbers like  0.000(15) we need to make sure all the zeroes are printed or it will look like   0(15)
    bool print_all_digits_if_zero = err.firstSigFigPow() < 0;
    cen.print(os, false,true, print_all_digits_if_zero);
    os << "(";
    err.printWithoutLeadingZeros(os);
    os << ")";
    if(cen.exponent() != 0) os << "\\times 10^{" << cen.exponent() << "}";
    
    return os;
  }
  
  cenErr(double c, double e): cen(c), err(e){
    assert(e>=0.);    
    //Put them under a common exponent
    int exp = cen.exponent(); // < err.exponent() ? cen.exponent() : err.exponent();
    cen = cen.setExp(exp);
    err = err.setExp(exp);
  }

  void setSigFigs(const int nsf, SigFigsSource sf = Largest){
    assert(err.exponent() == cen.exponent());

    if(sf == Largest)
      sf = fabs(cen.value()) >= fabs(err.value()) ? Central : Error;
    
    int sig_fig_pow_cen = cen.firstSigFigPow();
    int sig_fig_pow_err = err.firstSigFigPow();

    int round_pow = (sf == Central ? sig_fig_pow_cen : sig_fig_pow_err) - nsf + 1;

    //Round and truncate
    cen = cen.roundAtPow(round_pow);
    err = err.roundAtPow(round_pow);

    cen = cen.truncateAtPow(round_pow);
    err = err.truncateAtPow(round_pow);
  }
  void setRoundPower(const int round_pow){
    assert(err.exponent() == cen.exponent());

    cen = cen.roundAtPow(round_pow);
    err = err.roundAtPow(round_pow);

    cen = cen.truncateAtPow(round_pow);
    err = err.truncateAtPow(round_pow);
  }  
  
  void convertDecimal(){
    cen = cen.setExp(0);
    err = err.setExp(0);
  }
  void setExp(int p){
    cen = cen.setExp(p);
    err = err.setExp(p);
  }
};


SARLAC_END_NAMESPACE

#endif
