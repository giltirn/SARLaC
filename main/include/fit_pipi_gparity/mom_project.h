#ifndef _PIPI_MOM_PROJECT_H_
#define _PIPI_MOM_PROJECT_H_

struct PiPiProject{
  virtual bool operator()(std::complex<double> &coeff, const threeMomentum &mom) const = 0;
};
struct PiPiProjectA1: public PiPiProject{
  virtual bool operator()(std::complex<double> &coeff, const threeMomentum &mom) const{
    coeff = 1./8;
    return true;
  }
};

PiPiProject* getProjector(const PiPiProjector p){
  switch(p){
  case A1:
    return (PiPiProject*)(new PiPiProjectA1);
  default:
    error_exit(std::cout << "getProjector unknown projector " << p << std::endl);
  }
}

#endif
