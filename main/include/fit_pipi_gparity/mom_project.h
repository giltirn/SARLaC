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

//An extra discriminator for momentum pairs
struct PiPiMomAllow{
  virtual bool operator()(double &coeff, const threeMomentum &p1src, const threeMomentum &p1snk) const = 0;
};
struct PiPiMomAllowAll: public PiPiMomAllow{
  bool operator()(double &coeff, const threeMomentum &p1src, const threeMomentum &p1snk) const{
    coeff = 1.; return true;
  }
};
struct PiPiMomAllowSymmReduced: public PiPiMomAllow{
  typedef std::pair<threeMomentum, threeMomentum> momPair;
  std::map<momPair, int> allow;

  PiPiMomAllowSymmReduced(){
    //Reduced combination set obtained by utilizing parity and axis interchange symmetries
    allow[ momPair(threeMomentum({-1,-1,-1}), threeMomentum({-1,-1,-1}) ) ] = 2;
    allow[ momPair(threeMomentum({-1,-1,-1}), threeMomentum({-1,-1,1}) ) ] = 6;
    allow[ momPair(threeMomentum({-1,-1,-1}), threeMomentum({-1,1,1}) ) ] = 6;
    allow[ momPair(threeMomentum({-1,-1,-1}), threeMomentum({1,1,1}) ) ] = 2;
    allow[ momPair(threeMomentum({-1,-1,1}), threeMomentum({-1,-1,-1}) ) ] = 6;
    allow[ momPair(threeMomentum({-1,-1,1}), threeMomentum({-1,-1,1}) ) ] = 6;
    allow[ momPair(threeMomentum({-1,-1,1}), threeMomentum({-1,1,-1}) ) ] = 12;
    allow[ momPair(threeMomentum({-1,-1,1}), threeMomentum({-1,1,1}) ) ] = 12;
    allow[ momPair(threeMomentum({-1,-1,1}), threeMomentum({1,1,-1}) ) ] = 6;
    allow[ momPair(threeMomentum({-1,-1,1}), threeMomentum({1,1,1}) ) ] = 6;
  }    

  bool operator()(double &coeff, const threeMomentum &p1src, const threeMomentum &p1snk) const{
    std::map<momPair, int>::const_iterator it = allow.find(momPair(p1src,p1snk));
    if(it == allow.end()) return false;
    else{
      coeff = it->second;
      return true;
    }
  }
};



PiPiMomAllow* getMomPairFilter(const PiPiMomAllowed p){
  switch(p){
  case All:
    return (PiPiMomAllow*)(new PiPiMomAllowAll);
  case ReducedSymm:
    return (PiPiMomAllow*)(new PiPiMomAllowSymmReduced);
  default:
    error_exit(std::cout << "getMomPairFilter unknown filter " << p << std::endl);
  }
}


#endif
