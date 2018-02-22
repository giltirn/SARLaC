#ifndef _RANDOM_RNG_H_
#define _RANDOM_RNG_H_

//A global random number generator using c++11 Mersenne twister

#include <random>

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

class RNGstore{
public:
  typedef std::mt19937 RNGtype;
  typedef typename RNGtype::result_type seedType;
private:
  RNGtype* rng;
public:
  RNGstore(): rng(NULL){}

  void initialize(const seedType seed){ if(rng==NULL) rng = new RNGtype(seed); }
  void initialize(){ if(rng == NULL) rng = new RNGtype(); }
  
  RNGstore(const seedType seed){ initialize(seed); }

  bool isInitialized() const{ return rng != NULL; }
  
  RNGtype &operator()(){ return *rng; }
};

RNGstore RNG; //static instance

CPSFIT_END_NAMESPACE
#endif
