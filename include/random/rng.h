#ifndef _RANDOM_RNG_H_
#define _RANDOM_RNG_H_

//A global random number generator using c++11 Mersenne twister

#include <random>

#include<config.h>
#include<utils/macros.h>
#include<omp.h>

CPSFIT_START_NAMESPACE

class RNGstore{
public:
  typedef std::mt19937 RNGtype;
  typedef typename RNGtype::result_type seedType;
private:
  RNGtype* rng;
public:
  RNGstore(): rng(NULL){}
  RNGstore(const RNGstore &r): rng( r.rng == NULL ? NULL : new RNGtype(*r.rng) ){}
  RNGstore(RNGstore &&r): rng(r.rng){ r.rng = NULL; }  

  RNGstore & operator=(const RNGstore &r){
    if(rng != NULL) delete rng;
    rng = r.rng == NULL ? NULL : new RNGtype(*r.rng);
    return *this;
  }
  RNGstore & operator=(RNGstore &&r){
    if(rng != NULL) delete rng;
    rng = r.rng;
    r.rng = NULL;
    return *this;
  }

  void initialize(const seedType seed){ if(rng==NULL) rng = new RNGtype(seed); }
  void initialize(){ if(rng == NULL) rng = new RNGtype(); }
  
  RNGstore(const seedType seed): rng(NULL){ initialize(seed); }

  bool isInitialized() const{ return rng != NULL; }
  
  RNGtype &operator()(){ return *rng; }

  ~RNGstore(){ if(rng!=NULL) delete rng; }
};

RNGstore RNG; //static instance

struct threadRNGstore{
  std::vector<RNGstore> rngs;
  typedef RNGstore::seedType seedType;

  threadRNGstore(){}

  void initialize(const seedType seed, const int nthr = omp_get_max_threads()){
    rngs.resize(nthr);
    for(int i=0;i<rngs.size();i++)
      rngs[i].initialize( seed + 1 + i );
  }
    
  RNGstore &operator()(const int thr = omp_get_thread_num()){ return rngs[thr]; }
};

threadRNGstore threadRNG;


CPSFIT_END_NAMESPACE
#endif
