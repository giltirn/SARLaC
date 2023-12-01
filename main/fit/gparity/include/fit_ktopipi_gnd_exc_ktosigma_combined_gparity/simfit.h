#ifndef _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_SIMFIT_H
#define _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_SIMFIT_H

#include<config.h>
#include<utils/macros.h>

#include "input_param_args.h"

#include "simfit_base.h"
#include "simfit_jackknife_containers.h"
#include "simfit_bootstrap_containers.h"
#include "simfit_multistate.h"
#include "simfit_multistate_wavg.h"

SARLAC_START_NAMESPACE

template< template<typename, template<typename> class> class DistributionType > 
simultaneousFitBase<DistributionType>* getFitter(const SimFitFunction ff, const int nstate, const std::vector<PiPiOperator> &operators){
  switch(ff){
  case SimFitFunction::MultiState:
return new simultaneousFitMultiState<DistributionType>(nstate, operators);
  case SimFitFunction::MultiStateWavg:
return new simultaneousFitMultiStateWavg<DistributionType>(nstate, operators);
  default:
    assert(0);
  }
};


template<template<typename, template<typename> class> class DistributionType, typename Params, template<typename> class V>
std::vector<std::vector<DistributionType<double,V> > > convert7basisTo10basis(const int nstate,
									      const std::vector<DistributionType<Params,V> > &params){		    

  std::cout << "Converting 7 basis results to 10 basis" << std::endl;
  typedef DistributionType<Params,V> DistributionTypeP;
  typedef iterate<DistributionTypeP> iter_p;
  typedef DistributionType<double,V> DistributionTypeD;
  typedef iterate<DistributionTypeD> iter_d;
  DistributionTypeD zero(params[0].getInitializer()); //will initialize correctly jackknife or bootstrap
  zeroit(zero);

  std::vector<std::vector<DistributionTypeD> > params_10(10, std::vector<DistributionTypeD>(nstate, zero));    
    
  struct InContainer{
    size_t s;
    size_t idx;
    const std::vector<DistributionTypeP> &d;
    double operator[](const int q) const{ return iter_p::at(s, d[q])(idx); }
    InContainer(size_t s, size_t idx, const std::vector<DistributionTypeP> &d): s(s), idx(idx), d(d){}
  };
  struct OutContainer{
    size_t s;
    size_t state;
    std::vector<std::vector<DistributionTypeD> > &d;
    double & operator[](const int q){ return iter_d::at(s,d[q][state]); }
    OutContainer(size_t s, size_t state, std::vector<std::vector<DistributionTypeD> > &d): s(s), state(state), d(d){}
  };

  auto const* tag_map = params[0].sample(0).getTagMap();
  assert(tag_map != NULL);

  for(int i=0;i<nstate;i++){
    auto it = tag_map->find(stringize("M%d",i));
    assert(it != tag_map->end());
    int Midx = it->second;
      
    for(int s=0;s<iter_d::size(zero);s++){
      InContainer in(s, Midx, params);
      OutContainer out(s, i, params_10);
      convert7to10(out,in);
    }
  }
  return params_10;
}

SARLAC_END_NAMESPACE

#endif
