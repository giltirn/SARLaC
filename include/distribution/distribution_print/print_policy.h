#ifndef _CPSFIT_PRINTPOLICY_H_
#define _CPSFIT_PRINTPOLICY_H_

#include<utils/macros.h>

CPSFIT_START_NAMESPACE

//Use policies to control how the central value and error are defined for a particular distribution

//The default policy: the best field and standard error on the mean
template<typename DistributionType>
struct printStats{
  inline static auto centralValue(const DistributionType &d)->decltype(d.best()){ return d.best(); }
  inline static auto error(const DistributionType &d)->decltype(d.standardError()){ return d.standardError(); }
};

CPSFIT_END_NAMESPACE
#endif
