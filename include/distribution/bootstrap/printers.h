#ifndef _BOOTSTRAP_PRINTERS_H_
#define _BOOTSTRAP_PRINTERS_H_


#include<utils/macros.h>
#include<distribution/bootstrap/class.h>
#include<distribution/distribution_print/print_policy.h>

CPSFIT_START_NAMESPACE

template<typename T>
struct bootstrapAsymmetricErrorPrinter{};

template<typename T, template<typename> class U>
struct bootstrapAsymmetricErrorPrinter<bootstrapDistribution<T,U> >: public distributionPrinter<bootstrapDistribution<T,U> >{
  void print(std::ostream &os, const bootstrapDistribution<T,U> &dist) const{
    auto lh = dist.errorBounds();
    os << "(" << dist.best() << " + " << lh.second << " - " << lh.first << ")";
  }
};

CPSFIT_END_NAMESPACE
#endif
