#ifndef _FIT_KTOPIPI_COMPUTE_Q_AMPLITUDE_OPTS
#define _FIT_KTOPIPI_COMPUTE_Q_AMPLITUDE_OPTS

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

struct computeQamplitudeOpts{
  double alpha_scale; //multiply the coefficient alpha by some factor

  computeQamplitudeOpts(): alpha_scale(1.){}
};

CPSFIT_END_NAMESPACE

#endif
