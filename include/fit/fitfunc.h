#ifndef _FITFUNC_H
#define _FITFUNC_H

/*
Some predefined fit functions for commonly encountered scenarios
The user is of course free to define their own providing they conform to the following template:

Required typedefs:
  ValueType - the return type of the fit func evaluated at some coordinate
  ParameterType - the input type representing the parameters of the fit
  ValueDerivativeType - the type containing the first derivatives wrt the fit parameters
  GeneralizedCoordinate - the type representing the coordinate at which the function is evaluated

Required methods:
  ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const   - Evaluate the function
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const - The partial derivatives of the function wrt the parameters
  int Nparams() const - The number of parameters
*/

#include<fit/fitfunc/fitfunc_correlator.h>
#include<fit/fitfunc/fitfunc_linearmultidim.h>
#include<fit/fitfunc/fitfunc_other.h>
#include<fit/fitfunc/fitfunc_frozen.h>
#include<fit/fitfunc/fitfunc_mapping.h>

#endif
