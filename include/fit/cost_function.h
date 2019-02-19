#ifndef _COST_FUNCTION_H
#define _COST_FUNCTION_H

/*Cost functions for minimization

  Note regarding need for only first derivative:

  The Marquardt-Levenberg algorithm is a combination of the Gauss-Newton method and gradient descent.
  Gradient descent updates along the direction - grad(f) where f is the fitness function, and thus doesnt need the second derivatives.
  The Gauss-Newton method taylor expands the gradient vector
  grad f(x) = grad f(x0) + (x-x0)^T grad^2 f(x0) + higher order
  and solve for minimum grad f(x) = 0 at some x, then  x_i+1 = x_i - [grad^2 f(x_i)]^-1 grad f(x_i)
  
  as it cuts the grad off at first order, this is equivalent to approximating the function as a quadratic. 

  As  chisq = sum_i ( Y_i - fit( {X}_i ) )^2/sigma_i^2
  then approximating this as quadratic in the parameters {X} means that the fit function is approximated as linear.
  
  As a result we can treat the second derivatives of the fit function as zero without losing anything.
  
  We thus do not bother going to the effort of determining the second derivatives.
  cf. http://www.scribd.com/doc/10093320/Levenberg-Marquardt-Algorithm  (page 2)
*/

#include<fit/cost_function/invert_policy.h>
#include<fit/cost_function/uncorrelated_chisq.h>
#include<fit/cost_function/correlated_chisq.h>
#include<fit/cost_function/correlated_cov_chisq.h>
#include<fit/cost_function/correlated_chisq_terms.h>
#include<fit/cost_function/uncorrelated_chisq_terms.h>

#endif
