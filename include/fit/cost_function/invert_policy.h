#ifndef _CPSFIT_COST_FUNCTION_INVERT_POLICY_H_
#define _CPSFIT_COST_FUNCTION_INVERT_POLICY_H_

//User can optionally modify the procedure used by the minimizer to perform the matrix inversion

#include<config.h>
#include<tensors/numeric_square_matrix.h>

CPSFIT_START_NAMESPACE

template<typename CostType>
struct CostFunctionSVDinvert{
  inline static NumericSquareMatrix<CostType> invert(const NumericSquareMatrix<CostType> &in){
    NumericSquareMatrix<CostType> out(in.size()); svd_inverse(out,in); return out;
  }
};

CPSFIT_END_NAMESPACE
#endif
