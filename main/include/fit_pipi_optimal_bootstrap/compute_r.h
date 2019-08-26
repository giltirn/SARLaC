#ifndef FIT_PIPI_OPTIMAL_BOOTSTRAP_COMPUTE_R_H
#define FIT_PIPI_OPTIMAL_BOOTSTRAP_COMPUTE_R_H

CPSFIT_START_NAMESPACE

//idx are the indices for the couplings of this operator to each state
#define PARAM_ELEM_MEMBERS \
  (std::string, file)				\
  (std::vector<int>, idx)

struct ParamElem{
  GENERATE_MEMBERS(PARAM_ELEM_MEMBERS)
  ParamElem(): file("file.hdf5"), idx({0,1,2}){}
};
GENERATE_PARSER(ParamElem, PARAM_ELEM_MEMBERS)


std::vector<bootstrapDistributionD> computeR(const std::vector<ParamElem> &op_amplitudes, const double Ascale){
  int nop = op_amplitudes.size();
  std::cout << "Populating amplitude matrix for " << nop << " operators" << std::endl;

  //Load amplitude matrix A_ia    i=state  a=op
  NumericSquareMatrix<bootstrapDistributionD> A(nop);
  for(int a=0;a<nop;a++){
    std::cout << a << " " << op_amplitudes[a].file << std::endl;

    std::vector<bootstrapDistributionD> v;
    readParamsStandard(v, op_amplitudes[a].file);
    for(int i=0;i<nop;i++)
      A(i,a) = v[ op_amplitudes[a].idx[i] ] * sqrt(Ascale);
  }

  std::cout << "A: " << std::endl << A << std::endl;

  NumericSquareMatrix<bootstrapDistributionD> Ainv(A);
  svd_inverse(Ainv, A);
  
  std::cout << "A^-1: " << std::endl << Ainv << std::endl;

  std::vector<bootstrapDistributionD> r(nop);
  for(int a=0;a<nop;a++) r[a] = Ainv(a,0);
  
  std::cout << "r: " << r << std::endl;
  
  return r;
}

CPSFIT_END_NAMESPACE

#endif
