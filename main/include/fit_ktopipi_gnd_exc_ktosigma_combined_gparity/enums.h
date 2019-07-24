#ifndef _FIT_KTOSIGMA_KTOPIPI_GPARITY_ENUMS_H_
#define _FIT_KTOSIGMA_KTOPIPI_GPARITY_ENUMS_H_

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

GENERATE_ENUM_AND_PARSER(PiPiOperator, (PiPiGnd)(PiPiExc)(Sigma) );

inline void write(HDF5writer &wr, const PiPiOperator op, const std::string &tag){
  CPSfit::write(wr, (int)op, tag);
}
inline void read(HDF5reader &rd, PiPiOperator &op, const std::string &tag){
  int r;
  CPSfit::read(rd, r, tag);
  op = (PiPiOperator)r;
}


GENERATE_ENUM_AND_PARSER(SimFitFunction, (MultiState)(MultiStateWavg) );
GENERATE_ENUM_AND_PARSER(Basis, (Basis10)(Basis7) ); //choose the basis for the operators of the Weak Hamiltonian
GENERATE_ENUM_AND_PARSER(CovarianceMatrix, (Regular)(Block) );

CPSFIT_END_NAMESPACE

#endif
