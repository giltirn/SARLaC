#ifndef _FIT_KTOSIGMA_KTOPIPI_GPARITY_ENUMS_H_
#define _FIT_KTOSIGMA_KTOPIPI_GPARITY_ENUMS_H_

#include<config.h>
#include<utils/macros.h>

SARLAC_START_NAMESPACE

GENERATE_ENUM_AND_PARSER(PiPiOperator, (PiPiGnd)(PiPiExc)(Sigma) );

inline void write(HDF5writer &wr, const PiPiOperator op, const std::string &tag){
  SARLaC::write(wr, (int)op, tag);
}
inline void read(HDF5reader &rd, PiPiOperator &op, const std::string &tag){
  int r;
  SARLaC::read(rd, r, tag);
  op = (PiPiOperator)r;
}

inline std::string opAmplitudeParamFmt(PiPiOperator op){
  switch(op){
  case PiPiOperator::PiPiGnd:
    return "Apipi%d";
  case PiPiOperator::PiPiExc:
    return "Apipi_exc_%d";
  case PiPiOperator::Sigma:
    return "Asigma%d";
  }
  assert(0);
  return "";
}
inline std::string opDescr(PiPiOperator op){
  switch(op){
  case PiPiOperator::PiPiGnd:
    return "K->pipi(111)";
  case PiPiOperator::PiPiExc:
    return "K->pipi(311)";
  case PiPiOperator::Sigma:
    return "K->sigma";
  }
  assert(0);
  return "";
}

inline std::string opDescrFile(PiPiOperator op){
  switch(op){
  case PiPiOperator::PiPiGnd:
    return "kpipi_111";
  case PiPiOperator::PiPiExc:
    return "kpipi_311";
  case PiPiOperator::Sigma:
    return "ksigma";
  }
  assert(0);
  return "";
}



GENERATE_ENUM_AND_PARSER(SimFitFunction, (MultiState)(MultiStateWavg) );
GENERATE_ENUM_AND_PARSER(Basis, (Basis10)(Basis7) ); //choose the basis for the operators of the Weak Hamiltonian
GENERATE_ENUM_AND_PARSER(CovarianceMatrix, (Regular)(Block) );

SARLAC_END_NAMESPACE

#endif
