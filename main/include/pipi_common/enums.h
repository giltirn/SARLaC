#ifndef _FIT_PIPI_GPARITY_ENUMS_H_
#define _FIT_PIPI_GPARITY_ENUMS_H_

#include<parser.h>

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

GENERATE_ENUM_AND_PARSER(PiPiProjector, (A1)(A1momSet111)(A1momSet311)(Avg4)(Avg2)(Solo) );
GENERATE_ENUM_AND_PARSER(PiPiMomAllowed, (All)(Orig64)(ParityAxisPermSymmReduced)(AuxDiagSymmReduced)(AuxDiagParityAxisPermSymmReduced) );
GENERATE_ENUM_AND_PARSER(PiPiCorrSelector, (Basic)(SymmetrySubset) );
GENERATE_ENUM_AND_PARSER(PiPiFitFunction, (FCoshPlusConstant)(FCoshPlusConstantDoubleExp) );
GENERATE_ENUM_AND_PARSER(PiPiEffectiveEnergy, (TwoPoint)(TwoPointSubConstant)(ThreePoint) );

CPSFIT_END_NAMESPACE

#endif
