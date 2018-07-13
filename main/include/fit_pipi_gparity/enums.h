#ifndef _FIT_PIPI_GPARITY_ENUMS_H_
#define _FIT_PIPI_GPARITY_ENUMS_H_

GENERATE_ENUM_AND_PARSER(PiPiProjector, (A1)(Avg4)(Avg2)(Solo) );
GENERATE_ENUM_AND_PARSER(PiPiMomAllowed, (All)(Orig64)(ParityAxisPermSymmReduced)(AuxDiagSymmReduced)(AuxDiagParityAxisPermSymmReduced) );
GENERATE_ENUM_AND_PARSER(PiPiCorrSelector, (Basic)(SymmetrySubset) );
GENERATE_ENUM_AND_PARSER(PiPiFitFunction, (FCoshPlusConstant)(FCoshPlusConstantDoubleExp) );
GENERATE_ENUM_AND_PARSER(PiPiEffectiveEnergy, (TwoPoint)(TwoPointSubConstant)(ThreePoint) );


#endif
