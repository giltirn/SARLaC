#ifndef FIT_PIPI_GND_EXC_GPARITY_H_
#define FIT_PIPI_GND_EXC_GPARITY_H_

GENERATE_ENUM_AND_PARSER(FitFuncType, (FSimGenOneState)(FSimGenTwoState)(FSimGenThreeState)(FSimGenMultiState)(FSimGenMultiStateCparam)(FSimGenThreeStateLogEdiff)(FSimGenMultiStateLogEdiff)(FSimGenMultiStateSub)(FSimGenMultiStateTminSub)(FSimGenMultiStateTminSubForceZero) );

GENERATE_ENUM_AND_PARSER(CovarianceMatrix, (Regular)(Frozen)(BlockHybrid)(Block) );
#endif
