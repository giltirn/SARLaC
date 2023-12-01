#ifndef _FIT_PIPI_GPARITY_ENUMS_H_
#define _FIT_PIPI_GPARITY_ENUMS_H_

#include<parser.h>

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

GENERATE_ENUM_AND_PARSER(PiPiFitFunction, (FCoshPlusConstant)(FCoshPlusConstantDoubleExp) );
GENERATE_ENUM_AND_PARSER(PiPiEffectiveEnergy, (TwoPoint)(TwoPointSubConstant)(ThreePoint) );
GENERATE_ENUM_AND_PARSER(MomentumUnit, (PiOverL)(PiOverTwoL) );

CPSFIT_END_NAMESPACE

#endif
