#ifndef _FIT_KTOPIPI_GPARITY_ENUMS_H_
#define _FIT_KTOPIPI_GPARITY_ENUMS_H_


#include<config.h>
#include<utils/macros.h>

#include<parser.h>

CPSFIT_START_NAMESPACE

GENERATE_ENUM_AND_PARSER(KtoPiPiFitFunc, (FitSeparate)(FitSimultaneous)(FitSimultaneousChiralBasis)(FitSeparateWithConstant)(FitSeparateTwoExp)(FitSeparateTwoExpKaon)(FitSeparateExcPiK) );

CPSFIT_END_NAMESPACE

#endif
