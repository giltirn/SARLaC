#pragma once

GENERATE_ENUM_AND_PARSER(DataGenStrategy, (NormalUniform)(NormalTimeDep)(LogNormalUniform) );

GENERATE_ENUM_AND_PARSER(CovMatStrategy, (Correlated)(Uncorrelated)(Cutoff) );

GENERATE_ENUM_AND_PARSER(FitFuncType, (FConstant));

GENERATE_ENUM_AND_PARSER(preAnalysisType, (None)(CovMatEvals));
