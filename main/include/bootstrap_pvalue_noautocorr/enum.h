#pragma once

GENERATE_ENUM_AND_PARSER(DataGenStrategy, (NormalUniform)(NormalTimeDep)(LogNormalUniform)(NormalUniformPlusShift)(NormalUniformMixLeft)(Binned) );

GENERATE_ENUM_AND_PARSER(CovMatStrategy, (Correlated)(Uncorrelated)(Cutoff)(MCM) );

GENERATE_ENUM_AND_PARSER(FitFuncType, (FConstant));

GENERATE_ENUM_AND_PARSER(preAnalysisType, (None)(CovMatEvals)(CorrMatEvals)(StandardError));
