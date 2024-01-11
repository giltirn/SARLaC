#pragma once

GENERATE_ENUM_AND_PARSER(DataGenStrategy, (NormalUniform)(NormalTimeDep)(LogNormalUniform)(NormalUniformPlusShift)(NormalUniformMixLeft)(NormalTimeDepMixLeft)(Binned)(NormalUniformMetropolis) );

GENERATE_ENUM_AND_PARSER(CovMatStrategy, (Correlated)(Uncorrelated)(Cutoff)(MCM) );

GENERATE_ENUM_AND_PARSER(FitFuncType, (FConstant)(FLinear)(FConstantFrozen));

GENERATE_ENUM_AND_PARSER(preAnalysisType, (None)(CovMatEvals)(CorrMatEvals)(StandardError)(FitAutoCorrAvoid)(TauInt)(BlockBootstrapMeanBias)(BlockBootstrapQ2Bias)(BlockBootstrapStdErrBias));
