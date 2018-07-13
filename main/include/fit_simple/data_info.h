#ifndef _FIT_SIMPLE_DATAINFO_H_
#define _FIT_SIMPLE_DATAINFO_H_

GENERATE_ENUM_AND_PARSER(ParserType, (ParserStandard)(ParserMultiSourceAverage)(ParserMultiSourceAverageImag) );
GENERATE_ENUM_AND_PARSER(TimeDependence, (TimeDepNormal)(TimeDepReflect)(TimeDepFold)(TimeDepAntiFold) );
GENERATE_ENUM_AND_PARSER(Combination, (CombinationAverage)(CombinationAminusB)(CombinationAdivB) );
GENERATE_ENUM_AND_PARSER(FitFuncType, (FCosh)(FSinh)(FExp)(FConstant)(FTwoStateCosh) );

#define DATA_INFO_MEMBERS \
  ( ParserType, parser )	       \
  ( std::string, operation )   \
  ( TimeDependence, time_dep ) \
  ( std::string, file_fmt )

//file_fmt should contain a '%d' which is replaced by the trajectory index
//operation is a math expression. Use x to represent the data. Use an empty string to leave as-is

struct DataInfo{
  GENERATE_MEMBERS(DATA_INFO_MEMBERS);
  DataInfo(): parser(ParserType::ParserStandard), operation(""), time_dep(TimeDependence::TimeDepNormal), file_fmt("data.%d"){}
};
GENERATE_PARSER(DataInfo, DATA_INFO_MEMBERS);

#endif
