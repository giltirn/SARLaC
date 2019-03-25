#ifndef _PIPI_DATA_FILTERS_H
#define _PIPI_DATA_FILTERS_H

#include "enums.h"

GENERATE_ENUM_AND_PARSER(DataFilter, (TimeGreaterThan)(RelErrGreaterThan)(TimeOdd)(TimeEven) );

#define FILTER_MEMBERS \
  ( bool, all_ops)     \
  ( Operator, op1 )    \
  ( Operator, op2 )    \
  ( DataFilter, filter )			\
  ( std::string, value )

struct Filter{
  GENERATE_MEMBERS(FILTER_MEMBERS);

  bool filterOut(const Operator dop1, const Operator dop2, const int t, const jackknifeDistribution<double> &dval, std::string *reason = NULL){
    if(all_ops || (dop1 == op1 && dop2 == op2)){
      if(filter == DataFilter::TimeGreaterThan && t > strToAny<int>(value)){
	if(reason) *reason = stringize("Time %d > %s", t , value.c_str());
	return true;
      }
      else if(filter == DataFilter::RelErrGreaterThan && fabs(dval.standardError()/dval.mean()) > strToAny<double>(value)){
	if(reason) *reason = stringize("Rel.Err %f > %s", fabs(dval.standardError()/dval.mean()) , value.c_str());
	return true;
      }
      else if(filter == DataFilter::TimeOdd && t % 2 == 1){
	if(reason) *reason = stringize("Time %d is odd", t);
	return true;
      }
      else if(filter == DataFilter::TimeEven && t % 2 == 0){
	if(reason) *reason = stringize("Time %d is even", t);
	return true;
      }
    }
    return false;
  }

  Filter(): all_ops(false), op1(Operator::PiPiGnd), op2(Operator::PiPiExc), filter(DataFilter::TimeGreaterThan), value("12"){}
};

GENERATE_PARSER(Filter, FILTER_MEMBERS);

#define FILTERS_MEMBERS \
  ( std::vector<Filter>, filters )

struct Filters{
  GENERATE_MEMBERS(FILTERS_MEMBERS);

  Filters(): filters(1){}
};

GENERATE_PARSER(Filters, FILTERS_MEMBERS);


#endif
