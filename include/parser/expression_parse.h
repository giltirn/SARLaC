#ifndef _EXPRESSION_PARSE_H_
#define _EXPRESSION_PARSE_H_

#include<config.h>
#include<utils/macros.h>
#include<parser/expression_parse/AST.h>
#include<parser/expression_parse/tokenize.h>
#include<parser/expression_parse/shunting_yard.h>

SARLAC_START_NAMESPACE

inline expressionAST mathExpressionParse(const std::string &s){
  shuntingYardParser parser;
  return parser.parse(s);
};

SARLAC_END_NAMESPACE

#endif
