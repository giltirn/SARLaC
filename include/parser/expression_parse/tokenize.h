#ifndef _EXPRESSION_PARSE_TOKENIZE_H_
#define _EXPRESSION_PARSE_TOKENIZE_H_

#include<list>

#include<config.h>
#include<utils/macros.h>
#include<parser/expression_parse/AST.h>

CPSFIT_START_NAMESPACE

//A tokenize for math expressions
namespace _mathExpressionTokenize{
  namespace ascii = boost::spirit::x3::ascii;
  namespace x3 = boost::spirit::x3;
  
  auto addtoken_ = [&](auto& ctx){
    std::ostringstream os; os << x3::_attr(ctx);
    x3::_val(ctx).push_back(os.str());
  };
  auto const string_rule = x3::lexeme[+(x3::alnum)];

  x3::rule<struct parse_math_, std::list<std::string> > const parse_math_rule = "parse_math_rule";
  auto const parse_math_rule_def = x3::eps >> *(x3::char_("+")[addtoken_] | x3::char_("-")[addtoken_] | x3::char_("*/^()")[addtoken_] | x3::double_[addtoken_] | string_rule[addtoken_]);
  BOOST_SPIRIT_DEFINE(parse_math_rule);
};
  
std::list<std::string> mathExpressionTokenize(const std::string &s){
  using namespace _mathExpressionTokenize;
  std::list<std::string> tokens;
  bool r = x3::phrase_parse(s.begin(), s.end(), parse_math_rule, ascii::space, tokens);      
  if(!r){ std::cout << "Could not properly add spaces around tokens in string \"" << s << "\"\n"; exit(-1); }
  return tokens;
}

CPSFIT_END_NAMESPACE
#endif
