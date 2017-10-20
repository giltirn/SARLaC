#ifndef _EXPRESSION_PARSE_H_
#define _EXPRESSION_PARSE_H_

#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <cassert>
#include <sstream>
//#define BOOST_SPIRIT_X3_DEBUG
#include<string>
#include<cstdio>
#include<iostream>
#include<list>
#include <boost/spirit/home/x3.hpp>

//An AST representing a mathematical expression
class expressionAST{
  struct Leaf{
    virtual double value() const = 0;
    virtual ~Leaf(){}
  };
  struct UnOp: public Leaf{
    Leaf* a;
    ~UnOp(){ delete a; }
  };
  struct UnOpNegative: public UnOp{
    inline double value() const{ return -this->a->value(); }
  };
  struct UnOpSin: public UnOp{
    inline double value() const{ return sin(this->a->value()); }
  };
  struct UnOpCos: public UnOp{
    inline double value() const{ return cos(this->a->value()); }
  };
  struct UnOpTan: public UnOp{
    inline double value() const{ return tan(this->a->value()); }
  };
  struct UnOpSqrt: public UnOp{
    inline double value() const{ return sqrt(this->a->value()); }
  };
  struct UnOpExp: public UnOp{
    inline double value() const{ return exp(this->a->value()); }
  };
  struct UnOpLog: public UnOp{
    inline double value() const{ return log(this->a->value()); }
  };
  struct UnOpLog10: public UnOp{
    inline double value() const{ return log10(this->a->value()); }
  };

  struct BinOp: public Leaf{
    Leaf* a;
    Leaf* b;
    ~BinOp(){ delete a; delete b; }
  };

  struct BinOpPlus: public BinOp{
    inline double value() const{ return this->a->value() + this->b->value(); }
  };
  struct BinOpMinus: public BinOp{
    inline double value() const{ return this->a->value() - this->b->value(); }
  };
  struct BinOpTimes: public BinOp{
    inline double value() const{ return this->a->value() * this->b->value(); }
  };
  struct BinOpPow: public BinOp{
    inline double value() const{ return pow(this->a->value(),this->b->value()); }
  };
  struct Number: public Leaf{
    double v;
    Number(const std::string s){
      std::stringstream ss; ss<<s; ss>>v;
    }
    
    inline double value() const{ return v; }
  };
  struct Symbol: public Leaf{
    std::string sym;
    double v;
    bool value_set;
    
    Symbol(const std::string s): sym(s), value_set(false){
    }
    inline void setValue(const double _v){
      v = _v; value_set = true;
    }    
    inline double value() const{
      if(!value_set){ std::cout << "Error: No value has been assigned for symbol " << sym << std::endl; }
      return v;
    }
  };
  struct SymbolRedirect: public Leaf{
    Symbol* sym;
    SymbolRedirect(Symbol* _sym): sym(_sym){}
    inline double value() const{ return sym->value(); }
  };
    
  
  std::vector<Leaf*> vstack;
  std::map<std::string, Symbol*> symbols;
public:
  ~expressionAST(){
    for(int i=0;i<vstack.size();i++) delete vstack[i];
  }
  bool assignSymbol(const std::string &sym, const double val){
    std::map<std::string, Symbol*>::const_iterator it;
    if( (it=symbols.find(sym)) != symbols.end() ){
      it->second->setValue(val);
      return true;
    }else{
      return false;
    }
  }
  inline bool containsSymbol(const std::string &sym) const{
    return symbols.find(sym) != symbols.end();
  }
  inline int nSymbols() const{ return symbols.size(); }
  
  void stackOperand(const std::string &token){
    namespace ascii = boost::spirit::x3::ascii;
    namespace x3 = boost::spirit::x3;
    bool r = x3::phrase_parse(token.begin(), token.end(), x3::double_, ascii::space);      
    if(r){
      vstack.push_back(new Number(token));
    }else{
      std::map<std::string, Symbol*>::const_iterator it;
      if( (it = symbols.find(token)) != symbols.end()){
	vstack.push_back(new SymbolRedirect(it->second));
      }else{
	Symbol* s = new Symbol(token);
	vstack.push_back(s);
	symbols[token] = s;
      }	
    }      
  }
  void stackOperator(const std::string &op){
    BinOp* binop = NULL;
    UnOp* unop = NULL;
  
    if(op == "+"){
      binop = new BinOpPlus;
    }else if(op == "-"){
      binop = new BinOpMinus;
    }else if(op == "*"){
      binop = new BinOpTimes;
    }else if(op == "^"){
      binop = new BinOpPow;
    }else if(op == "neg"){
      unop = new UnOpNegative;
    }else if(op == "sin"){
      unop = new UnOpSin;
    }else if(op == "cos"){
      unop = new UnOpCos;
    }else if(op == "tan"){
      unop = new UnOpTan;
    }else if(op == "sqrt"){
      unop = new UnOpSqrt;
    }else if(op == "exp"){
      unop = new UnOpExp;
    }else if(op == "log"){
      unop = new UnOpLog;
    }else if(op == "log10"){
      unop = new UnOpLog10; 
    }else{
      std::cout << "stackOperator: Unknown operator " << op << std::endl;
      exit(-1);
    }

    if(binop){
      assert(vstack.size() >= 2);
      binop->b = vstack.back(); vstack.pop_back();
      binop->a = vstack.back(); vstack.pop_back();
      vstack.push_back(binop);
    }else if(unop){
      assert(vstack.size() >= 1);
      unop->a = vstack.back(); vstack.pop_back();
      vstack.push_back(unop);
    }  
  }

  double value() const{
    assert(vstack.size() == 1);
    return vstack[0]->value();
  }
  
};



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



//An implementation of the shunting-yard algorithm for parsing math expressions
class shuntingYardParser{
  //Shunting-yard algorithm
  // while there are tokens to be read:
  // 	read a token.
  // 	if the token is a number, then push it to the output queue.
  // 	if the token is an operator, then:
  // 		while there is an operator at the top of the operator stack with
  // 			greater than or equal to precedence and the operator is left associative:
  // 				pop operators from the operator stack, onto the output queue.
  // 		push the read operator onto the operator stack.
  // 	if the token is a left bracket (i.e. "("), then:
  // 		push it onto the operator stack.
  // 	if the token is a right bracket (i.e. ")"), then:
  // 		while the operator at the top of the operator stack is not a left bracket:
  // 			pop operators from the operator stack onto the output queue.
  // 		pop the left bracket from the stack.
  // 		/* if the stack runs out without finding a left bracket, then there are
  // 		mismatched parentheses. */
  // if there are no more tokens to read:
  // 	while there are still operator tokens on the stack:
  // 		/* if the operator token on the top of the stack is a bracket, then
  // 		there are mismatched parentheses. */
  // 		pop the operator onto the output queue.
  // exit.


  // Precedence:
  // ^ 	4 	Right
  // × 	3 	Left
  // ÷ 	3 	Left
  // + 	2 	Left
  // − 	2 	Left


  //Alternative summary of the Rules
  //     If the incoming symbols is an operand, print it..
  //     If the incoming symbol is a left parenthesis, push it on the stack.
  //     If the incoming symbol is a right parenthesis: discard the right parenthesis, pop and print the stack symbols until you see a left parenthesis. Pop the left parenthesis and discard it.
  //     If the incoming symbol is an operator and the stack is empty or contains a left parenthesis on top, push the incoming operator onto the stack.
  //     If the incoming symbol is an operator and has either higher precedence than the operator on the top of the stack, or has the same precedence as the operator on the top of the stack and is right associative -- push it on the stack.
  //     If the incoming symbol is an operator and has either lower precedence than the operator on the top of the stack, or has the same precedence as the operator on the top of the stack and is left associative -- continue to pop the stack until this is not true. Then, push the incoming operator.
  //     At the end of the expression, pop and print all operators on the stack. (No parentheses should remain.)

  typedef std::map<std::string, std::pair<int,int> > operatorMapType;
  typedef typename operatorMapType::const_iterator map_const_iterator;

  enum {Left=0, Right=1, NotAnOperator=-1};
  
  const std::map<std::string, std::pair<int,int> > & getMap() const{
    static bool initted = false;
    static std::map<std::string, std::pair<int,int> > operators;    //operator -> {precedence, associativity}
    if(!initted){
      operators["+"] = std::pair<int,int>(2,Left);
      operators["-"] = std::pair<int,int>(2,Left);
      operators["*"] = std::pair<int,int>(3,Left);
      operators["/"] = std::pair<int,int>(3,Left);
      operators["^"] = std::pair<int,int>(4,Right);
      operators["neg"] = std::pair<int,int>(4,Right); //unary negation
      operators["("] = std::pair<int,int>(NotAnOperator,NotAnOperator); //not an operator but does get pushed onto stack at times
      operators["sin"] = std::pair<int,int>(NotAnOperator,NotAnOperator);
      operators["cos"] = std::pair<int,int>(NotAnOperator,NotAnOperator);
      operators["tan"] = std::pair<int,int>(NotAnOperator,NotAnOperator);
      operators["sqrt"] = std::pair<int,int>(NotAnOperator,NotAnOperator);
      operators["exp"] = std::pair<int,int>(NotAnOperator,NotAnOperator);
      operators["log"] = std::pair<int,int>(NotAnOperator,NotAnOperator);
      operators["log10"] = std::pair<int,int>(NotAnOperator,NotAnOperator);
      initted = true;
    }
    return operators;
  }

  static inline const std::string & Op(std::map<std::string, std::pair<int,int> >::const_iterator mit){ return mit->first; }
  static inline const int Prec(std::map<std::string, std::pair<int,int> >::const_iterator mit){ return mit->second.first; }
  static inline const int Assoc(std::map<std::string, std::pair<int,int> >::const_iterator mit){ return mit->second.second; }
  static inline const bool isFunction(std::map<std::string, std::pair<int,int> >::const_iterator mit){ return (Prec(mit) == -1 && Op(mit) != "("); }

  std::vector<map_const_iterator> stack;

  void handleUnaryMinus(map_const_iterator &opit, std::string &token, const std::string &prev_token, const std::list<std::string> &tok, std::list<std::string>::const_iterator tokit){
    //Handle binary/unary minus
    //A minus sign is always binary if it immediately follows an operand or a right parenthesis, and it is always unary if it immediately follows another operator
    //or a left parenthesis, or if it occurs at the very beginning of the input. The algorithm must be modified in order to distinguish between the two.
    bool is_unary = false;
    if(tokit == tok.begin()){
      is_unary = true;	  
    }else{
      //Check if right parenthesis or operand came previously
      if( prev_token == "(" ) is_unary = true;
      else{
	map_const_iterator popit = getMap().find(prev_token);
	if(popit != getMap().end() && Prec(popit) != NotAnOperator) is_unary = true;
      }
    }
    if(is_unary){
      token = "neg";
      opit = getMap().find(token);
    }
  }
  
  inline void actionNotAnOperator(map_const_iterator opit, std::list<std::string>::const_iterator tokit, const std::list<std::string> &tok){
    if(isFunction(opit)){
      auto next_tok = tokit; ++next_tok;
      if(next_tok == tok.end()){ std::cout << "Error: Got a function " << Op(opit) << " but nothing follows\n"; exit(-1); }
      else if( (*next_tok) != "("){  std::cout << "Error: Got a function " << Op(opit) << " not followed by an open parentheses\n"; exit(-1); }
    }
    stack.push_back(opit);
  }

  void actionOperator(expressionAST &AST, std::ostringstream &output, map_const_iterator opit, std::list<std::string>::const_iterator tokit, const std::list<std::string> &tok, std::string token, const std::string &prev_token){
    if(token == "-") handleUnaryMinus(opit,token,prev_token,tok,tokit);
	  	  
    if(stack.size() == 0 || Op(stack.back()) == "("){
      //If the incoming symbol is an operator and the stack is empty or contains a left parenthesis on top, push the incoming operator onto the stack.
      stack.push_back(opit);
    }else if(Prec(opit) > Prec(stack.back()) ||  ( Prec(opit) == Prec(stack.back()) && Assoc(opit) == Right)  ){
      //If the incoming symbol is an operator and has either higher precedence than the operator on the top of the stack, or has the same precedence as the operator on the top of the stack
      //and is right associative -- push it on the stack.
      stack.push_back(opit);
    }else{
      //If the incoming symbol is an operator and has either lower precedence than the operator on the top of the stack, or has the same precedence as the operator on the top of the stack
      //and is left associative -- continue to pop the stack until this is not true. Then, push the incoming operator.
      bool pop_cond = Prec(opit) < Prec(stack.back()) || ( Prec(opit) == Prec(stack.back()) && Assoc(opit) == Left);
      assert(pop_cond);

      while(pop_cond){
	output << Op(stack.back()) << " ";
	AST.stackOperator(Op(stack.back()));
	stack.pop_back();

	pop_cond = stack.size() > 0 && ( Prec(opit) < Prec(stack.back()) || ( Prec(opit) == Prec(stack.back()) && Assoc(opit) == Left) );
      }
      stack.push_back(opit);
    }
  }

  void actionRightBracket(expressionAST &AST, std::ostringstream &output, map_const_iterator opit, std::list<std::string>::const_iterator tokit, const std::list<std::string> &tok){
    //If the token is a right bracket (i.e. ")"), then:
    //   while the operator at the top of the operator stack is not a left bracket:
    // 	pop operators from the operator stack onto the output queue.
    //   pop the left bracket from the stack.
    //If the stack runs out without finding a left bracket, then there are
    //mismatched parentheses.
    while(stack.size() > 0){
      if(Op(stack.back()) == "(") break;
      output << Op(stack.back()) << " ";
      AST.stackOperator(Op(stack.back()));
      stack.pop_back();
    }
    if(stack.size() == 0){
      std::cout << "Error: mismatched parentheses (internal)\n";
      exit(-1);
    }
    assert(Op(stack.back()) == "(");
    stack.pop_back();

    //Check for a function
    if(stack.size() > 0 && isFunction(stack.back()) ){
      output << Op(stack.back()) << " ";
      AST.stackOperator(Op(stack.back()));
	
      stack.pop_back();
    }
  }

  void actionOperand(expressionAST &AST, std::ostringstream &output, const std::string &token){
    //If the token is a number or symbol, then push it to the output queue.
    output << token << " ";
    AST.stackOperand(token);  
  }
  void finalize(expressionAST &AST, std::ostringstream &output){
    // if there are no more tokens to read:
    // 	while there are still operator tokens on the stack:
    // 		/* if the operator token on the top of the stack is a bracket, then
    // 		there are mismatched parentheses. */
    // 		pop the operator onto the output queue.
    while(stack.size() > 0){
      if(Op(stack.back()) == "("){
	std::cout << "Error: mismatched parentheses (post)\n";
	exit(-1);
      }
      output << Op(stack.back()) << " ";
      AST.stackOperator(Op(stack.back()));
      stack.pop_back();
    }
  }
    
public:
  expressionAST parse(const std::string &s){
    std::list<std::string> tok = mathExpressionTokenize(s);
    std::ostringstream output; //for debugging
    
    expressionAST AST;
    
    map_const_iterator opit;
    std::string prev_token;
    for (std::list<std::string>::const_iterator tokit = tok.begin(); tokit != tok.end(); ++tokit){
      std::string token = *tokit;
      if((opit = getMap().find(token)) != getMap().end()){ //Is an operator / function / left bracket
	// If the token is a left bracket (i.e. "("), then push it onto the operator stack.
	// Same is true of functions
	if(Prec(opit) == NotAnOperator){
	  actionNotAnOperator(opit,tokit,tok);
	}else{
	  actionOperator(AST,output,opit,tokit,tok,token,prev_token);
	}
      }else if(token == ")"){
	actionRightBracket(AST,output,opit,tokit,tok);
      }else{
	actionOperand(AST,output,token);
      }
      prev_token = token;
    }
    finalize(AST,output);
    
    std::cout << "Expression: " << s << std::endl;
    std::cout << "Reordered: " << output.str() << std::endl;
    
    return AST;
  }
};


inline expressionAST mathExpressionParse(const std::string &s){
  shuntingYardParser parser;
  return parser.parse(s);
};

#endif
