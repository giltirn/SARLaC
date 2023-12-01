#ifndef _EXPRESSION_PARSE_SHUNTING_YARD_H_
#define _EXPRESSION_PARSE_SHUNTING_YARD_H_

#include<config.h>
#include<utils/macros.h>

#include<parser/expression_parse/AST.h>
#include<parser/expression_parse/tokenize.h>

SARLAC_START_NAMESPACE

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
      operators["pow"] = std::pair<int,int>(NotAnOperator,NotAnOperator);
      initted = true;
    }
    return operators;
  }

  static inline const std::string & Op(std::map<std::string, std::pair<int,int> >::const_iterator mit){ return mit->first; }
  static inline int Prec(std::map<std::string, std::pair<int,int> >::const_iterator mit){ return mit->second.first; }
  static inline int Assoc(std::map<std::string, std::pair<int,int> >::const_iterator mit){ return mit->second.second; }
  static inline bool isFunction(std::map<std::string, std::pair<int,int> >::const_iterator mit){ return (Prec(mit) == -1 && Op(mit) != "("); }

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


SARLAC_END_NAMESPACE
#endif
