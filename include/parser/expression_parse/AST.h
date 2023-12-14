#ifndef _EXPRESSION_PARSE_AST_H_
#define _EXPRESSION_PARSE_AST_H_

#include<iostream>
#include<map>
#include<sstream>

#include<boost/spirit/home/x3.hpp>

#include<config.h>
#include<utils/macros.h>

SARLAC_START_NAMESPACE

//An abstract syntax tree (AST) representing a mathematical expression
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
  struct BinOpDivide: public BinOp{
    inline double value() const{ return this->a->value() / this->b->value(); }
  };
  struct BinOpPow: public BinOp{
    inline double value() const{ return ::pow(this->a->value(),this->b->value()); }
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
    }else if(op == "/"){
      binop = new BinOpDivide;
    }else if(op == "^" || op == "pow"){
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


SARLAC_END_NAMESPACE
#endif
