#include<parser/expression_parse.h>

using namespace CPSfit;

//This is just a basic calculator with assignable symbols!
int main(const int argc, const char* argv[]){
  assert(argc >= 2);
  assert(argc % 2 == 0);
  std::string s(argv[1]);

  std::vector<std::pair<std::string,double> > symbols;
  for(int i=2;i<argc;i+=2){
    std::string sym = argv[i];
    double val;
    std::stringstream ss; ss << argv[i+1]; ss >> val;
    symbols.push_back(std::pair<std::string,double>(sym,val));
  }
  
  expressionAST AST =  mathExpressionParse(s);

  for(int i=0;i<symbols.size();i++){
    AST.assignSymbol(symbols[i].first,symbols[i].second);
  }

  std::cout << "Value: " << AST.value() << std::endl;
  
}
