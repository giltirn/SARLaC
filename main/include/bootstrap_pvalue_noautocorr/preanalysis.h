#pragma once

struct preAnalysisBase{
  virtual void run(const Args &args, const covMatStrategyBase &covgen, const randomDataBase &datagen) const = 0;
  virtual ~preAnalysisBase(){}
};
struct preAnalysisNone: public preAnalysisBase{
  void run(const Args &args, const covMatStrategyBase &covgen, const randomDataBase &datagen) const override{};
};
struct preAnalysisCovMatEvals: public preAnalysisBase{
  void run(const Args &args, const covMatStrategyBase &covgen, const randomDataBase &datagen) const override{}
};



std::unique_ptr<preAnalysisBase> preAnalysisFactory(preAnalysisType type){
  if(type == preAnalysisType::None){
    return std::unique_ptr<preAnalysisBase>(new preAnalysisNone);
  }else if(type == preAnalysisType::CovMatEvals){
    return std::unique_ptr<preAnalysisBase>(new preAnalysisCovMatEvals);
  }else{
    error_exit(std::cout << "Invalid pre-analysis type" << std::endl);
  }
}
