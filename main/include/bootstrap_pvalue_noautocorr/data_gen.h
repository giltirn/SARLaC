#pragma once

class randomDataBase{
public:
  //Assumed to be thread safe
  virtual rawDataDistributionD generate(const int t, const int nsample) const = 0;
  virtual ~randomDataBase(){}
};
class randomDataGaussian: public randomDataBase{
  std::vector<double> mu;
  std::vector<double> sigma;
public:
  randomDataGaussian(int Lt, const std::vector<double> &mu, const std::vector<double> &sigma): mu(mu), sigma(sigma){
    if(mu.size() != Lt || sigma.size() != Lt) error_exit(std::cout << "mu, sigma size must equal Lt" << std::endl);
  }
  randomDataGaussian(int Lt, double mu_all, double sigma_all): mu(Lt, mu_all), sigma(Lt, sigma_all){}
  
  rawDataDistributionD generate(const int t, const int nsample) const override{
    rawDataDistributionD out(nsample);
    gaussianRandom(out, mu[t], sigma[t], threadRNG());
    return out;
  }
};
class randomDataLogNormal: public randomDataBase{
  std::vector<double> mu;
  std::vector<double> sigma;
public:
  randomDataLogNormal(int Lt, const std::vector<double> &mu, const std::vector<double> &sigma): mu(mu), sigma(sigma){
    if(mu.size() != Lt || sigma.size() != Lt) error_exit(std::cout << "mu, sigma size must equal Lt" << std::endl);
  }
  randomDataLogNormal(int Lt, double mu_all, double sigma_all): mu(Lt, mu_all), sigma(Lt, sigma_all){}
  
  rawDataDistributionD generate(const int t, const int nsample) const override{
    rawDataDistributionD mult(nsample);
    gaussianRandom(mult, 0, 1, threadRNG());
    mult = exp(mu[t] + sigma[t]*mult);
    return mult;
  }
};


//Parameters for random data generators with just mu, sigma shared for all timeslices
#define RDATA_UNIFORM_NRMLIKE (double, mu)(double, sigma)
struct RdataUniformNrmLikeArgs{
  GENERATE_MEMBERS(RDATA_UNIFORM_NRMLIKE); 
  RdataUniformNrmLikeArgs(): mu(0.), sigma(1.){  }
};
GENERATE_PARSER( RdataUniformNrmLikeArgs, RDATA_UNIFORM_NRMLIKE);

#define RDATA_TIMEDEP_NRMLIKE (std::vector<double>, mu)(std::vector<double>, sigma)
struct RdataTimeDepNrmLikeArgs{
  GENERATE_MEMBERS(RDATA_TIMEDEP_NRMLIKE); 
  RdataTimeDepNrmLikeArgs(): mu(10,0.), sigma(10,1.){  }
};
GENERATE_PARSER( RdataTimeDepNrmLikeArgs, RDATA_TIMEDEP_NRMLIKE);


std::unique_ptr<randomDataBase> dataGenStrategyFactory(DataGenStrategy strat, const std::string &params_file, const int Lt){
  if(strat == DataGenStrategy::NormalUniform){
    RdataUniformNrmLikeArgs args; parseOrTemplate(args, params_file, "datagen_template.args");
    return std::unique_ptr<randomDataBase>(new randomDataGaussian(Lt, args.mu, args.sigma));
  }else if(strat == DataGenStrategy::NormalTimeDep){
    RdataTimeDepNrmLikeArgs args; parseOrTemplate(args, params_file, "datagen_template.args");
    return std::unique_ptr<randomDataBase>(new randomDataGaussian(Lt, args.mu, args.sigma));
  }else if(strat == DataGenStrategy::LogNormalUniform){
    RdataUniformNrmLikeArgs args; parseOrTemplate(args, params_file, "datagen_template.args");
    return std::unique_ptr<randomDataBase>(new randomDataLogNormal(Lt, args.mu, args.sigma));
  }else{
    error_exit(std::cout << "Invalid data generation strategy" << std::endl);
  }
}

