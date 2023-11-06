#pragma once

class randomDataBase{
public:
  //Assumed to be thread safe
  virtual correlationFunction<double, rawDataDistributionD> generate(const int Lt, const int nsample) const = 0;
  virtual ~randomDataBase(){}
};

std::unique_ptr<randomDataBase> dataGenStrategyFactory(DataGenStrategy strat, const std::string &params_file, const int Lt);

class randomDataGaussian: public randomDataBase{
  std::vector<double> mu;
  std::vector<double> sigma;
public:
  randomDataGaussian(int Lt, const std::vector<double> &mu, const std::vector<double> &sigma): mu(mu), sigma(sigma){
    if(mu.size() != Lt || sigma.size() != Lt) error_exit(std::cout << "mu, sigma size must equal Lt" << std::endl);
  }
  randomDataGaussian(int Lt, double mu_all, double sigma_all): mu(Lt, mu_all), sigma(Lt, sigma_all){}
  
  correlationFunction<double, rawDataDistributionD> generate(const int Lt, const int nsample) const override{
    correlationFunction<double, rawDataDistributionD> out(Lt, nsample);
    for(int t=0;t<Lt;t++){
      out.coord(t) = t;
      gaussianRandom(out.value(t), mu[t], sigma[t], threadRNG());
    }
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
  
  correlationFunction<double, rawDataDistributionD> generate(const int Lt, const int nsample) const override{
    correlationFunction<double, rawDataDistributionD> out(Lt, nsample);
    for(int t=0;t<Lt;t++){
      out.coord(t) = t;
      rawDataDistributionD &mult = out.value(t);
      gaussianRandom(mult, 0, 1, threadRNG());
      mult = exp(mu[t] + sigma[t]*mult);
    }
    return out;
  }
};
class randomDataGaussianPlusShift: public randomDataBase{
  std::vector<double> mu;
  std::vector<double> sigma;
  double shift_mu;
  double shift_sigma;
  double shift_alpha; //amount of the common shift added to each timeslice
public:
  randomDataGaussianPlusShift(int Lt, const std::vector<double> &mu, const std::vector<double> &sigma,
			      double shift_mu, double shift_sigma, double shift_alpha): mu(mu), sigma(sigma), shift_mu(shift_mu), shift_sigma(shift_sigma), shift_alpha(shift_alpha){
    if(mu.size() != Lt || sigma.size() != Lt) error_exit(std::cout << "mu, sigma size must equal Lt" << std::endl);
  }
  randomDataGaussianPlusShift(int Lt, double mu_all, double sigma_all, double shift_mu, double shift_sigma, double shift_alpha): mu(Lt, mu_all), sigma(Lt, sigma_all), shift_mu(shift_mu), shift_sigma(shift_sigma), shift_alpha(shift_alpha){}
  
  correlationFunction<double, rawDataDistributionD> generate(const int Lt, const int nsample) const override{
    rawDataDistributionD shift(nsample);
    gaussianRandom(shift, shift_mu, shift_sigma, threadRNG());

    correlationFunction<double, rawDataDistributionD> out(Lt, nsample);
    for(int t=0;t<Lt;t++){
      out.coord(t) = t;
      rawDataDistributionD &ov = out.value(t);
      gaussianRandom(ov, mu[t], sigma[t], threadRNG());
      ov = ov + shift_alpha * shift;
    }
    return out;
  }
};



class randomDataBinned: public randomDataBase{
  int bin_size;
  std::unique_ptr<randomDataBase> base_gen;
public:
  //nsample in input args should be size *after* binning
  randomDataBinned(int Lt, int bin_size, DataGenStrategy strat, const std::string &base_params_file): bin_size(bin_size), 
												      base_gen(dataGenStrategyFactory(strat,base_params_file,Lt)){}

  correlationFunction<double, rawDataDistributionD> generate(const int Lt, const int nsample) const override{
    correlationFunction<double, rawDataDistributionD> b = base_gen->generate(Lt,nsample*bin_size);
    correlationFunction<double, rawDataDistributionD> out(Lt);
    for(int t=0;t<Lt;t++){
      out.coord(t) = t;
      out.value(t) = b.value(t).bin(bin_size);
    }
    return out;
  }
};

//Parameters for random data generators with just mu, sigma shared for all timeslices
#define MEMBERS (double, mu)(double, sigma)
struct RdataUniformNrmLikeArgs{
  GENERATE_MEMBERS(MEMBERS); 
  RdataUniformNrmLikeArgs(): mu(0.), sigma(1.){  }
};
GENERATE_PARSER( RdataUniformNrmLikeArgs, MEMBERS);
#undef MEMBERS

#define MEMBERS (std::vector<double>, mu)(std::vector<double>, sigma)
struct RdataTimeDepNrmLikeArgs{
  GENERATE_MEMBERS(MEMBERS); 
  RdataTimeDepNrmLikeArgs(): mu(10,0.), sigma(10,1.){  }
};
GENERATE_PARSER( RdataTimeDepNrmLikeArgs, MEMBERS);
#undef MEMBERS

#define MEMBERS (double, mu)(double, sigma)(double, shift_mu)(double, shift_sigma)(double, shift_alpha)
struct RdataUniformNrmLikePlusShiftArgs{
  GENERATE_MEMBERS(MEMBERS); 
  RdataUniformNrmLikePlusShiftArgs(): mu(0.), sigma(1.), shift_mu(0.), shift_sigma(1.), shift_alpha(0.1){  }
};
GENERATE_PARSER( RdataUniformNrmLikePlusShiftArgs, MEMBERS);
#undef MEMBERS

#define MEMBERS (int, bin_size)(DataGenStrategy, base_strat)(std::string, base_params_file)
struct BinnedDataArgs{
  GENERATE_MEMBERS(MEMBERS); 
  BinnedDataArgs(): bin_size(1), base_strat(DataGenStrategy::NormalUniform), base_params_file("datagen_base.args"){  }
};
GENERATE_PARSER( BinnedDataArgs, MEMBERS);
#undef MEMBERS

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
  }else if(strat == DataGenStrategy::NormalUniformPlusShift){
    RdataUniformNrmLikePlusShiftArgs args; parseOrTemplate(args, params_file, "datagen_template.args");
    return std::unique_ptr<randomDataBase>(new randomDataGaussianPlusShift(Lt, args.mu, args.sigma, args.shift_mu, args.shift_sigma, args.shift_alpha));
  }else if(strat == DataGenStrategy::Binned){
    BinnedDataArgs args; parseOrTemplate(args, params_file, "datagen_template.args");
    return std::unique_ptr<randomDataBase>(new randomDataBinned(Lt, args.bin_size, args.base_strat, args.base_params_file));
  }else{
    error_exit(std::cout << "Invalid data generation strategy" << std::endl);
  }
}

