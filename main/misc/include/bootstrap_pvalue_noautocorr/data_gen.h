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
class randomDataGaussianMixLeft: public randomDataBase{
  std::vector<double> mu;
  std::vector<double> sigma;
  std::vector<double> alpha; //fraction of data point at t-1 to include

public:
  randomDataGaussianMixLeft(int Lt, const std::vector<double> &mu, const std::vector<double> &sigma, const std::vector<double> &alpha): mu(mu), sigma(sigma), alpha(alpha){
    if(mu.size() != Lt || sigma.size() != Lt || alpha.size() != Lt) error_exit(std::cout << "mu, sigma, alpha size must equal Lt" << std::endl);
  }
  randomDataGaussianMixLeft(int Lt, double mu_all, double sigma_all, const std::vector<double> &alpha): mu(Lt, mu_all), sigma(Lt, sigma_all), alpha(alpha){}
  
  correlationFunction<double, rawDataDistributionD> generate(const int Lt, const int nsample) const override{
    correlationFunction<double, rawDataDistributionD> out(Lt, nsample);
    for(int t=0;t<Lt;t++){
      out.coord(t) = t;
      rawDataDistributionD &ov = out.value(t);
      gaussianRandom(ov, mu[t], sigma[t], threadRNG());

      if(t>0) ov = alpha[t]*out.value(t-1) + (1-alpha[t])*ov;
    }
    return out;
  }
};


class randomDataBinned: public randomDataBase{
  int bin_size;
  int nsample_unbinned;
  std::unique_ptr<randomDataBase> base_gen;
public:
  //nsample in input args should be size *after* binning
  randomDataBinned(int Lt, int bin_size, int nsample_unbinned, DataGenStrategy strat, const std::string &base_params_file): bin_size(bin_size), nsample_unbinned(nsample_unbinned),
												      base_gen(dataGenStrategyFactory(strat,base_params_file,Lt)){}

  correlationFunction<double, rawDataDistributionD> generate(const int Lt, const int nsample) const override{
    if(nsample_unbinned / bin_size != nsample) error_exit(std::cout << "Expect nsample_unbinned ("<<nsample_unbinned<<") / bin_size (" << bin_size <<") == nsample (" << nsample << 
							  ") got" << nsample_unbinned / bin_size << std::endl); //allow for truncation

    correlationFunction<double, rawDataDistributionD> b = base_gen->generate(Lt,nsample_unbinned); //this way we can be sure we always get the same data even if the bin size is not an exact divisor of nsample_unbinned
    correlationFunction<double, rawDataDistributionD> out(Lt);
    for(int t=0;t<Lt;t++){
      out.coord(t) = t;
      out.value(t) = b.value(t).bin(bin_size,true);
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

#define MEMBERS (double, mu)(double, sigma)(std::vector<double>, alpha)
struct RdataUniformNrmLikeMixLeftArgs{
  GENERATE_MEMBERS(MEMBERS); 
  RdataUniformNrmLikeMixLeftArgs(): mu(0.), sigma(1.), alpha(10,0.1){  }
};
GENERATE_PARSER( RdataUniformNrmLikeMixLeftArgs, MEMBERS);
#undef MEMBERS

#define MEMBERS (std::vector<double>, mu)(std::vector<double>, sigma)(std::vector<double>, alpha)
struct RdataTimeDepNrmLikeMixLeftArgs{
  GENERATE_MEMBERS(MEMBERS); 
  RdataTimeDepNrmLikeMixLeftArgs(): mu(10,0.), sigma(10,1.), alpha(10,0.1){  }
};
GENERATE_PARSER( RdataTimeDepNrmLikeMixLeftArgs, MEMBERS);
#undef MEMBERS


#define MEMBERS (int, bin_size)(int, nsample_unbinned)(DataGenStrategy, base_strat)(std::string, base_params_file)
struct BinnedDataArgs{
  GENERATE_MEMBERS(MEMBERS); 
  BinnedDataArgs(): bin_size(1), nsample_unbinned(100),base_strat(DataGenStrategy::NormalUniform), base_params_file("datagen_base.args"){  }
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
  }else if(strat == DataGenStrategy::NormalUniformMixLeft){
    RdataUniformNrmLikeMixLeftArgs args; parseOrTemplate(args, params_file, "datagen_template.args");
    return std::unique_ptr<randomDataBase>(new randomDataGaussianMixLeft(Lt, args.mu, args.sigma, args.alpha));
  }else if(strat == DataGenStrategy::NormalTimeDepMixLeft){
    RdataTimeDepNrmLikeMixLeftArgs args; parseOrTemplate(args, params_file, "datagen_template.args");
    return std::unique_ptr<randomDataBase>(new randomDataGaussianMixLeft(Lt, args.mu, args.sigma, args.alpha));
  }else if(strat == DataGenStrategy::Binned){
    BinnedDataArgs args; parseOrTemplate(args, params_file, "datagen_template.args");
    return std::unique_ptr<randomDataBase>(new randomDataBinned(Lt, args.bin_size, args.nsample_unbinned, args.base_strat, args.base_params_file));
  }else{
    error_exit(std::cout << "Invalid data generation strategy" << std::endl);
  }
}

