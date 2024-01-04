#pragma once

class randomDataBase{
  std::map<int, std::vector<double> > timeslice_means_cache; //cache computation of timeslice means to avoid needing to recompute
public:
  //Assumed to be thread safe
  virtual correlationFunction<double, rawDataDistributionD> generate(const int Lt, const int nsample) const = 0;
  virtual std::vector<double> populationTimesliceMeans() const = 0;
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

  std::vector<double> populationTimesliceMeans() const override{
    return mu;
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

  std::vector<double> populationTimesliceMeans() const override{
    std::vector<double> err(mu.size());
    for(int t=0;t<mu.size();t++)
      err[t] = exp( mu[t] + sigma[t]*sigma[2]/2. );
    return err;
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
  std::vector<double> populationTimesliceMeans() const override{
    std::vector<double> err(mu.size());
    for(int t=0;t<mu.size();t++)
      err[t] = mu[t] + shift_alpha * shift_mu;
    return err;
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
  std::vector<double> populationTimesliceMeans() const override{
    std::vector<double> err(mu.size());
    for(int t=0;t<mu.size();t++){
      err[t] = mu[t];
      if(t>0) err[t] = alpha[t]*err[t-1] + (1-alpha[t])*err[t]; //TODO: Check this
    }
    return err;
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
  std::vector<double> populationTimesliceMeans() const override{
    return base_gen->populationTimesliceMeans();
  }
};

template<typename Ffunc, typename UpdateFunc>
class MetropolisDataGen{
  Ffunc ffunc;
  UpdateFunc update;

  int nsample_therm; //thermalization, paid once
  int nsample_decorr; //decorrelation, paid every timeslice (note, we assume that we wish each timeslice's data to be independent)

  mutable std::vector<double> t_state; //thread state
public:
  MetropolisDataGen(const Ffunc &ffunc, const UpdateFunc &update, int nsample_therm, int nsample_decorr, double init): ffunc(ffunc), update(update), nsample_therm(nsample_therm), nsample_decorr(nsample_decorr), t_state(omp_get_max_threads()){
    int nthr = omp_get_max_threads();
    rawDataDistribution<double> setup = MetropolisHastings(nthr*nsample_decorr, nsample_therm, init, ffunc, update, RNG); //main thread
    std::cout << "MetropolisDataGen initialized for " << nthr << " thread streams with starts:";
    for(int i=0;i<nthr;i++){
      t_state[i] = setup.sample(i*nsample_decorr);
      std::cout << " " << t_state[i] << std::endl;
    }
    std::cout << std::endl;
  }

  //Compute the acceptance by running a new Metropolis chain of a certain length with a new RNG (so as not to change state)
  double acceptance() const{
    assert(!omp_in_parallel());
    RNGstore rng(RNG);
    double a;
    rawDataDistributionD tmp = MetropolisHastings(0,200000, t_state[0], ffunc, update, rng, &a);
    return a;
  }

  rawDataDistributionD generate(const int nsample, const int traj_inc = 1) const{
    int me = omp_get_thread_num();
    if(me > t_state.size()) error_exit(std::cout << "Thread index is larger than the initial max_threads" << std::endl);
    rawDataDistributionD tmp = MetropolisHastings(nsample*traj_inc, nsample_decorr, t_state[me], ffunc, update, threadRNG());
    t_state[me] = tmp.sample(nsample*traj_inc-1); //update thread state
    if(traj_inc > 1){
      rawDataDistributionD out(nsample);
      for(int i=0;i<nsample;i++) out.sample(i) = tmp.sample(traj_inc*i);
      return out;
    }else{
      return tmp;
    }
  }
};


//tParamsType is a struct containing the model and update func params for a given timeslice. Must have operator<,  FFunc gen_func(), UpdateFunc gen_update()
template<typename tParamsType>
class randomDataMetropolisBase: public randomDataBase{  
  typedef typename tParamsType::Ffunc Ffunc;
  typedef typename tParamsType::UpdateFunc UpdateFunc;
  typedef MetropolisDataGen<Ffunc,UpdateFunc> genType;
  std::map<tParamsType, std::unique_ptr<genType> > gen_m; 
  std::vector<genType*> gen_t;
  int traj_inc;

public:
  randomDataMetropolisBase(const std::vector<tParamsType> &tslice_params, int nsample_therm, int nsample_decorr, int traj_inc, double init): traj_inc(traj_inc){
    int Lt=tslice_params.size();
    gen_t.resize(Lt,nullptr);
    for(int t=0;t<Lt;t++){
      const tParamsType &tp = tslice_params[t];
      auto it = gen_m.find(tp);
      if(it == gen_m.end()){
	std::cout << "For timeslice t=" << t << " starting NEW Metropolis chain with params " << tp << std::endl;
	gen_m[tp].reset(new genType(tp.gen_func(), tp.gen_update(), nsample_therm, nsample_decorr, init));
	std::cout << "Got acceptance " << gen_m[tp]->acceptance() << std::endl;
      }else{
	std::cout << "For timeslice t=" << t << " reusing existing Metropolis chain with params " << tp << std::endl;
      }
      gen_t[t] = gen_m[tp].get();
    }
  }

  //nsample is the number of samples in the output distribution, independent of the sample frequency traj_inc
  correlationFunction<double, rawDataDistributionD> generate(const int Lt, const int nsample) const override{
    if(Lt > gen_t.size()) error_exit(std::cout << "Lt is larger than #timeslices for which the generator has been initialized" << std::endl);

    int me = omp_get_thread_num();
    correlationFunction<double, rawDataDistributionD> out(Lt, nsample);
    for(int t=0;t<Lt;t++){
      out.coord(t) = t;
      out.value(t) = gen_t[t]->generate(nsample,traj_inc);
    }
    return out;
  }  
};

struct GaussianMetropolisTimesliceParams{
  typedef ProbGaussian Ffunc;
  typedef UpdateFuncGaussian UpdateFunc;
  std::array<double,3> p;
  GaussianMetropolisTimesliceParams(){}
  GaussianMetropolisTimesliceParams(double mu, double sigma, double update_width): p({mu,sigma,update_width}){}

  inline Ffunc gen_func() const{ return Ffunc(p[0],p[1]); }
  inline UpdateFunc gen_update() const{ return UpdateFunc(p[2]); }
  bool operator<(const GaussianMetropolisTimesliceParams &r) const{ return p<r.p; }
};      
std::ostream & operator<<(std::ostream &os, const GaussianMetropolisTimesliceParams &p){
  os << "(mu=" << p.p[0] <<",sigma=" << p.p[1] << ",width=" << p.p[2] <<")"; 
  return os;
}

class randomDataMetropolisGaussian: public randomDataMetropolisBase<GaussianMetropolisTimesliceParams>{
  std::vector<double> _mu;

  static inline std::vector<GaussianMetropolisTimesliceParams> genParams(int Lt, const std::vector<double> &mu, const std::vector<double> &sigma, const std::vector<double> &update_width){
    if(mu.size() != Lt || sigma.size() != Lt || update_width.size() != Lt) error_exit(std::cout << "mu, sigma, update_width size must equal Lt" << std::endl);
    std::vector<GaussianMetropolisTimesliceParams> out(Lt); 
    for(int t=0;t<Lt;t++) out[t] = GaussianMetropolisTimesliceParams(mu[t],sigma[t],update_width[t]);
    return out;
  }

public:
  randomDataMetropolisGaussian(int Lt, const std::vector<double> &mu, const std::vector<double> &sigma, const std::vector<double> &update_width, 
			       int nsample_therm, int nsample_decorr, int traj_inc, double init): randomDataMetropolisBase<GaussianMetropolisTimesliceParams>(genParams(Lt,mu,sigma,update_width),nsample_therm,nsample_decorr,traj_inc,init),  _mu(mu){}

  randomDataMetropolisGaussian(int Lt, const double mu_all, const double sigma_all, const double update_width_all, 
			       int nsample_therm, int nsample_decorr, int traj_inc, double init): randomDataMetropolisBase<GaussianMetropolisTimesliceParams>(genParams(Lt,std::vector<double>(Lt,mu_all),std::vector<double>(Lt,sigma_all),std::vector<double>(Lt,update_width_all)),nsample_therm,nsample_decorr,traj_inc,init),  _mu(Lt,mu_all){}

  std::vector<double> populationTimesliceMeans() const override{ return _mu; }
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

#define MEMBERS (double, mu)(double, sigma)(double, update_width)(int, nsample_therm)(int, nsample_decorr)(int, traj_inc)(double, init)
struct RdataUniformNrmLikeMetropolisArgs{
  GENERATE_MEMBERS(MEMBERS); 
  RdataUniformNrmLikeMetropolisArgs(): mu(0.), sigma(1.), update_width(1.0), nsample_therm(50000), nsample_decorr(1000), traj_inc(4), init(0.0){  }
};
GENERATE_PARSER( RdataUniformNrmLikeMetropolisArgs, MEMBERS);
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
  }else if(strat == DataGenStrategy::NormalUniformMetropolis){
    RdataUniformNrmLikeMetropolisArgs args; parseOrTemplate(args, params_file, "datagen_template.args");
    return std::unique_ptr<randomDataBase>(new randomDataMetropolisGaussian(Lt, args.mu, args.sigma, args.update_width, args.nsample_therm, args.nsample_decorr, args.traj_inc, args.init));
  }else{
    error_exit(std::cout << "Invalid data generation strategy" << std::endl);
  }
}

