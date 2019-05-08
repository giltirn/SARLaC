#ifndef _CPSFIT_METROPOLIS_HASTINGS_H
#define _CPSFIT_METROPOLIS_HASTINGS_H

//An implementation of the Metropolis-Hastings Markov chain Monte Carlo algorithm

#include<config.h>
#include<utils/template_wizardry.h>
#include<random/random_number.h>
#include<distribution/raw_data_distribution.h>

CPSFIT_START_NAMESPACE

//Example update function: update x -> x' where x' is drawn from a Gaussian with x as central value and some chosen width
struct UpdateFuncGaussian{
  double width;

  UpdateFuncGaussian(const double width): width(width){}  

  double operator()(const double x) const{ return gaussianRandom<double>(x,width); }
};

//Example probability function - Gaussian
struct ProbGaussian{
  double mu;
  double sigma;
  double _2s2;

  ProbGaussian(const double mu, const double sigma): mu(mu), sigma(sigma), _2s2(2.*sigma*sigma){}

  inline double operator()(const double x) const{ double d = x-mu;  return exp(-d*d/_2s2); }
};

//Update func generates a new candidate x_i+1 based on current x_i. Must be reversible  g(x|y) = g(y|x)
//Ffunc is a function proportional to the target distribution
template<typename Ffunc, typename UpdateFunc>
rawDataDistribution<double> MetropolisHastings(const int Nkeep, const int Nwarmup, const double x0, const Ffunc &f, const UpdateFunc &g){
  double x = x0;

  rawDataDistribution<double> out(Nkeep);

  int naccept = 0;

  for(int i=0;i<Nwarmup+Nkeep;i++){
    if(i>=Nwarmup) out.sample(i-Nwarmup) = x;

    double xp = g(x);
    double alpha = f(xp)/f(x);
    double u = uniformRandom<double>(0.0,1.0);

    bool accept = u <= alpha;
    if(accept) naccept++;

    x = accept ? xp : x;
  }

  std::cout << "Acceptance " << double(naccept)/(Nwarmup + Nkeep) << std::endl;

  return out;
} 


CPSFIT_END_NAMESPACE

#endif
