#include<physics/alpha_s.h>

CPSFIT_START_NAMESPACE

struct as_params{
  double mu;
  double match_to;
  int Nf;
  int Nc;
};

struct as_accessor{
  static double alpha_s(const as_params& p, const double Lambda){
    return ComputeAlphaS::alpha_s(p.mu, Lambda, p.Nf, p.Nc);
  }
};

inline double as_f (double Lambda, void * p) {
  const as_params& params = *((as_params*)p);
  double as = as_accessor::alpha_s(params,Lambda);
  double diff = as - params.match_to;
  double diff2 = diff*diff; //to minimize

  //printf("alpha_s(Lambda=%f)=%f, target %f,  diff2 = %f, target 0\n",Lambda,as,params.match_to,diff2);
  return diff2;
}

//Compute Lambda_QCD by matching alpha_s in Nf-flavor theory to an input value
//match_at is an energy scale, match_to is a value of alpha_s
double ComputeAlphaS::computeLambda(const int Nf, const double match_at, const double match_to, const int Nc, bool vrb){
  if(vrb) printf("Starting computeLambda with Nf=%d, match_at=%f match_to=%f, Nc=%d\n",Nf,match_at,match_to,Nc);
  as_params p;
  p.mu = match_at;
  p.match_to = match_to;
  p.Nf = Nf;
  p.Nc = Nc;

  gsl_function F;
  F.function = &as_f;
  F.params = &p;

  //Note, due to log(log(mu^2/Lambda^2)),   Lambda must always be smaller than mu or else we get a log of a negative number
  double start = 1e-6;
  double guess = 1e-5;
  double end = 0.99999*match_at;

  const gsl_min_fminimizer_type * T = gsl_min_fminimizer_brent;
  gsl_min_fminimizer * s = gsl_min_fminimizer_alloc (T);
  gsl_min_fminimizer_set (s, &F, guess, start, end);

  int status;
  int iter = 0, max_iter = 100;
  double m,a,b;
  do{
    iter++;
    status = gsl_min_fminimizer_iterate (s);

    m = gsl_min_fminimizer_x_minimum (s);
    a = gsl_min_fminimizer_x_lower (s);
    b = gsl_min_fminimizer_x_upper (s);

    status = gsl_min_test_interval (a, b, 1e-7, 0.0);

    if(vrb && status == GSL_SUCCESS){
      printf ("Converged:\n");

      printf("%5d [%.7f, %.7f] "
	     "%.7f %.7f\n",
	     iter, a, b,
	     m, b - a);
    }
      
  }
  while (status == GSL_CONTINUE && iter < max_iter);  
  gsl_min_fminimizer_free(s);

  if(status != GSL_SUCCESS){
    printf("Converge failed\n");
    exit(-1);
  }

  return m;
}


//Set to values from Qi's thesis
void PerturbativeInputs::setQisValues(){
  mcmc = 1.27;
  mbmb = 4.19;
  mWmW = 80.399;
  mtmt = 170;
  a_e = 1./128;    
  thetaW = asin(sqrt(0.23116));
  Nc = 3;

  double Lambda4 = 0.3298655;
  double as4mb = ComputeAlphaS::alpha_s(mbmb, Lambda4, 4, Nc);

  a_s.value = as4mb;
  a_s.scale = mbmb;
  a_s.Nf = 4;
}
//Set to values from ep' PRL
void PerturbativeInputs::setepPRLValues(){
  mcmc = 1.275;
  mbmb = 4.18;
  mWmW = 80.385;
  mtmt = 160;
  a_e = 1./127.94;    
  thetaW = asin(sqrt(0.23126));
  Nc = 3;

  double Lambda4 = 0.331416;
  double as4mb = ComputeAlphaS::alpha_s(mbmb, Lambda4, 4, Nc);

  a_s.value = as4mb;
  a_s.scale = mbmb;
  a_s.Nf = 4;
}
//Set to values from Ziyuan Bai's thesis pg 67
void PerturbativeInputs::setepZiyuanValues(){
  mcmc = 1.275; //unknown, leave as PRL
  mbmb = 4.19;
  mWmW = 80.4;
  mtmt = 172.2;
  a_e = 1./128; //unknown    
  thetaW = asin(sqrt(0.23126)); //unknown
  Nc = 3;

  double mz = 91.1876;
  double Lambda5 = 0.1184;
  double as5mz = ComputeAlphaS::alpha_s(mz, Lambda5, 5, Nc);

  a_s.value = as5mz;
  a_s.scale = mz;
  a_s.Nf = 5;
}
//Set to values from epsilon' 2019 paper
void PerturbativeInputs::setep2019values(){
  Nc = 3;
  mcmc = 1.27;
  mbmb = 4.18;
  mWmW = 80.379;
  mtmt = 160;
  thetaW = 0.50162777299552883387;
  a_e = 0.00781524754796608183;
  a_s.value = 0.1181;
  a_s.scale = 91.1876;
  a_s.Nf = 5;
}



void PerturbativeVariables::initialize(const PerturbativeInputs &_inputs){    
  inputs = _inputs;

  if(vrb) printf("Initializing perturbative variables\n");
     
  if(inputs.a_s.Nf == 3){
    Lambda[3 -3] = ComputeAlphaS::computeLambda(3, inputs.a_s.scale, inputs.a_s.value, inputs.Nc, vrb);
       
    double as_mc = ComputeAlphaS::alpha_s(inputs.mcmc, Lambda[3 -3], 3, inputs.Nc);

    Lambda[4 -3] = ComputeAlphaS::computeLambda(4, inputs.mcmc, as_mc, inputs.Nc, vrb);
       
    double as_mb = ComputeAlphaS::alpha_s(inputs.mbmb, Lambda[4 -3], 4, inputs.Nc);

    Lambda[5 -3] = ComputeAlphaS::computeLambda(5, inputs.mbmb, as_mb, inputs.Nc, vrb);
  }else if(inputs.a_s.Nf == 4){
    Lambda[4 -3] = ComputeAlphaS::computeLambda(4, inputs.a_s.scale, inputs.a_s.value, inputs.Nc, vrb);
       
    double as_mc = ComputeAlphaS::alpha_s(inputs.mcmc, Lambda[4 -3], 4, inputs.Nc);

    Lambda[3 -3] = ComputeAlphaS::computeLambda(3, inputs.mcmc, as_mc, inputs.Nc, vrb);

    double as_mb = ComputeAlphaS::alpha_s(inputs.mbmb, Lambda[4 -3], 4, inputs.Nc);

    Lambda[5 -3] = ComputeAlphaS::computeLambda(5, inputs.mbmb, as_mb, inputs.Nc, vrb);
  }else if(inputs.a_s.Nf == 5){
    Lambda[5 -3] = ComputeAlphaS::computeLambda(5, inputs.a_s.scale, inputs.a_s.value, inputs.Nc, vrb);

    double as_mb = ComputeAlphaS::alpha_s(inputs.mbmb, Lambda[5 -3], 5, inputs.Nc);

    Lambda[4 -3] = ComputeAlphaS::computeLambda(4, inputs.mbmb, as_mb, inputs.Nc, vrb);

    double as_mc = ComputeAlphaS::alpha_s(inputs.mcmc, Lambda[4 -3], 4, inputs.Nc);
       
    Lambda[3 -3] = ComputeAlphaS::computeLambda(3, inputs.mcmc, as_mc, inputs.Nc, vrb);
  }else assert(0);


  if(vrb){
    printf("Computed Lambda_5 = %.7f\n", Lambda[5 -3]);
    printf("Computed Lambda_4 = %.7f\n", Lambda[4 -3]);
    printf("Computed Lambda_3 = %.7f\n", Lambda[3 -3]);      
  }
}

CPSFIT_END_NAMESPACE
