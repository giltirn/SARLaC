#include<tensors.h>
#include<random.h>
#include<plot.h>
#include<minimizer.h>

using namespace SARLaC;

typedef NumericSquareMatrix<double> Matrix;
typedef NumericRectMatrix<double> RectMatrix;
typedef NumericVector<double> Vector;

Vector randMVN(const Vector &mu, const Matrix &Sigma){
  int N = mu.size(); assert(Sigma.size() == N);
  Vector v(N);
  for(int i=0;i<N;i++) v(i) = gaussianRandom<double>(0,1);
  return mu + Cholesky(Sigma) * v;
}
inline double cov(double x1, double x2, double sigma, double lambda){ return sigma*sigma*exp(-pow( (x1-x2)/lambda, 2 )); }

Vector vrange(double start, double inc, double lessthan){
  Vector out;
  for(double v=start; v<lessthan; v+=inc) out.push_back(v);
  return out;
}

struct VectorAccessor : public CurveDataAccessorBase<double>{
  const Vector &xx;
  const Vector &yy;
  
  VectorAccessor(const Vector &xx, const Vector &yy): xx(xx), yy(yy){}
  double x(const int i) const override{ return xx(i); }
  double y(const int i) const override{ return yy(i); }
  int size() const override{ return xx.size(); }
};
struct VectorAccessorErrBand{
  const Vector &xx;
  const Vector &yy;
  const Vector &ss;
  
  VectorAccessorErrBand(const Vector &xx, const Vector &yy, const Vector &ss): xx(xx), yy(yy), ss(ss){}
  double upper(const int i) const{ return yy[i] + ss[i]; };
  double lower(const int i) const{ return yy[i] - ss[i]; };
  double x(const int i) const{ return xx[i]; }

  int size() const{ return xx.size(); }
};


typedef MatPlotLibScriptGenerateBase::kwargsType kwargsType;

//Gaussian process "regression" for points xp using training data (xt,yt)
//Prior correlations sigma, lambda and a function mu for the prior means
template<typename MuFuncType>
void GPfit(Vector &yp, Matrix &Sp, const Vector &xp, const Vector &yt, const Vector &xt, const MuFuncType &mu, const double sigma, const double lambda, double epsilon = 1e-10){
  int ntrain = xt.size(); assert(yt.size() == xt.size());
  int npred = xp.size();
 
  Matrix Stt(ntrain, [&](const int i,const int j){ return cov(xt[i],xt[j],sigma,lambda); });
  RectMatrix Spt(npred,ntrain, [&](const int i,const int j){ return cov(xp[i],xt[j],sigma,lambda); });
  RectMatrix Spp(npred,npred, [&](const int i,const int j){ return cov(xp[i],xp[j],sigma,lambda); });

  Matrix Sttinv(ntrain);
  svd_inverse(Sttinv,Stt);

  Vector mu_xt(ntrain);
  for(int i=0;i<ntrain;i++) mu_xt(i) = mu(xt(i));

  yp = Spt*(Sttinv*(yt - mu_xt));
  for(int i=0;i<npred;i++) yp(i) += mu(xp(i));

  RectMatrix Sp_r = Spp - (Spt*Sttinv)*Spt.transpose();
  Sp.resize(npred);
  for(int i=0;i<npred;i++)
    for(int j=0;j<npred;j++)
      Sp(i,j) = Sp_r(i,j);
}

void normalizeVec(Vector &v){
  double s = 0;
  double s2 = 0.;
  for(int i=0;i<v.size();i++){
    s += v(i);
    s2 += v(i)*v(i);
  }
  double N = v.size();
  double mean = s/N;
  double var = s2/N - mean*mean;
  double stddev = sqrt(var);
  
  for(int i=0;i<v.size();i++) 
    v(i) = (v(i) - mean)/stddev;
}

struct LOO_SSE_Loss{
  typedef double CostType;
  typedef singleValueContainer<double> ParameterType;
  typedef singleValueContainer<double> CostDerivativeType;
  typedef NumericSquareMatrix<double> CostSecondDerivativeMatrixType;

  std::vector<Vector> xt_loo;
  std::vector<Vector> yt_loo;
  Vector xt;
  Vector yt;
  int np;
  int npm1;

  LOO_SSE_Loss(const Vector &xt, const Vector &yt): xt(xt), yt(yt){
    np = xt.size();
    npm1 = np - 1;
    xt_loo.resize(np, Vector(npm1));
    yt_loo.resize(np, Vector(npm1));
    for(int p=0;p<np;p++){
      int ii=0;
      for(int j=0;j<np;j++){
	if(j != p){
	  xt_loo[p][ii] = xt[j];
	  yt_loo[p][ii] = yt[j];
	  ii++;
	}
      }
    }
  };

  double cost(const ParameterType &fparams) const{
    double lambda = *fparams;
    double sigma = 1.;
    Vector yp(np);
    Matrix Sp(np);
    
    double out = 0.;
    for(int p=0;p<np;p++){
      GPfit(yp,Sp,xt,yt_loo[p],xt_loo[p],[](double x){ return 0.; },sigma,lambda);
      out += pow(yp(p) - yt(p), 2.);
    }
    return out;
  }
  void derivatives(CostDerivativeType &derivs, CostSecondDerivativeMatrixType &second_derivs, const ParameterType &params) const{ }

  int Nparams() const{ return 1; }
};



int main(void){
  RNG.initialize(1234);
  
  {
    //Demonstrate how random multivariate normal data with a square exponential covariance function results in a smooth curve
    Vector xpred = vrange(0.,0.01,10.);
    int npred = xpred.size();

    double sigma = 1;
    double lambda = 1;

    Vector mu(npred); mu.zero();
    Matrix Sigma(npred);
    for(int i=0;i<npred;i++)
      for(int j=0;j<npred;j++)
	Sigma(i,j) = cov(xpred[i],xpred[j],sigma,lambda) + (i==j ? 1e-10 : 0.);
  
    Vector v = randMVN(mu,Sigma);
  
    {
      MatPlotLibScriptGenerate plot;
      plot.errorBand(VectorAccessor(xpred,v));
      plot.write("randMVN.py","randMVN.pdf");
    }
  }

  {
    //Demonstrate how we can solve exactly for a smooth curve that passes through an arbitrary number of training data points
    double sigma = 1;
    double lambda = 1;

    //Training data
    Vector xt({4.5, 7, 9.5, 10}); 
    Vector yt({3, 1.5, 2.8, 3.5});

    //Positions for predictions
    Vector xp = vrange(0.,0.01,15);
  
    int ntrain = xt.size();
    int npred = xp.size();

    Matrix Stt(ntrain, [&](const int i,const int j){ return cov(xt[i],xt[j],sigma,lambda); });
    RectMatrix Spt(npred,ntrain, [&](const int i,const int j){ return cov(xp[i],xt[j],sigma,lambda); });

    Matrix Sttinv(ntrain);
    svd_inverse(Sttinv,Stt);

    Vector yp = Spt*(Sttinv*yt);

    {
      MatPlotLibScriptGenerate plot;
      plot.plotData(VectorAccessor(xt,yt));

      kwargsType kwargs;
      kwargs["color"] = 'b';
      plot.errorBand(VectorAccessor(xp,yp),kwargs);

      plot.write("GPinterp.py","GPinterp.pdf");
    }

    //We can also solve for the errors on the interpolating function
    RectMatrix Spp(npred,npred, [&](const int i,const int j){ return cov(xp[i],xp[j],sigma,lambda); });
    RectMatrix Sp = Spp - (Spt*Sttinv)*Spt.transpose();
    Vector sp(npred, [&](const int i){ return Sp(i,i)>= 0. ? sqrt(Sp(i,i)) : 0.; });
    std::cout << "sp: " << sp << std::endl;

    {
      MatPlotLibScriptGenerate plot;
      plot.plotData(VectorAccessor(xt,yt));

      kwargsType kwargs;
      kwargs["color"] = 'b';
      plot.errorBand(VectorAccessor(xp,yp),kwargs);

      kwargs["alpha"] = 0.3;
      plot.errorBand(VectorAccessorErrBand(xp,yp,sp),kwargs);
      plot.write("GPinterp_werr.py","GPinterp_werr.pdf");
    }
  }


  {
    //Repeat the above using the general fit function
    double sigma = 1;
    double lambda = 1;

    //Training data
    Vector xt({4.5, 7, 9.5, 10}); 
    Vector yt({3, 1.5, 2.8, 3.5});

    //Positions for predictions
    Vector xp = vrange(0.,0.01,15);
  
    //Normalize vectors
    normalizeVec(xt);
    normalizeVec(yt);
    normalizeVec(xp);
    
    Vector yp;
    Matrix Sp;
    GPfit(yp,Sp,xp,yt,xt,[](double x){ return 0.; },sigma,lambda);

    Vector sp(xp.size(), [&](const int i){ return Sp(i,i)>= 0. ? sqrt(Sp(i,i)) : 0.; });

    {
      MatPlotLibScriptGenerate plot;
      plot.plotData(VectorAccessor(xt,yt));

      kwargsType kwargs;
      kwargs["color"] = 'b';
      plot.errorBand(VectorAccessor(xp,yp),kwargs);

      kwargs["alpha"] = 0.3;
      plot.errorBand(VectorAccessorErrBand(xp,yp,sp),kwargs);
      plot.write("GPgen.py","GPgen.pdf");
    }
  }

  {
    //Repeat the above but try to optimize lambda for a more realistic interpolation

    //Training data
    Vector xt({4.5, 7, 9.5, 10}); 
    Vector yt({3, 1.5, 2.8, 3.5});

    //Positions for predictions
    Vector xp = vrange(0.,0.01,15);
  
    //Normalize vectors
    normalizeVec(xt);
    normalizeVec(yt);
    normalizeVec(xp);

    LOO_SSE_Loss loss(xt, yt);
    GSLmultidimMinimizerParams min_p;
    min_p.algorithm = GSLmultiminAlgorithm::NMsimplex2;
    min_p.step_size = std::vector<double>(1,0.05);
    min_p.verbose = true;
    min_p.stopping_conditions = std::vector<GSLmultiminStoppingConditionAndValue>({ {GSLmultiminStoppingCondition::StopSimplexSize, 1e-4} });

    GSLmultidimMinimizer<LOO_SSE_Loss> min(loss,min_p);
    singleValueContainer<double> opt_lambda(1.);
    double min_cost = min.fit(opt_lambda);
    std::cout << "best lambda " << *opt_lambda << std::endl;

    Vector yp;
    Matrix Sp;
    GPfit(yp,Sp,xp,yt,xt,[](double x){ return 0.; },1.,*opt_lambda);

    Vector sp(xp.size(), [&](const int i){ return Sp(i,i)>= 0. ? sqrt(Sp(i,i)) : 0.; });

    {
      MatPlotLibScriptGenerate plot;
      plot.plotData(VectorAccessor(xt,yt));

      kwargsType kwargs;
      kwargs["color"] = 'b';
      plot.errorBand(VectorAccessor(xp,yp),kwargs);

      kwargs["alpha"] = 0.3;
      plot.errorBand(VectorAccessorErrBand(xp,yp,sp),kwargs);
      plot.write("GPgen_opt.py","GPgen_opt.pdf");
    }
  }


  return 0;
}
