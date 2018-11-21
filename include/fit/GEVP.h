#ifndef _CPSFIT_GEVP_H_
#define _CPSFIT_GEVP_H_

//Generalized eigenvalue problem solver for lattice data

#include<config.h>
#include<utils/macros.h>
#include<utils/utils.h>

#include<tensors/numeric_vector.h>
#include<data_series/data_series.h>

CPSFIT_START_NAMESPACE

template<typename T>
class GEVPsolver{
  std::map< std::pair<int,int>, std::vector<T> > evals_t0_t; //only populate elements with t > t0 and t0,t < tmax providing the solver is able to compute the GEVP
  std::map< std::pair<int,int>, std::vector<NumericVector<T> > > evecs_t0_t;
  int tmax;
  bool verbose;
public:
  //Get the eigenvalues iff the element with t0,t is available, else return NULL
  std::vector<T> const* evals(const int t0, const int t) const{ 
    auto it = evals_t0_t.find({t0,t});
    return it == evals_t0_t.end() ? NULL : it->second;
  }
  //Get the eigenvectors iff the element with t0,t is available, else return NULL
  std::vector<NumericVector<T> > const* evecs(const int t0, const int t) const{ 
    auto it = evecs_t0_t.find({t0,t});
    return it == evecs_t0_t.end() ? NULL : it->second;
  }

  //Solve for a general time series of real symmetric positive-definite matrices for t0 in {0..tmax-1} and with t={t0+1..tmax}. Results that are solvable are added to the map
  //Accessor should be lambda-like returning matrix with time t, e.g.  const MatrixType & operator()(const MatrixSeriesType &s, const int t){ return s[t]; }
  template<typename MatrixSeriesType, typename Accessor>
  void solve(const MatrixSeriesType &matrix_series, const Accessor &acc, const int _tmax){
    tmax = _tmax;
    for(int t0=0; t0<_tmax-1;t0++){
      for(int t=t0+1; t<_tmax; t++){
	bool fail = false;
	std::vector<NumericVector<T> > evecs_t;
	std::vector<T> evals_t;
	try{
	  auto A = acc(matrix_series,t);
	  auto B = acc(matrix_series,t0);
	  auto resid = symmetricMatrixGEVPsolve(evecs_t, evals_t,A,B);
	}catch(const std::exception & ex){
	  if(verbose) std::cout << t0 << " " << t << " GEVP solver fail with error \"" << ex.what() << "\" (\"input domain error\" usually means one or both matrices are not positive definite), skipping\n";
	  fail = true;
	}
	if(!fail){
	  evals_t0_t[{t0,t}] = std::move(evals_t);
	  evecs_t0_t[{t0,t}] = std::move(evecs_t);
	}
      }
    }
  }
  template<typename _GeneralizedCoordinate, typename _MatrixType, template<typename,typename> class _PairType = std::pair>
  inline void solve(const dataSeries<_GeneralizedCoordinate,_MatrixType,_PairType> &matrix_series, const int _tmax){
    solve(matrix_series, [&](const dataSeries<_GeneralizedCoordinate,_MatrixType,_PairType> &s, const int t){ return matrix_series.value(t); }, _tmax);
  }
  
  //Compute the effective energies for time t. Vector size is zero if failed to compute
  std::vector<T> effectiveEnergy(const int t0, const int t) const{
    auto it_tp1 = evals_t0_t.find({t0,t+1});
    auto it_t = evals_t0_t.find({t0,t});
    if(it_tp1 != evals_t0_t.end() && it_t != evals_t0_t.end()){ 
      const int N = it_t->second.size();
      std::vector<T> E_eff(N);
      for(int n=0;n<N;n++){
	E_eff[n] = log( it_t->second[n] ) - log( it_tp1->second[n] );
      }      
      return E_eff;
    }else{
      return std::vector<T>();
    }
  }
  // template<typename MatrixSeriesType, typename Accessor>
  // std::vector<T> effectiveAmplitude(const int t0, const int t, const MatrixSeriesType &matrix_series, const Accessor &acc){) const{
  //   //\sum_j C_ij(t)u^n_j = e^{-E_n t}\psi_i^n     where \psi_i^n = A_ni  is the amplitude of operator i with state n
  //   //u^n_j \propto v(t,t0)^n_j  i.e.   u^n_j = \alpha^n(t,t0)v_j^n(t,t0)
  //   //Thus   \sum_j \alpha^n(t,t0) C_ij(t)v_j^n(t,t0) = e^{-E_n t}\psi_i^n 
  //   //Define \chi_i^n(t,t0) = \psi_i^n / \alpha^n(t,t0) = e^{E_n t}\sum_j C_ij(t)v_j^n(t,t0)
  //   //We also have C_ij(t) = \sum_n e^{-E_n t} \psi_i^n \psi_j^n = \sum_n e^{-E_n t} [\alpha^n(t,t0)]^2 \chi_i^n(t,t0) \chi_j^n(t,t0)

  GEVPsolver(bool verbose = false): tmax(0), verbose(verbose){}
};

CPSFIT_END_NAMESPACE
#endif
