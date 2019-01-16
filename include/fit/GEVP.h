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
    return it == evals_t0_t.end() ? NULL : &it->second;
  }
  //Get the eigenvectors iff the element with t0,t is available, else return NULL
  //Eigenvector index is equivalent to the state index, and the coordinates are *operator* indices
  std::vector<NumericVector<T> > const* evecs(const int t0, const int t) const{ 
    auto it = evecs_t0_t.find({t0,t});
    return it == evecs_t0_t.end() ? NULL : &it->second;
  }

  //Solve for a general time series of real symmetric positive-definite matrices for t0 in {0..tmax-1} and with t={t0+1..tmax}. Results that are solvable are added to the map
  //Accessor should be lambda-like returning matrix with time t, e.g.  const MatrixType & operator()(const MatrixSeriesType &s, const int t){ return s[t]; }
  template<typename MatrixSeriesType, typename Accessor>
  void solve(const MatrixSeriesType &matrix_series, const Accessor &acc, const int _tmax){
    tmax = _tmax;
    for(int t0=0; t0<=_tmax-1;t0++){
      for(int t=t0+1; t<=_tmax; t++){
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

  //Indexing of output matrix is [operator index][state index]
  //Note: Requires evecs at (t0,t) and evals at  (t0+t/2,t0)  (t0+t,t0)  - returns empty output if these are not available 
  template<typename MatrixSeriesType, typename Accessor>
  std::vector<std::vector<T> > effectiveAmplitude(const int t0, const int t, const MatrixSeriesType &matrix_series, const Accessor &acc) const{
    //This is the operator A_eff from Eq. 2.15-2.17 of http://arxiv.org/pdf/0902.1265.pdf    between state <n| and |0>
    //Note that this paper put the later time (t as opposed to t0) as the row index whereas in this code it is the column index
    //As we have only values at integer times we must make a small adjustment to 2.17:
    //Given that 
    //          v(t0,t) is independent of t (iff no excited state contamination)
    //          C(t) ~ exp(-Et)   [2.5]    
    //          lambda(t0,t) ~ exp(-E(t-t0))
    //Thus for even t=2k   (v(t0,2k), C(2k)v(t0,2k))^{-1/2}  ~ ( exp(-2Ek) )^{-1/2} ~  exp(Ek)
    //lambda(t0, t0+(2k)/2)/lambda(t0,t0+2k) ~ exp(-Ek)/exp(-2Ek) ~  exp(Ek)
    //and the time dependence of Rn ~ exp(2Ek)
    
    //However consider odd times  t=2k+1   and use integer division
    //(v(t0,2k+1), C(2k+1)v(t0,2k+1))^{-1/2}  ~ ( exp(-E (2k+1) ) )^{-1/2} ~ sqrt( exp(E (2k+1) ) )  ~ exp(Ek)sqrt(exp(E)) ~ exp(E (k+1))
    //lambda(t0, t0+(2k+1)/2)/lambda(t0,t0+2k+1) = lambda(t0, t0+k)/lambda(t0,t0+2k+1) ~ exp(-Ek)/exp(-E(2k+1)) ~ exp(Ek)exp(E)
    //and the time dependence of Rn ~ exp(2Ek) exp(E)^{3/2} ~ exp(E (2k+1)) exp(E)^{1/2}

    //Thus for odd timeslices we must remove an additional factor of exp(E)^{1/2}

    std::vector<NumericVector<T> > const* evecs_t0_t = evecs(t0,t);
    std::vector<T> const* evals_t0_t0ptd2 = evals(t0,t0+t/2);
    std::vector<T> const* evals_t0_t0pt = evals(t0,t0+t);
    std::vector<T> Eeff = effectiveEnergy(t0,t); //indexed in state

    std::vector<std::vector<T> > out;
    if(evecs_t0_t == NULL || evals_t0_t0ptd2 == NULL || evals_t0_t0pt == NULL || Eeff.size() == 0) return out; //return empty vector if not possible

    T zero = Eeff[0]; zeroit(zero);

    int nstate = evecs_t0_t->size();
    int nop = nstate;
    
    //Construct a matrix from the eigenvector with row index the state and column index the operator
    NumericSquareMatrix<T> v_n_i(nstate);
    for(int n=0;n<nstate;n++)
      for(int i=0;i<nop;i++)
	v_n_i(n,i) = (*evecs_t0_t)[n](i);
    NumericSquareMatrix<T> v_n_i_inv(v_n_i);

    svd_inverse(v_n_i_inv, v_n_i); //inverse has row idx = operator index and column idx = state index
    
    out.resize(nop, std::vector<T>(nstate));

    for(int n=0;n<nstate;n++){ //loop over state index
      //Construct inner product
      //( v_n,  C(t)v_n )
      
      std::vector<T> Ctvn(nop, zero); 
      for(int i=0;i<nop;i++)
	for(int j=0;j<nop;j++)
	  Ctvn[i] = Ctvn[i] + acc(matrix_series,t)(i,j) *  (*evecs_t0_t)[n](j);
      
      T v_dot_Ctvn(zero);
      for(int i=0;i<nop;i++)
	v_dot_Ctvn = v_dot_Ctvn + (*evecs_t0_t)[n](i) * Ctvn[i];
	
      T Rn = (*evals_t0_t0ptd2)[n]/(*evals_t0_t0pt)[n]/sqrt(v_dot_Ctvn);
      
      //Apply correction for odd t
      if(t % 2 == 1){
	Rn = Rn / sqrt(exp(Eeff[n]));
      }
      for(int i=0;i<nop;i++){
	out[i][n] = v_n_i_inv(i,n) * exp(Eeff[n] * t) / Rn;
      }
    }
    return out;
  }

  template<typename _GeneralizedCoordinate, typename _MatrixType, template<typename,typename> class _PairType = std::pair>
  std::vector<std::vector<T> > effectiveAmplitude(const int t0, const int t, const dataSeries<_GeneralizedCoordinate,_MatrixType,_PairType> &matrix_series) const{
    return effectiveAmplitude(t0,t, matrix_series, [&](const dataSeries<_GeneralizedCoordinate,_MatrixType,_PairType> &s, const int t){ return matrix_series.value(t); });
  }

  GEVPsolver(bool verbose = false): tmax(0), verbose(verbose){}
};

CPSFIT_END_NAMESPACE
#endif
