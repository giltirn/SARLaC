#ifndef _SARLAC_GEVP_H_
#define _SARLAC_GEVP_H_

//Generalized eigenvalue problem solver for lattice data

#include<config.h>
#include<utils/macros.h>
#include<utils/utils.h>
#include<serialize/hdf5_serialize.h>

#include<tensors/numeric_vector.h>
#include<distribution/jackknife.h>
#include<data_series/data_series.h>

SARLAC_START_NAMESPACE

template<typename T>
class GEVPsolverBase{
public:
  typedef T DistributionType;

protected:
  typedef std::pair<int,int> t0_t_pair;

  std::map< t0_t_pair, std::vector<T> > evals_t0_t; //only populate elements with t > t0 and t0,t < tmax providing the solver is able to compute the GEVP
  std::map< t0_t_pair, std::vector<NumericVector<T> > > evecs_t0_t;
  std::map< t0_t_pair, std::vector<T> > resid_t0_t;
  int tmax;
  bool verbose;
public:
  //Get the eigenvalues iff the element with t0,t is available, else return NULL
  std::vector<T> const* evals(const int t0, const int t) const{ 
    auto it = evals_t0_t.find({t0,t});
    return it == evals_t0_t.end() ? NULL : &it->second;
  }
  std::vector<T> * evals(const int t0, const int t){ 
    auto it = evals_t0_t.find({t0,t});
    return it == evals_t0_t.end() ? NULL : &it->second;
  }

  //Get the eigenvectors iff the element with t0,t is available, else return NULL
  //Eigenvector index is equivalent to the state index, and the coordinates are *operator* indices
  std::vector<NumericVector<T> > const* evecs(const int t0, const int t) const{ 
    auto it = evecs_t0_t.find({t0,t});
    return it == evecs_t0_t.end() ? NULL : &it->second;
  }
  std::vector<NumericVector<T> > * evecs(const int t0, const int t){ 
    auto it = evecs_t0_t.find({t0,t});
    return it == evecs_t0_t.end() ? NULL : &it->second;
  }

  //GEVP residuals
  std::vector<T> const* residuals(const int t0, const int t) const{ 
    auto it = resid_t0_t.find({t0,t});
    return it == resid_t0_t.end() ? NULL : &it->second;
  }

  //Solve for a general time series of real symmetric positive-definite matrices for t0 in {t0_start..t0_end} and with t={t0+1..tmax}. Results that are solvable are added to the map
  //Accessor should be lambda-like returning matrix with time t, e.g.  const MatrixType & operator()(const MatrixSeriesType &s, const int t){ return s[t]; }
  template<typename MatrixSeriesType, typename Accessor>
  void solve(const MatrixSeriesType &matrix_series, const Accessor &acc, const int _tmax, const int t0_start, const int t0_end){
    tmax = _tmax;
    for(int t0=t0_start; t0<=t0_end;t0++){
      for(int t=t0+1; t<=_tmax; t++){
	bool fail = false;
	std::vector<NumericVector<T> > evecs_t;
	std::vector<T> evals_t;	
	std::vector<T> resid_t;

	try{
	  auto A = acc(matrix_series,t);
	  auto B = acc(matrix_series,t0);
	  resid_t = symmetricMatrixGEVPsolve(evecs_t, evals_t,A,B);
	}catch(const std::exception & ex){
	  if(verbose) std::cout << "t0=" << t0 << " t=" << t << " GEVP solver fail with error \"" << ex.what() << "\" (\"input domain error\" usually means one or both matrices are not positive definite), skipping\n";
	  fail = true;
	}

	if(!fail){
	  evals_t0_t[{t0,t}] = std::move(evals_t);
	  evecs_t0_t[{t0,t}] = std::move(evecs_t);
	  resid_t0_t[{t0,t}] = std::move(resid_t);
	}
      }
    }
  }
  template<typename _GeneralizedCoordinate, typename _MatrixType, template<typename,typename> class _PairType = std::pair>
  inline void solve(const dataSeries<_GeneralizedCoordinate,_MatrixType,_PairType> &matrix_series, const int _tmax, const int t0_start, const int t0_end){
    solve(matrix_series, [&](const dataSeries<_GeneralizedCoordinate,_MatrixType,_PairType> &s, const int t){ return matrix_series.value(t); }, _tmax, t0_start, t0_end);
  }
  

  GEVPsolverBase(bool verbose = false): tmax(0), verbose(verbose){}
  GEVPsolverBase(const GEVPsolverBase &c) = default;
  GEVPsolverBase(GEVPsolverBase &&c) = default;

  void write(HDF5writer &wr, const std::string &tag) const{
    wr.enter(tag);
    SARLaC::write(wr, evals_t0_t, "evals_t0_t");
    SARLaC::write(wr, evecs_t0_t, "evecs_t0_t");
    SARLaC::write(wr, resid_t0_t, "resid_t0_t");
    SARLaC::write(wr, tmax, "tmax");
    wr.leave();
  }
  void write(const std::string &filename, const std::string &tag = "GEVP_solutions"){
    HDF5writer wr(filename);
    write(wr, tag);
  }

  void read(HDF5reader &rd, const std::string &tag){
    rd.enter(tag);
    SARLaC::read(rd, evals_t0_t, "evals_t0_t");
    SARLaC::read(rd, evecs_t0_t, "evecs_t0_t");
    SARLaC::read(rd, resid_t0_t, "resid_t0_t");
    SARLaC::read(rd, tmax, "tmax");
    rd.leave();
  }
  void read(const std::string &filename, const std::string &tag = "GEVP_solutions"){
    HDF5reader rd(filename);
    read(rd, tag);
  }

};


//GEVP for standard correlators C(t)
template<typename T>
class GEVPsolver: public GEVPsolverBase<T>{
public:
  GEVPsolver(bool verbose = false):  GEVPsolverBase<T>(verbose){}
  GEVPsolver(const GEVPsolver &c) = default;
  GEVPsolver(GEVPsolver &&c) = default;

  //Compute the effective energies for time t. Vector size is zero if failed to compute
  std::vector<T> effectiveEnergy(const int t0, const int t) const{
    std::vector<T> const* evals_tp1 = this->evals(t0,t+1);
    std::vector<T> const* evals_t = this->evals(t0,t);

    if(evals_tp1 != NULL && evals_t != NULL){
      const int N = evals_t->size();
      std::vector<T> E_eff(N);
      for(int n=0;n<N;n++){
	E_eff[n] = log( (*evals_t)[n] ) - log( (*evals_tp1)[n] );
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

    std::vector<NumericVector<T> > const* evecs_t0_t = this->evecs(t0,t);
    std::vector<T> const* evals_t0_t0ptd2 = this->evals(t0,t0+t/2);
    std::vector<T> const* evals_t0_t0pt = this->evals(t0,t0+t);
    std::vector<T> Eeff = this->effectiveEnergy(t0,t); //indexed in state

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

};



//For C(t) - C(t+1).  Definition of the effective energy is the same but that of the amplitude will need to be modified (currently not implemented)
template<typename T>
class GEVPsubNeighborTslice: public GEVPsolverBase<T>{
public:
  GEVPsubNeighborTslice(bool verbose = false):  GEVPsolverBase<T>(verbose){}
  GEVPsubNeighborTslice(const GEVPsubNeighborTslice &c) = default;
  GEVPsubNeighborTslice(GEVPsubNeighborTslice &&c) = default;

  //Compute the effective energies for time t. Vector size is zero if failed to compute
  std::vector<T> effectiveEnergy(const int t0, const int t) const{
    std::vector<T> const* evals_tp1 = this->evals(t0,t+1);
    std::vector<T> const* evals_t = this->evals(t0,t);

    if(evals_tp1 != NULL && evals_t != NULL){
      const int N = evals_t->size();
      std::vector<T> E_eff(N);
      for(int n=0;n<N;n++){
	E_eff[n] = log( (*evals_t)[n] ) - log( (*evals_tp1)[n] );
      }      
      return E_eff;
    }else{
      return std::vector<T>();
    }
  }

  //Not yet implemented
  template<typename MatrixSeriesType, typename Accessor>
  std::vector<std::vector<T> > effectiveAmplitude(const int t0, const int t, const MatrixSeriesType &matrix_series, const Accessor &acc) const{
    return std::vector<std::vector<T> >();
  }
  template<typename _GeneralizedCoordinate, typename _MatrixType, template<typename,typename> class _PairType = std::pair>
  std::vector<std::vector<T> > effectiveAmplitude(const int t0, const int t, const dataSeries<_GeneralizedCoordinate,_MatrixType,_PairType> &matrix_series) const{
    return effectiveAmplitude(t0,t, matrix_series, [&](const dataSeries<_GeneralizedCoordinate,_MatrixType,_PairType> &s, const int t){ return matrix_series.value(t); });
  }

};


//For C(t) - C(tsub)  for fixed tsub. Requires tsub > t > t0 
class GEVPsubFixedTslice: public GEVPsolverBase<jackknifeDistribution<double> >{
  typedef jackknifeDistribution<double> jackknifeDistributionD;
  int tsub;
  
  class TsubEnergySolveFunc{
  public:
    typedef double ValueType;
    typedef singleValueContainer<double> ParameterType;
    typedef singleValueContainer<double> ValueDerivativeType;
    typedef double GeneralizedCoordinate; 

    int t;
    int t0;
    int tsub;

    TsubEnergySolveFunc(int t0, int t, int tsub): t0(t0), t(t), tsub(tsub){}

    template<typename T, typename Boost>
    T eval(const GeneralizedCoordinate &c, const ParameterType &p, const Boost &boost) const{
      T E = boost(*p, 0);
      T num = exp(-E*t) - exp(-E*tsub);
      T den = exp(-E*t0) - exp(-E*tsub);
      T f= num/den;
      return f;
    }

    inline ValueType value(const GeneralizedCoordinate &x, const ParameterType &p) const{
      return eval<double>(x,p,[&](const double v, const int i){ return v; });
    }
    inline ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &x, const ParameterType &p) const{
      return ValueDerivativeType(eval<dual>(x,p,[&](const double v, const int j){ return dual(v, j==0 ? 1.:0.); }).xp);
    }
  
    static inline int Nparams(){ return 1; }
  };

  double solveEffectiveEnergyFixedTsub(double lambda, const int t0, const int t, bool verbose = false) const{
    typedef correlationFunction<double, double> DataType;
    DataType d(1); d.coord(0) = 0; d.value(0) = lambda;
    
    typedef UncorrelatedChisqCostFunction<TsubEnergySolveFunc,DataType> CostFunction;
    TsubEnergySolveFunc fitfunc(t0,t,tsub);
    CostFunction cost(fitfunc, d, std::vector<double>(1,1.0));
    MarquardtLevenbergParameters<double> mlp;
    //mlp.verbose = verbose;
    mlp.exit_on_convergence_fail = false;
    mlp.delta_cost_min = 1e-10;
    mlp.max_iter = 100;
    MarquardtLevenbergMinimizer<CostFunction> min(cost, mlp);
    singleValueContainer<double> E(0.3);
    
    //if(verbose) std::cout << "Starting fit for "  << "t0=" << t0 << " t=" << t << " tsub=" << tsub << " to obtain lambda=" << lambda << std::endl;
    
    double chisq = min.fit(E);
    
    double fvalue = fitfunc.value(0., E);
    
    if(!min.hasConverged()){
      if(verbose) std::cout << "t0=" << t0 << " t=" << t << " tsub=" << tsub<< " did not converge in search for solution of lambda=" << lambda << std::endl;
      return nan("");
    }else if( fabs( (fvalue - lambda)/lambda ) > 1e-4 ){
      //if(verbose) 
      //std::cout << "t0=" << t0 << " t=" << t << " tsub=" << tsub<< " could not find solution: expect=" <<  lambda << " best= " << fvalue << " at E=" << *E << std::endl;
      return nan("");
    }else{
      //if(verbose) std::cout << "t0=" << t0 << " t=" << t << " tsub=" << tsub<< " inverse E=" << *E << " lambda(E)= " << fvalue << " expect " << lambda << std::endl;
      return *E;
    }
  }

public:
  GEVPsubFixedTslice(const int tsub, bool verbose = false):  tsub(tsub), GEVPsolverBase<jackknifeDistributionD>(verbose){}
  GEVPsubFixedTslice(const GEVPsubFixedTslice &c) = default;
  GEVPsubFixedTslice(GEVPsubFixedTslice &&c) = default;

  //Compute the effective energies for time t. Vector size is zero if failed to compute
  std::vector<jackknifeDistributionD> effectiveEnergy(const int t0, const int t) const{
    if(t0 == tsub || t == tsub || t0 == t) return std::vector<jackknifeDistributionD>();

    std::vector<jackknifeDistributionD> const* eval = this->evals(t0,t);
    if(eval!=NULL){ 
      const int N = eval->size();
      const int nsample = (*eval)[0].size();
      
      std::vector<jackknifeDistributionD> E_eff(N, jackknifeDistributionD(nsample));
      for(int n=0;n<N;n++)
	for(int s=0;s<nsample;s++)
	  E_eff[n].sample(s) = solveEffectiveEnergyFixedTsub( (*eval)[n].sample(s),t0,t,s==0);
      
      return E_eff;
    }else{
      return std::vector<jackknifeDistributionD>();
    }
  }

  //Not yet implemented
  template<typename MatrixSeriesType, typename Accessor>
  std::vector<std::vector<jackknifeDistributionD> > effectiveAmplitude(const int t0, const int t, const MatrixSeriesType &matrix_series, const Accessor &acc) const{
    return std::vector<std::vector<jackknifeDistributionD> >();
  }
  template<typename _GeneralizedCoordinate, typename _MatrixType, template<typename,typename> class _PairType = std::pair>
  std::vector<std::vector<jackknifeDistributionD> > effectiveAmplitude(const int t0, const int t, const dataSeries<_GeneralizedCoordinate,_MatrixType,_PairType> &matrix_series) const{
    return effectiveAmplitude(t0,t, matrix_series, [&](const dataSeries<_GeneralizedCoordinate,_MatrixType,_PairType> &s, const int t){ return matrix_series.value(t); });
  }
  
};




//This more relaxed algorithm will work for non-symmetric and/or non positive-definite matrices
//The eigenvalues and eigenvectors will in general be complex
template<typename T>
class GEVPnonSymmSolverBase{
protected:
  typedef std::pair<int,int> t0_t_pair;
  typedef Complexify<T> complex_T;

  std::map< t0_t_pair, std::vector<complex_T> > evals_t0_t; //only populate elements with t > t0 and t0,t < tmax providing the solver is able to compute the GEVP
  std::map< t0_t_pair, std::vector<NumericVector<complex_T> > > evecs_t0_t;
  std::map< t0_t_pair, std::vector<T> > resid_t0_t;
  int tmax;
  bool verbose;
public:
  typedef T DataType;

  //Get the eigenvalues iff the element with t0,t is available, else return NULL
  std::vector<complex_T> const* evals(const int t0, const int t) const{ 
    auto it = evals_t0_t.find({t0,t});
    return it == evals_t0_t.end() ? NULL : &it->second;
  }
  std::vector<complex_T> * evals(const int t0, const int t){ 
    auto it = evals_t0_t.find({t0,t});
    return it == evals_t0_t.end() ? NULL : &it->second;
  }

  //Get the eigenvectors iff the element with t0,t is available, else return NULL
  //Eigenvector index is equivalent to the state index, and the coordinates are *operator* indices
  std::vector<NumericVector<complex_T> > const* evecs(const int t0, const int t) const{ 
    auto it = evecs_t0_t.find({t0,t});
    return it == evecs_t0_t.end() ? NULL : &it->second;
  }
  std::vector<NumericVector<complex_T> > * evecs(const int t0, const int t){ 
    auto it = evecs_t0_t.find({t0,t});
    return it == evecs_t0_t.end() ? NULL : &it->second;
  }

  //GEVP residuals
  std::vector<T> const* residuals(const int t0, const int t) const{ 
    auto it = resid_t0_t.find({t0,t});
    return it == resid_t0_t.end() ? NULL : &it->second;
  }

  //Solve for a general time series of real symmetric positive-definite matrices for t0 in {t0_start..t0_end} and with t={t0+1..tmax}. Results that are solvable are added to the map
  //Accessor should be lambda-like returning matrix with time t, e.g.  const MatrixType & operator()(const MatrixSeriesType &s, const int t){ return s[t]; }
  template<typename MatrixSeriesType, typename Accessor>
  void solve(const MatrixSeriesType &matrix_series, const Accessor &acc, const int _tmax, const int t0_start, int t0_end){
    tmax = _tmax;
    for(int t0=t0_start; t0<=t0_end;t0++){
      for(int t=t0+1; t<=_tmax; t++){
	bool fail = false;
	std::vector<NumericVector<complex_T> > evecs_t;
	std::vector<complex_T> evals_t;	
	std::vector<T> resid_t;
	try{
	  auto A = acc(matrix_series,t);
	  auto B = acc(matrix_series,t0);
	  resid_t = nonSymmetricMatrixGEVPsolve(evecs_t, evals_t,A,B);
	}catch(const std::exception & ex){
	  if(verbose) std::cout << "t0=" << t0 << " t=" << t << " GEVP solver fail with error \"" << ex.what() << "\"\n";
	  fail = true;
	}
	if(!fail){
	  evals_t0_t[{t0,t}] = std::move(evals_t);
	  evecs_t0_t[{t0,t}] = std::move(evecs_t);
	  resid_t0_t[{t0,t}] = std::move(resid_t);
	}
      }
    }
  }
  template<typename _GeneralizedCoordinate, typename _MatrixType, template<typename,typename> class _PairType = std::pair>
  inline void solve(const dataSeries<_GeneralizedCoordinate,_MatrixType,_PairType> &matrix_series, const int _tmax, const int t0_start, int t0_end){
    solve(matrix_series, [&](const dataSeries<_GeneralizedCoordinate,_MatrixType,_PairType> &s, const int t){ return matrix_series.value(t); }, _tmax, t0_start, t0_end);
  }
  

  explicit GEVPnonSymmSolverBase(bool verbose = false): tmax(0), verbose(verbose){}
  GEVPnonSymmSolverBase(const GEVPnonSymmSolverBase &c) = default;
  GEVPnonSymmSolverBase(GEVPnonSymmSolverBase &&c) = default;

  void write(HDF5writer &wr, const std::string &tag) const{
    wr.enter(tag);
    SARLaC::write(wr, evals_t0_t, "evals_t0_t");
    SARLaC::write(wr, evecs_t0_t, "evecs_t0_t");
    SARLaC::write(wr, resid_t0_t, "resid_t0_t");
    SARLaC::write(wr, tmax, "tmax");
    wr.leave();
  }
  void write(const std::string &filename, const std::string &tag = "GEVP_solutions"){
    HDF5writer wr(filename);
    write(wr, tag);
  }

  void read(HDF5reader &rd, const std::string &tag){
    rd.enter(tag);
    SARLaC::read(rd, evals_t0_t, "evals_t0_t");
    SARLaC::read(rd, evecs_t0_t, "evecs_t0_t");
    SARLaC::read(rd, resid_t0_t, "resid_t0_t");
    SARLaC::read(rd, tmax, "tmax");
    rd.leave();
  }
  void read(const std::string &filename, const std::string &tag = "GEVP_solutions"){
    HDF5reader rd(filename);
    read(rd, tag);
  }

};


SARLAC_END_NAMESPACE
#endif
