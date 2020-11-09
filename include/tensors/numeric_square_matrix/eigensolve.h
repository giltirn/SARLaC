#ifndef _CPSFIT_NUMERIC_SQUARE_MATRIX_EIGENSOLVE_H_
#define _CPSFIT_NUMERIC_SQUARE_MATRIX_EIGENSOLVE_H_

#include<mutex>
#include<distribution/distribution_iterate.h>
#include<tensors/gsl_eigensolve.h>
#include<tensors/numeric_vector.h>
#include<tensors/numeric_square_matrix/class.h>
#include<utils/template_wizardry/complexify.h>

CPSFIT_START_NAMESPACE

template<typename T, int is_distribution>
struct _SquareMatrixEigensolve{};

//For regular matrices
template<typename T>
struct _SquareMatrixEigensolve<T,0>{
  //Symmetric matrix solve
  inline static std::vector<double> solveSymm(std::vector<NumericVector<T> > &evecs, std::vector<T> &evals, const NumericSquareMatrix<T> &A, bool sort = true){
    return GSLsymmEigenSolver<NumericVector<T>, NumericSquareMatrix<T> >::symmetricMatrixSolve(evecs,evals,A,sort);
  }
  //Non-symmetric matrix solve
  inline static std::vector<double> solveNonSymm(std::vector<NumericVector<std::complex<T> > > &evecs, std::vector<std::complex<T> > &evals, const NumericSquareMatrix<T> &A, bool sort = true){
    return GSLnonSymmEigenSolver<NumericVector<std::complex<T> >, NumericSquareMatrix<T> >::nonSymmetricMatrixSolve(evecs,evals,A,sort);
  }

  //Symmetric matrix GEVP solve
  inline static std::vector<double> solveSymmGEVP(std::vector<NumericVector<T> > &evecs, std::vector<T> &evals, const NumericSquareMatrix<T> &A, const NumericSquareMatrix<T> &B, bool sort = true){
    return GSLsymmEigenSolver<NumericVector<T>, NumericSquareMatrix<T> >::symmetricGEVPsolveGen(evecs,evals,A,B,sort);
  }

  //NonSymmetric matrix GEVP solve
  inline static std::vector<double> solveNonSymmGEVP(std::vector<NumericVector<std::complex<T> > > &evecs, std::vector<std::complex<T> > &evals, const NumericSquareMatrix<T> &A, const NumericSquareMatrix<T> &B, bool sort = true){
    return GSLnonSymmEigenSolver<NumericVector<T>, NumericSquareMatrix<T> >::nonSymmetricGEVPsolve(evecs,evals,A,B,sort);
  }
};


//For matrices of distributions
template<typename T>
struct _SquareMatrixEigensolve<T,1>{
  static std::vector<T> solveSymm(std::vector<NumericVector<T> > &evecs, std::vector<T> &evals, const NumericSquareMatrix<T> &A, bool sort = true){
    assert(A.size() > 0);
    const int size = A.size();
    const int nsample = A(0,0).size();

    for(int i=0;i<size;i++){
      evals[i].resize(nsample);
      for(int j=0;j<size;j++)
	evecs[i][j].resize(nsample);
    }
    
    typedef typename iterate<T>::type type;
    NumericSquareMatrix<type> A_s(size);
    std::vector<NumericVector<type> > evecs_s(size, NumericVector<type>(size));
    std::vector<type> evals_s(size);

    std::vector<T> residuals(size, T(nsample));
    
    const int nit = iterate<T>::size(A(0,0));
    for(int s=0;s<nit;s++){
      for(int i=0;i<size;i++)
	for(int j=0;j<size;j++)
	  A_s(i,j) = iterate<T>::at(s, A(i,j));
	
      std::vector<double> r = GSLsymmEigenSolver<NumericVector<type>, NumericSquareMatrix<type> >::symmetricMatrixSolve(evecs_s,evals_s,A_s, sort);

      for(int i=0;i<size;i++){
	iterate<T>::at(s, residuals[i]) = r[i];
	iterate<T>::at(s, evals[i]) = evals_s[i];
	
	for(int j=0;j<size;j++)
	  iterate<T>::at(s, evecs[i](j)) = evecs_s[i](j);
      }
    }
    return residuals;
  }

  static std::vector<T> solveNonSymm(std::vector<NumericVector<Complexify<T> > > &evecs, 
				     std::vector<Complexify<T> > &evals, 
				     const NumericSquareMatrix<T> &A, bool sort = true){
    assert(A.size() > 0);
    const int size = A.size();
    const int nsample = A(0,0).size();

    for(int i=0;i<size;i++){
      evals[i].resize(nsample);
      for(int j=0;j<size;j++)
	evecs[i][j].resize(nsample);
    }
    
    typedef Complexify<T> TC;
    typedef typename iterate<T>::type type;
    typedef std::complex<type> complex_type;

    NumericSquareMatrix<type> A_s(size);
    std::vector<NumericVector<complex_type> > evecs_s(size, NumericVector<complex_type>(size));
    std::vector<complex_type> evals_s(size);

    std::vector<T> residuals(size, T(nsample));
    
    const int nit = iterate<T>::size(A(0,0));
    for(int s=0;s<nit;s++){
      for(int i=0;i<size;i++)
	for(int j=0;j<size;j++)
	  A_s(i,j) = iterate<T>::at(s, A(i,j));
	
      std::vector<double> r = GSLnonSymmEigenSolver<NumericVector<complex_type>, NumericSquareMatrix<type> >::nonSymmetricMatrixSolve(evecs_s,evals_s,A_s, sort);

      for(int i=0;i<size;i++){
	iterate<T>::at(s, residuals[i]) = r[i];
	iterate<TC>::at(s, evals[i]) = evals_s[i];
	
	for(int j=0;j<size;j++)
	  iterate<TC>::at(s, evecs[i](j)) = evecs_s[i](j);
      }
    }
    return residuals;
  }




  //GEVP: solve Av = \lambda B v
  //evals = \lambda
  //evecs = v
  static std::vector<T> solveSymmGEVP(std::vector<NumericVector<T> > &evecs, 
				      std::vector<T> &evals, 
				      const NumericSquareMatrix<T> &A, 
				      const NumericSquareMatrix<T> &B, 
				      bool sort = true){
    assert(A.size() > 0);
    const int size = A.size();
    auto dist_init = A(0,0).getInitializer();
    assert(B.size() == size && B(0,0).getInitializer() == dist_init);

    for(int i=0;i<size;i++){
      evals[i].resize(dist_init);
      for(int j=0;j<size;j++)
	evecs[i][j].resize(dist_init);
    }
    
    typedef typename iterate<T>::type type;
    std::vector<T> residuals(size, T(dist_init));
    
    const int nit = iterate<T>::size(A(0,0));
    std::mutex m;
    bool exception = false;
    std::string err_str;
    
#pragma omp parallel for
    for(int s=0;s<nit;s++){ //sample loop
      if(exception) continue;

      NumericSquareMatrix<type> A_s(size);
      NumericSquareMatrix<type> B_s(size);
      std::vector<NumericVector<type> > evecs_s(size, NumericVector<type>(size));
      std::vector<type> evals_s(size);

      for(int i=0;i<size;i++)
	for(int j=0;j<size;j++){
	  A_s(i,j) = iterate<T>::at(s, A(i,j));
	  B_s(i,j) = iterate<T>::at(s, B(i,j));
	}
      std::vector<double> r;
      try{
	r = GSLsymmEigenSolver<NumericVector<type>, NumericSquareMatrix<type> >::symmetricGEVPsolveGen(evecs_s,evals_s,A_s,B_s, sort); //doesn't require positive-definite B matrix
      }catch(const std::exception &exc){
	std::lock_guard<std::mutex> lock(m);
	exception = true;
	err_str = exc.what();
	continue;
      }

      for(int i=0;i<size;i++){
	iterate<T>::at(s, residuals[i]) = r[i];
	iterate<T>::at(s, evals[i]) = evals_s[i];
	
	for(int j=0;j<size;j++)
	  iterate<T>::at(s, evecs[i](j)) = evecs_s[i](j);
      }
    }
    if(exception) throw std::runtime_error(err_str);
    
    return residuals;
  }



  //GEVP: solve Av = \lambda B v for non-symmetric A,B
  //evals = \lambda
  //evecs = v
  static std::vector<T> solveNonSymmGEVP(std::vector<NumericVector<Complexify<T> > > &evecs, 
					 std::vector<Complexify<T> > &evals, 
					 const NumericSquareMatrix<T> &A, 
					 const NumericSquareMatrix<T> &B, 
					 bool sort = true){
    assert(A.size() > 0);
    const int size = A.size();
    auto dist_init = A(0,0).getInitializer();
    assert(B.size() == size && B(0,0).getInitializer() == dist_init);

    for(int i=0;i<size;i++){
      evals[i].resize(dist_init);
      for(int j=0;j<size;j++)
	evecs[i][j].resize(dist_init);
    }
    
    typedef typename iterate<T>::type type;
    std::vector<T> residuals(size, T(dist_init));
    
    const int nit = iterate<T>::size(A(0,0));
    std::mutex m;
    bool exception = false;
    std::string err_str;
    
    typedef Complexify<T> TC; //complex distribution type
    typedef Complexify<type> complex_type; //complex base data type

#pragma omp parallel for
    for(int s=0;s<nit;s++){ //sample loop
      if(exception) continue;

      NumericSquareMatrix<type> A_s(size);
      NumericSquareMatrix<type> B_s(size);
      std::vector<NumericVector<complex_type> > evecs_s(size, NumericVector<complex_type>(size));
      std::vector<complex_type> evals_s(size);

      for(int i=0;i<size;i++)
	for(int j=0;j<size;j++){
	  A_s(i,j) = iterate<T>::at(s, A(i,j));
	  B_s(i,j) = iterate<T>::at(s, B(i,j));
	}
      std::vector<double> r;
      try{
	r = GSLnonSymmEigenSolver<NumericVector<complex_type>, NumericSquareMatrix<type> >::nonSymmetricGEVPsolve(evecs_s,evals_s,A_s,B_s, sort);
      }catch(const std::exception &exc){
	std::lock_guard<std::mutex> lock(m);
	exception = true;
	err_str = exc.what();
	continue;
      }

      for(int i=0;i<size;i++){
	iterate<T>::at(s, residuals[i]) = r[i];
	iterate<TC>::at(s, evals[i]) = evals_s[i];
	
	for(int j=0;j<size;j++)
	  iterate<TC>::at(s, evecs[i](j)) = evecs_s[i](j);
      }
    }
    if(exception) throw std::runtime_error(err_str);
    
    return residuals;
  }



};


//Returns residuals
template<typename T>
std::vector<T> symmetricMatrixEigensolve(std::vector<NumericVector<T> > &evecs, std::vector<T> &evals, const NumericSquareMatrix<T> &A, bool sort=true){
  const int size = A.size();
  evecs.resize(size, NumericVector<T>(size));
  evals.resize(size);
  return _SquareMatrixEigensolve<T, hasSampleMethod<T>::value>::solveSymm(evecs,evals,A,sort);
}
template<typename T>
std::vector<T> nonSymmetricMatrixEigensolve(std::vector<NumericVector<Complexify<T> > > &evecs, 
					    std::vector<Complexify<T> > &evals, 
					    const NumericSquareMatrix<T> &A, bool sort=true){
  const int size = A.size();
  evecs.resize(size, NumericVector<Complexify<T> >(size));
  evals.resize(size);
  return _SquareMatrixEigensolve<T, hasSampleMethod<T>::value>::solveNonSymm(evecs,evals,A,sort);
}
template<typename T>
std::vector<T> symmetricMatrixGEVPsolve(std::vector<NumericVector<T> > &evecs, std::vector<T> &evals, const NumericSquareMatrix<T> &A, const NumericSquareMatrix<T> &B, bool sort=true){
  const int size = A.size(); assert(B.size() == A.size());
  evecs.resize(size, NumericVector<T>(size));
  evals.resize(size);
  return _SquareMatrixEigensolve<T, hasSampleMethod<T>::value>::solveSymmGEVP(evecs,evals,A,B,sort);
}

//A, B any square matrices
template<typename T>
std::vector<T> nonSymmetricMatrixGEVPsolve(std::vector<NumericVector<Complexify<T> > > &evecs, std::vector<Complexify<T> > &evals, 
					   const NumericSquareMatrix<T> &A, const NumericSquareMatrix<T> &B, bool sort=true){
  const int size = A.size(); assert(B.size() == A.size());
  evecs.resize(size, NumericVector<Complexify<T> >(size));
  evals.resize(size);
  return _SquareMatrixEigensolve<T, hasSampleMethod<T>::value>::solveNonSymmGEVP(evecs,evals,A,B,sort);
}



CPSFIT_END_NAMESPACE

#endif
