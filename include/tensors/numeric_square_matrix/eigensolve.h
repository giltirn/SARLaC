#ifndef _CPSFIT_NUMERIC_SQUARE_MATRIX_EIGENSOLVE_H_
#define _CPSFIT_NUMERIC_SQUARE_MATRIX_EIGENSOLVE_H_

#include<distribution/distribution_iterate.h>
#include<tensors/gsl_eigensolve.h>
#include<tensors/numeric_vector.h>
#include<tensors/numeric_square_matrix/class.h>

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
    return GSLsymmEigenSolver<NumericVector<T>, NumericSquareMatrix<T> >::symmetricGEVPsolve(evecs,evals,A,B,sort);
  }
};

template<typename T, int is_distribution>
struct _complexify{
};
template<typename T>
struct _complexify<T,0>{
  typedef std::complex<T> type;
  typedef NumericVector<type> VectorType;
};
template<typename T>
struct _complexify<T,1>{
  typedef typename T::template rebase<std::complex<typename T::DataType> > type;
  typedef NumericVector<type> VectorType;
};

template<typename T>
struct complexify{
  typedef typename _complexify<T, hasSampleMethod<T>::value>::type type;
  typedef typename _complexify<T, hasSampleMethod<T>::value>::VectorType VectorType;
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

  static std::vector<T> solveNonSymm(std::vector<typename complexify<T>::VectorType> &evecs, 
				     std::vector<typename complexify<T>::type> &evals, 
				     const NumericSquareMatrix<T> &A, bool sort = true){
    assert(A.size() > 0);
    const int size = A.size();
    const int nsample = A(0,0).size();

    for(int i=0;i<size;i++){
      evals[i].resize(nsample);
      for(int j=0;j<size;j++)
	evecs[i][j].resize(nsample);
    }
    
    typedef typename complexify<T>::VectorType VectorDistributionType;
    typedef typename complexify<T>::type ComplexDistributionType;
    
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
	iterate<ComplexDistributionType>::at(s, evals[i]) = evals_s[i];
	
	for(int j=0;j<size;j++)
	  iterate<ComplexDistributionType>::at(s, evecs[i](j)) = evecs_s[i](j);
      }
    }
    return residuals;
  }




  //GEVP
  static std::vector<T> solveSymmGEVP(std::vector<NumericVector<T> > &evecs, std::vector<T> &evals, const NumericSquareMatrix<T> &A, const NumericSquareMatrix<T> &B, bool sort = true){
    assert(A.size() > 0);
    const int size = A.size();
    const int nsample = A(0,0).size();
    assert(B.size() == size && B(0,0).size() == nsample);

    for(int i=0;i<size;i++){
      evals[i].resize(nsample);
      for(int j=0;j<size;j++)
	evecs[i][j].resize(nsample);
    }
    
    typedef typename iterate<T>::type type;
    NumericSquareMatrix<type> A_s(size);
    NumericSquareMatrix<type> B_s(size);
    std::vector<NumericVector<type> > evecs_s(size, NumericVector<type>(size));
    std::vector<type> evals_s(size);

    std::vector<T> residuals(size, T(nsample));
    
    const int nit = iterate<T>::size(A(0,0));
    for(int s=0;s<nit;s++){
      for(int i=0;i<size;i++)
	for(int j=0;j<size;j++){
	  A_s(i,j) = iterate<T>::at(s, A(i,j));
	  B_s(i,j) = iterate<T>::at(s, B(i,j));
	}

      std::vector<double> r = GSLsymmEigenSolver<NumericVector<type>, NumericSquareMatrix<type> >::symmetricGEVPsolve(evecs_s,evals_s,A_s,B_s, sort);

      for(int i=0;i<size;i++){
	iterate<T>::at(s, residuals[i]) = r[i];
	iterate<T>::at(s, evals[i]) = evals_s[i];
	
	for(int j=0;j<size;j++)
	  iterate<T>::at(s, evecs[i](j)) = evecs_s[i](j);
      }
    }
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
std::vector<T> nonSymmetricMatrixEigensolve(std::vector<typename complexify<T>::VectorType> &evecs, 
					    std::vector<typename complexify<T>::type> &evals, 
					    const NumericSquareMatrix<T> &A, bool sort=true){
  const int size = A.size();
  evecs.resize(size, typename complexify<T>::VectorType(size));
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



CPSFIT_END_NAMESPACE

#endif
