//Check the wrapper around the general GEVP solve by comparing to one in which B is positive definite
#include<tensors.h>
#include<random.h>
using namespace SARLaC;

void randomizeSymmetric(NumericSquareMatrix<double> &M, double range_start, double range_end){
  int N = M.size();
  for(int i=0;i<N;i++)
    for(int j=0;j<=i;j++){
      double v = uniformRandom<double>(range_start, range_end, RNG);  
      M(i,j) = v;
      M(j,i) = v;
    }
}



bool isPositiveDefinite(const NumericSquareMatrix<double> &M){
  int N = M.size();
  std::vector<double> evals(N);
  std::vector<NumericVector<double> > evecs(N, NumericVector<double>(N));
  GSLsymmEigenSolver<NumericVector<double>, NumericSquareMatrix<double> >::symmetricMatrixSolve(evecs, evals, M);

  // std::cout << "Evals: ";
  // for(int i=0;i<N;i++) std::cout << evals[i] << " ";
  // std::cout << std::endl;

  for(int i=0;i<N;i++) if(evals[i]<=0) return false;
  return true;
}

int main(void){
  RNG.initialize(1234);

  int N=4;
  NumericSquareMatrix<double> A(N),B(N);

  randomizeSymmetric(A, 0.1,1);
  std::cout << "A: " << A.print() << std::endl;

  int attempts = 1;
  while(1){
    randomizeSymmetric(B, 0.1,1);
    if(isPositiveDefinite(B)){
      std::cout << "B: " << B.print() << std::endl << "obtained after " << attempts << " attempts" << std::endl;
      break;
    }
    attempts++;
  }

  std::vector<double> evals_1(N);
  std::vector<NumericVector<double> > evecs_1(N, NumericVector<double>(N));
  GSLsymmEigenSolver<NumericVector<double>, NumericSquareMatrix<double> >::symmetricGEVPsolve(evecs_1, evals_1, A,B);

  std::vector<std::complex<double> > evals_2(N);
  std::vector<NumericVector<std::complex<double> > > evecs_2(N, NumericVector<std::complex<double> >(N));
  GSLnonSymmEigenSolver<NumericVector<std::complex<double> >, NumericSquareMatrix<double> >::nonSymmetricGEVPsolve(evecs_2, evals_2, A,B);

  std::vector<double> evals_3(N);
  std::vector<NumericVector<double> > evecs_3(N, NumericVector<double>(N));
  GSLsymmEigenSolver<NumericVector<double>, NumericSquareMatrix<double> >::symmetricGEVPsolveGen(evecs_3, evals_3, A,B);
  

  //Check evals are real for second
  for(int i=0;i<N;i++)
    assert( imag(evals_2[i])/abs(evals_2[i]) < 1e-5);
  
  //The second sorts by the absolute complex value whereas the first sorts by the value itself. Thus if the matrix has negative eigenvalues it will sort differently
  //Compute a map to fix the ordering
  std::vector<int> map_2(N);
  for(int i=0;i<N;i++) map_2[i] = i;
  std::sort(map_2.begin(),map_2.end(), [&](const int i, const int j){ return real(evals_2[i])> real(evals_2[j]); });

  std::cout << "Compare 1,2" << std::endl;
  for(int i=0;i<N;i++){
    double a = evals_1[i];
    double b = real(evals_2[map_2[i]]);
    double diff = 2*(b-a)/(b+a);

    std::cout << a << " " << b << " " << diff  << std::endl;

    assert(diff < 1e-5);
  }

  std::cout << "Compare 1,3" << std::endl;
  for(int i=0;i<N;i++){
    double a = evals_1[i];
    double b = evals_3[i];
    double diff = 2*(b-a)/(b+a);

    std::cout << a << " " << b << " " << diff  << std::endl;

    assert(diff < 1e-5);
  }
  


  return 0;
}
