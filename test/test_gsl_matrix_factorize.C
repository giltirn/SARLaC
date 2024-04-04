#include<tensors.h>
#include<random.h>
using namespace SARLaC;

typedef std::complex<double> complexD;

void randomizeSymmetric(NumericSquareMatrix<double> &M, double range_start, double range_end){
  int N = M.size();
  for(int i=0;i<N;i++)
    for(int j=0;j<=i;j++){
      double v = uniformRandom<double>(range_start, range_end, RNG);  
      M(i,j) = v;
      M(j,i) = v;
    }
}
void randomizeComplexHermitian(NumericSquareMatrix<complexD > &M, double range_start, double range_end){
  int N = M.size();
  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      M(i,j) = complexD( uniformRandom<double>(range_start, range_end, RNG), uniformRandom<double>(range_start, range_end, RNG) );
  M = 0.5*M + 0.5*M.dagger();
}
bool isPositiveDefinite(const NumericSquareMatrix<double> &M){
  int N = M.size();
  std::vector<double> evals(N);
  std::vector<NumericVector<double> > evecs(N, NumericVector<double>(N));
  GSLsymmEigenSolver<NumericVector<double>, NumericSquareMatrix<double> >::symmetricMatrixSolve(evecs, evals, M);
  for(int i=0;i<N;i++) if(evals[i]<=0) return false;
  return true;
}
bool isPositiveDefinite(const NumericSquareMatrix<complexD > &M){
  int N = M.size();
  std::vector<double> evals(N);
  std::vector<NumericVector<complexD> > evecs(N, NumericVector<complexD>(N));
  GSLhermEigenSolver<NumericVector<complexD>, NumericSquareMatrix<complexD> >::HermitianMatrixSolve(evecs, evals, M);
  for(int i=0;i<N;i++) if(evals[i]<=0) return false;
  return true;
}

int main(void){
  RNG.initialize(1234);

  int N=4;
  NumericSquareMatrix<double> Mr(N);  
  randomizeSymmetric(Mr,0.1,1.0);
  while(!isPositiveDefinite(Mr))
    randomizeSymmetric(Mr,0.1,1.0);

  NumericSquareMatrix<double> Lr(N);
  GSLmatrixFactorize<NumericSquareMatrix<double> >::Cholesky(Lr,Mr);
  std::cout << Mr << std::endl;
  std::cout << Lr << std::endl;

  NumericSquareMatrix<double> diffr = Lr * Lr.transpose() - Mr;
  double dd = modE(diffr);
  std::cout << "Diff from reconstructed real matrix (expect 0): " << dd << std::endl;
  assert(dd < 1e-5);


  NumericSquareMatrix<complexD> Mz(N);  
  randomizeComplexHermitian(Mz,0.1,1.0);
  while(!isPositiveDefinite(Mz))
    randomizeComplexHermitian(Mz,0.1,1.0);

  NumericSquareMatrix<complexD> Lz(N);
  GSLmatrixFactorize<NumericSquareMatrix<complexD> >::Cholesky(Lz,Mz);
  std::cout << Mz << std::endl;
  std::cout << Lz << std::endl;

  NumericSquareMatrix<complexD> diffz = Lz * Lz.dagger() - Mz;
  dd = modE(diffz).real();
  std::cout << "Diff from reconstructed real matrix (expect 0): " << dd << std::endl;
  assert(dd < 1e-5);

  
  return 0;
}
