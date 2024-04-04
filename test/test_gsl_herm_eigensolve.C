#include<tensors.h>
#include<random.h>
using namespace SARLaC;


void randomizeComplexHermitian(NumericSquareMatrix<std::complex<double> > &M, double range_start, double range_end){
  int N = M.size();
  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      M(i,j) = std::complex<double>( uniformRandom<double>(range_start, range_end, RNG), uniformRandom<double>(range_start, range_end, RNG) );
  M = 0.5*M + 0.5*M.dagger();
}
int main(void){
  RNG.initialize(1234);

  typedef std::complex<double> complexD;

  int N=4;
  NumericSquareMatrix<complexD> A(N);
  randomizeComplexHermitian(A, 0.1,1);
  std::cout << "A: " << A.print() << std::endl;
  
  //Test matrix is Hermitian
  NumericSquareMatrix<complexD> diff = A - A.dagger();
  double dd = modE(diff).real();
  std::cout << "Diff from Hermitian matrix (expect 0): " << dd << std::endl;
  assert(dd < 1e-5);

  std::vector<double> evals(N);
  std::vector<NumericVector<complexD> > evecs(N, NumericVector<complexD>(N));
  std::vector<double> resid = GSLhermEigenSolver<NumericVector<complexD>, NumericSquareMatrix<complexD> >::HermitianMatrixSolve(evecs, evals, A, true);
  std::cout << "Residuals:";
  for(int i=0;i<N;i++) std::cout << " " << resid[i];
  std::cout << std::endl;

  //Check are evecs
  for(int i=0;i<N;i++){
    NumericVector<complexD> Mv_m_lv = A * evecs[i] - evals[i] * evecs[i];
    std::cout << "Resid vec " << i << ": " << Mv_m_lv << std::endl;
  }


  //Check evals in descending order
  for(int i=1;i<N;i++) assert( evals[i] <= evals[i-1] );

  //Check eigendecomposition
  NumericSquareMatrix<complexD> Ar(N);
  Ar.zero();

  for(int l=0;l<N;l++)
    for(int i=0;i<N;i++)
      for(int j=0;j<N;j++)
	Ar(i,j) = Ar(i,j) + evecs[l](i) * evals[l] * std::conj(evecs[l](j));

  diff = Ar - A;
  
  dd = modE(diff).real();
  std::cout << "Diff from reconstructed matrix (expect 0): " << dd << std::endl;
  assert(dd < 1e-5);
  return 0;
}
  
