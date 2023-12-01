#include<array>
#include<vector>
#include<complex>
#include<iostream>
#include<fstream>

#include<tensors.h>

using namespace SARLaC;

int main(void){
  {
    NumericSquareMatrix<double> m(2, [&](const int i,const int j){ return i == j ? 2. : 1.; });
    NumericVector<double> v(2);  v(0) = 1; v(1) = 2;

    NumericVector<double> mv = m*v;

    std::cout << "m*v : \n" << m << "\n*\n" << v << "\n=\n" << mv << std::endl;

    std::cout << "v.v = " << dot(v,v) << std::endl;
    std::cout << "mv.mv = " << dot(mv,mv) << std::endl;
  }
  
  typedef std::initializer_list<int> IL;

  {
    std::cout << "Tensor rank 2, size (4,4):\n";
    size_t vol = _tensor_helper<double,2,1>::vol(IL({4,4}).begin() );
    std::cout << "Volume: " << vol << std::endl;
    std::cout << "Mapping:\n";
    int unmap[2];
    for(int i=0;i<4;i++)
      for(int j=0;j<4;j++){
	size_t off = _tensor_helper<double,2,1>::map(0,IL({i,j}).begin() ,IL({4,4}).begin() );      
	std::cout << "(" << i << "," << j << ")->" << off;	
	_tensor_helper<double,2,1>::unmap(unmap,off,IL({4,4}).begin(),vol);
	std::cout << "->[" << unmap[0] << "," << unmap[1] << "]\n";
      }
      
  }

  {
    NumericTensor<double,2> t({2,2});
    t({0,0}) = 0;
    t({0,1}) = 1;
    t({1,0}) = 2;
    t({1,1}) = 3;

    std::cout << "Test 2x2 matrix:\n" << t << std::endl;

    NumericTensor<double,2> u({2,2});
    u({0,0}) = 4;
    u({0,1}) = 5;
    u({1,0}) = 6;
    u({1,1}) = 7;

    std::cout << "Second 2x2 matrix:\n" << u << std::endl;
  
    NumericTensor<double,2> v = t + u;
    std::cout << "Test 2x2 matrix sum:\n" << v << std::endl;
  }
  {
    std::cout << "Test transform\n";
    NumericTensor<double,2> w({2,2});
    w({0,0}) = 1.123;
    w({0,1}) = 6.79;
    w({1,0}) = -5.777;
    w({1,1}) = 9.009;
    
    NumericTensor<unsigned int,2> wt = w.transform( [&](const int *coord, const double &from){ return (unsigned int)(fabs(from)); } );
    
    std::cout << w << "\n->\n" << wt << std::endl;

    std::cout << "Test reduce by summing over rows\n";
    NumericTensor<double,1> wr = w.reduce(1, [](double &into, int const* coord, const NumericTensor<double,2> &m){ into = m({coord[0],0}) + m({coord[0],1}); } );
    
    std::cout << w << "\n->\n" << wr << std::endl;
  }
  {
    std::cout << "Test rank-2 tensor multiply:\n";

    NumericTensor<double,2> t({2,3});
    t({0,0}) = 0;
    t({0,1}) = 1;
    t({0,2}) = 2;
    t({1,0}) = 3;
    t({1,1}) = 4;
    t({1,2}) = 5;

    NumericTensor<double,2> u({3,2});
    u({0,0}) = 0;
    u({0,1}) = 1;
    u({1,0}) = 2;
    u({1,1}) = 3;
    u({2,0}) = 4;
    u({2,1}) = 5;

    NumericTensor<double,2> v = t * u;

    std::cout << t << "\n*\n" << u << "\n=\n" << v << std::endl;

    std::cout << "Test rank-2 tensor contract:\n";    
    NumericTensor<double,2> vv = contract(t,u,1,0);

    std::cout << t << "\n*\n" << u << "\n=\n" << v << std::endl;

    std::cout << "Test matrix-vector multiplication:\n";
    NumericTensor<double,1> vec({3});
    vec({0}) = 1;
    vec({1}) = 2;
    vec({2}) = 3;

    NumericTensor<double,1> ovec = t*vec;
    std::cout << t << "\n*\n" << vec << "\n=\n" << ovec << std::endl;

    std::cout << "Test matrix-vector multiplication via contraction:\n";
    NumericTensor<double,1> ovec2 = contract(t,vec,1,0);
    std::cout << t << "\n*\n" << vec << "\n=\n" << ovec2 << std::endl;
  }

  //Test 2x2 determinant
  {
    NumericSquareMatrix<double> M({1.3,-4.5, 0.3, 72.4});
    double expect = M(0,0)* M(1,1) - M(0,1)*M(1,0);
    double got = determinant(M);
    std::cout << "2x2 determinant got " << got << " expect " << expect << std::endl;
    assert(fabs(expect-got) < 1e-10);
  }
  //Test 3x3 determinant
  {
    NumericSquareMatrix<double> M(
				  {1.3,-4.5, 0.3, 
				   72.4, 31.2, -0.07,
				   12.1, 4e-2, 5e2 }
				  );
    double expect = M(0,0) * ( M(1,1) * M(2,2) - M(2,1) * M(1,2) )
      - M(0,1) * ( M(1,0)*M(2,2) - M(2,0)*M(1,2) )
      + M(0,2) * ( M(1,0)*M(2,1) - M(2,0)*M(1,1) );
    
    double got = determinant(M);
    std::cout << "3x3 determinant got " << got << " expect " << expect << std::endl;
    assert(fabs(expect-got) < 1e-10);
  }



  
  return 0;
}
