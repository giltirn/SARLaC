#include<physics/luscherzeta.h>

using namespace CPSfit;

//This is just a basic calculator with assignable symbols!
int main(){
  //Test by comparing with numbers computed by Daiqian (based off Qi's MatLab code) 
  std::array<double,3> d = {0,0,0}; double gamma = 1.;
    
  {
    //q  dx dy dz phi
    //0.100000	0	0	1	0.790456
    LuscherZeta zeta({0,0,1},d);
    printf("0.100000	0	0	1	Expect 0.790456 Got %.6f\n", zeta.calcPhi(0.1, gamma));
  }      
  {
    //0.300000	1	1	1	1.079671
    LuscherZeta zeta({1,1,1},d);
    printf("0.300000	1	1	1	Expect 1.079671 Got %.6f\n", zeta.calcPhi(0.3, gamma));
  }

  {
    //0.400000	0	0	0	0.571791
    LuscherZeta zeta({0,0,0},d);
    printf("0.400000	0	0	0	Expect 0.571791  Got %.6f\n", zeta.calcPhi(0.4, gamma));
  }

  {
    //2.000000	1	1	1	-1.315088
    LuscherZeta zeta({1,1,1},d);
    printf("2.000000	1	1	1	Expect -1.315088 Got %.6f\n", zeta.calcPhi(2.0, gamma));
  }
  return 0;
}
