#ifndef KTOPIPI_BASIS_CONVERT_H
#define KTOPIPI_BASIS_CONVERT_H

#include<vector>
#include<data_series.h>

SARLAC_START_NAMESPACE

//Both container types should have an operator[](const int q)  to return the data with operator index q  (starting 0)
//Output should be of correct size
template<typename InputContainer, typename OutputContainer>
void convert10to7(OutputContainer &out,
		  const InputContainer &in){
  assert((void*)&out!=(void const*)&in);
  for(int q=0;q<7;q++) out[q] = in[0];

  //Convert Q123 -> Q'123
  static const double Q123rot[3][3] = {  { 3    ,  2,    -1     },
					 { 2./5 , -2./5,  1./5  },
					 {-3./5,   3./5,  1./5  } };

#define Q(i,j) Q123rot[i-1][j-1]      
#define MO(i) out[i-1]
#define MI(i) in[i-1]

  for(int i=1;i<=3;i++)
    MO(i) = Q(i,1)*MI(1) + Q(i,2)*MI(2) + Q(i,3)*MI(3);    
  for(int i=4;i<=7;i++)
    MO(i) = MI(i+1); //5->4  6->5 etc

#undef MO
#undef MI
#undef Q
}


template<typename InputContainer, typename OutputContainer>
void convert7to10(OutputContainer &out,
		  const InputContainer &in){
  assert((void*)&out!=(void const*)&in);
  for(int q=0;q<10;q++) out[q] = in[0];

  //Convert Q'123 -> Q123
  static const double Q123invrot[3][3] = {  {1./5,   1,   0},
					    {1./5,   0,   1},
					    {  0 ,   3,   2} };

#define Qinv(i,j) Q123invrot[i-1][j-1]      
#define MO(i) out[i-1]
#define MI(i) in[i-1]

  for(int i=1;i<=3;i++)
    MO(i) = Qinv(i,1)*MI(1) + Qinv(i,2)*MI(2) + Qinv(i,3)*MI(3);
    
  MO(4) = MO(2) + MO(3) - MO(1); //Q4 = Q2 + Q3 - Q1    [Lehner, Sturm, arXiv:1104.4948 eq 9]
  for(int i=5;i<=8;i++) MO(i) = MI(i-1); //4->5 5->6 etc
    
  MO(9) = 3./2*MO(1)  -1./2*MO(3); //Q9 = 3/2 Q1 - 1/2 Q3
  MO(10) = 1./2*MO(1)  -1./2*MO(3) + MO(2); //Q10 = 1/2(Q1 - Q3) + Q2
    
#undef MO
#undef MI
#undef Q
}

SARLAC_END_NAMESPACE

#endif
