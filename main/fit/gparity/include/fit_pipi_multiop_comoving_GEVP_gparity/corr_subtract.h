#ifndef _FIT_PIPI_COMOVING_GEVP_GPARITY_CORR_SUBTRACT_H
#define _FIT_PIPI_COMOVING_GEVP_GPARITY_CORR_SUBTRACT_H

#define SUBDATA_MEMBERS \
  (int, state0)					\
  (int, state1)					\
  (std::string, file)				\
  (std::vector<int>, idx)

struct SubData{
  GENERATE_MEMBERS(SUBDATA_MEMBERS);
  SubData(): state0(0), state1(1), file("myfit.hdf5"), idx({4}){}
};
GENERATE_PARSER(SubData, SUBDATA_MEMBERS);

#define SUBARGS_MEMBERS \
  (std::vector<SubData>, sub) \
  (double, scale)

struct SubArgs{
  GENERATE_MEMBERS(SUBARGS_MEMBERS);
  SubArgs(): sub(1), scale(1e13){}
};
GENERATE_PARSER(SubArgs, SUBARGS_MEMBERS);

void correlatorSubtract(correlationFunction<double, NumericSquareMatrix<jackknifeDistributionD> > &C,
			std::string subtract_from_data_file){
  const int nsample = C.value(0)(0,0).size();
  const int nop = C.value(0).size();
  SubArgs subargs;
  parse(subargs, subtract_from_data_file);
    
  NumericSquareMatrix<jackknifeDistributionD> M(nop, jackknifeDistributionD(nsample,0.));
    
  for(int i=0;i<subargs.sub.size();i++){
    jackknifeDistributionD c;
    readHDF5file(c, subargs.sub[i].file, subargs.sub[i].idx);
    M(subargs.sub[i].state0, subargs.sub[i].state1) = M(subargs.sub[i].state1, subargs.sub[i].state0) = c * subargs.scale;
  }
  std::cout << "Subtracting constant matrix:\n" << M << std::endl;
    
  for(int t=0;t<C.size();t++){
    //std::cout << t << "\n" << C << "\n ----> \n";
    C.value(t) = C.value(t) - M;
    //std::cout << C << std::endl;
  }
}

void correlatorSubtract(correlationFunction<double, NumericSquareMatrix<doubleJackknifeDistributionD> > &C,
			std::string subtract_from_data_file){
  const int nsample = C.value(0)(0,0).size();
  const int nop = C.value(0).size();
  SubArgs subargs;
  parse(subargs, subtract_from_data_file);
    
  NumericSquareMatrix<jackknifeDistributionD> M(nop, jackknifeDistributionD(nsample,0.));
    
  for(int i=0;i<subargs.sub.size();i++){
    jackknifeDistributionD c;
    readHDF5file(c, subargs.sub[i].file, subargs.sub[i].idx);
    M(subargs.sub[i].state0, subargs.sub[i].state1) = M(subargs.sub[i].state1, subargs.sub[i].state0) = c * subargs.scale;
  }
  std::cout << "Subtracting constant matrix:\n" << M << std::endl;
    
  for(int t=0;t<C.size();t++){
    for(int a=0;a<nop;a++)
      for(int b=0;b<nop;b++)
	for(int o=0;o<nsample;o++)
	  for(int i=0;i<nsample-1;i++)
	    C.value(t)(a,b).sample(o).sample(i) = C.value(t)(a,b).sample(o).sample(i) - M(a,b).sample(o); //same value subtracted for each inner jackknife sample	
  }
}

//C(t) -> C(t) - C(t+1)
template<typename DistributionType>
void correlatorSubtractNeighbor(correlationFunction<double, NumericSquareMatrix<DistributionType> > &C){
  correlationFunction<double, NumericSquareMatrix<DistributionType> > out(C.size()-1);
  for(int t=0;t<C.size()-1;t++){
    out.coord(t) = t;
    out.value(t) = C.value(t) - C.value(t+1);
  }
  C = out;
}

template<typename DistributionType>
void correlatorSubtractFixedT(correlationFunction<double, NumericSquareMatrix<DistributionType> > &C, const int tsub){
  int nop = C.value(0).size();
  for(int t=0;t<C.size();t++)
    for(int i=0;i<nop;i++)
      for(int j=i;j<nop;j++)
	C.value(t)(i,j) = C.value(t)(j,i) = C.value(t)(i,j) - C.value(tsub)(i,j);
}


#endif
