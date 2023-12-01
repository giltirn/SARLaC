#ifndef _FIT_PIPI_COMOVING_GEVP_RESAMPLED_DATA_H
#define _FIT_PIPI_COMOVING_GEVP_RESAMPLED_DATA_H

void saveCheckpoint(const std::map<threeMomentum, ResampledData<jackknifeCorrelationFunctionD> > &data_j,
		    const std::string &file){
  std::cout << "Saving data checkpoint to " << file << std::endl;
  HDF5writer wr(file);
  write(wr, data_j, "j_data");
}


void loadCheckpoint(std::map<threeMomentum, ResampledData<jackknifeCorrelationFunctionD> > &data_j,
		    const std::string &file){
  std::cout << "Reading data checkpoint from " << file << std::endl;
  HDF5reader rd(file);
  read(rd, data_j, "j_data");
}

template<typename DistributionType>
correlationFunction<double, NumericSquareMatrix<DistributionType> > createCorrelatorMatrix(const std::map<threeMomentum, ResampledData<correlationFunction<double,DistributionType> > > &data,
											   const int Lt, const std::vector<Operator> &ops, const std::vector<std::array<int,3> > &p_tot){
  const int nop = ops.size();
  correlationFunction<double, NumericSquareMatrix<DistributionType> > C(Lt);
  for(int t=0;t<Lt;t++){
    C.coord(t) = t;
    C.value(t).resize(nop);

    for(int i=0;i<nop;i++){
      for(int j=i;j<nop;j++){
	//Average over the total momentum values (it is assumed these are equivalent)
	auto const &data_p = data.find(p_tot[0])->second;
	DistributionType value = data_p.correlator(ops[i],ops[j]).value(t);

	for(int p=1;p<p_tot.size();p++){
	  auto const &data_p = data.find(p_tot[p])->second;
	  value = value + data_p.correlator(ops[i],ops[j]).value(t);
	}
	value = value / double(p_tot.size());
	
	C.value(t)(i,j) = C.value(t)(j,i) = value;
      }
    }
  }
  return C;
}

#endif
