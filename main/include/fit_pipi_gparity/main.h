#ifndef _PIPI_MAIN_H__
#define _PIPI_MAIN_H__

template<typename DataAllMomentumType>
typename DataAllMomentumType::ContainerType projectA2(const char fig, const DataAllMomentumType &raw_data){
  std::cout << "Computing A2 projection of figure " << fig << " with DataAllMomentumType = " << printType<DataAllMomentumType>() << "\n"; 
  boost::timer::auto_cpu_timer t(std::string("Report: Computed A2 projection of figure ") + fig + " with DataAllMomentumType = " + printType<DataAllMomentumType>() + " in %w s\n");
  threeMomentum R[8] = { {1,1,1}, {-1,-1,-1},
			 {1,1,-1}, {-1,-1,1},
			 {1,-1,1}, {-1,1,-1},
			 {-1,1,1}, {1,-1,-1} };
  
  typename DataAllMomentumType::ContainerType out(raw_data.getLt(), raw_data.getNsample()); out.zero();

  for(int psnk=0;psnk<8;psnk++)
    for(int psrc=0;psrc<8;psrc++)
      out = out + raw_data(fig, momComb(R[psnk], R[psrc]));

  out = out / double(64);
  return out;
}

template<typename DataAllMomentumType, typename BubbleDataType>
void computeV(DataAllMomentumType &raw_data, const BubbleDataType &raw_bubble_data, const int tsep_pipi){
  std::cout << "Computing V diagrams with BubbleDataType = " << printType<BubbleDataType>() << "\n"; 
  boost::timer::auto_cpu_timer t(std::string("Report: Computed V diagrams with BubbleType = ") + printType<BubbleDataType>() + " in %w s\n");
				 
  threeMomentum R[8] = { {1,1,1}, {-1,-1,-1},
			 {1,1,-1}, {-1,-1,1},
			 {1,-1,1}, {-1,1,-1},
			 {-1,1,1}, {1,-1,-1} };
  const int Lt = raw_bubble_data.getLt();
  const int Nsample = raw_bubble_data.getNsample();
  raw_data.setup(Lt,Nsample);

  for(int psnk=0;psnk<8;psnk++)
    for(int psrc=0;psrc<8;psrc++){  
      //B(tsrc + tsep + tsep_pipi, -p1_snk) B(tsrc, p1_src)
      const auto &Bmp1_snk = raw_bubble_data( -R[psnk] );
      const auto &Bp1_src  = raw_bubble_data(  R[psrc] );

      auto &into = raw_data('V',momComb(R[psnk],R[psrc]));
      
      for(int tsrc=0;tsrc<Lt;tsrc++)
	for(int tsep=0;tsep<Lt;tsep++)
	  into(tsrc,tsep) = Bmp1_snk( (tsrc + tsep + tsep_pipi) % Lt ) * Bp1_src( tsrc );
    }
}



  
template<typename FigureDataType>
auto realSourceAverage(const FigureDataType & data)->correlationFunction<typename std::decay<decltype(data(0,0).real())>::type>{
  typedef typename std::decay<decltype(data(0,0).real())>::type RealDistributionType;


  int Lt = data.getLt();
  int nsample = data.getNsample();
  correlationFunction<RealDistributionType> into(data.getLt(),
						  [nsample](int i) {  return typename correlationFunction<RealDistributionType>::ElementType(i, RealDistributionType(nsample,0.)); }
						  );
  std::vector<int> tsrc_include;
  for(int tsrc=0;tsrc<Lt;tsrc++){
    bool is_nonzero = !data.isZero(tsrc);
    if(is_nonzero)
      tsrc_include.push_back(tsrc);
  }
  double N(tsrc_include.size());

  for(int tsep=0;tsep<Lt;tsep++){
    auto & v = into.value(tsep);
    v.zero();
    for(int i=0;i<tsrc_include.size();i++){
      //auto re = data(tsrc_include[i],tsep).real();
      //v = v + re;

      v = v + data(tsrc_include[i],tsep).real();      
    }
    v = v/N;
  }
  return into;
}


bubbleDataDoubleJackAllMomenta doubleJackknifeResampleBubble(const bubbleDataAllMomenta &bubbles){
  int Lt = bubbles.getLt();
  bubbleDataDoubleJackAllMomenta out(Lt,bubbles.getNsample());

  for(auto it = bubbles.begin(); it != bubbles.end(); it++){
    const threeMomentum & mom = it->first;
    const typename bubbleDataAllMomenta::ContainerType & raw = it->second;
    for(int t=0;t<Lt;t++)
      out(mom)(t).resample(raw(t));
  }
  return out;
}

template<typename T>
inline correlationFunction<T> fold(const correlationFunction<T> &f, const int tsep_pipi){
  const int Lt = f.size();
  correlationFunction<T> out(Lt);
  const int Tref = Lt-2*tsep_pipi;
  for(int t=0;t<Lt;t++){
    out.coord(t) = f.coord(t);
    out.value(t) = ( f.value(t) + f.value( (Tref-t+Lt) % Lt ) )/2.;
  }
  return out;
}


#endif
