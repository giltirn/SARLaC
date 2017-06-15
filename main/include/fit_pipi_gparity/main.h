#ifndef _PIPI_MAIN_H__
#define _PIPI_MAIN_H__

figureData projectA2(const char fig, const figureDataAllMomenta &raw_data){
  threeMomentum R[8] = { {1,1,1}, {-1,-1,-1},
			 {1,1,-1}, {-1,-1,1},
			 {1,-1,1}, {-1,1,-1},
			 {-1,1,1}, {1,-1,-1} };
  
  figureData out(raw_data.getLt(), raw_data.getNsample()); out.zero();
  for(int psnk=0;psnk<8;psnk++)
    for(int psrc=0;psrc<8;psrc++)
      out = out + raw_data(fig, momComb(R[psnk], R[psrc]));
  out = out / double(64);
  return out;
}

void computeV(figureDataAllMomenta &raw_data, const bubbleDataAllMomenta &raw_bubble_data, const int tsep_pipi){
  threeMomentum R[8] = { {1,1,1}, {-1,-1,-1},
			 {1,1,-1}, {-1,-1,1},
			 {1,-1,1}, {-1,1,-1},
			 {-1,1,1}, {1,-1,-1} };
  const int Lt = raw_bubble_data.getLt();
  
  for(int psnk=0;psnk<8;psnk++)
    for(int psrc=0;psrc<8;psrc++){  
      //B(tsrc + tsep + tsep_pipi, -p1_snk) B(tsrc, p1_src)
      const bubbleData Bmp1_snk = raw_bubble_data( -R[psnk] );
      const bubbleData Bp1_src  = raw_bubble_data(  R[psrc] );

      figureData &into = raw_data('V',momComb(R[psnk],R[psrc]));
      
      for(int tsrc=0;tsrc<Lt;tsrc++)
	for(int tsep=0;tsep<Lt;tsep++)
	  into(tsrc,tsep) = Bmp1_snk( (tsrc + tsep + tsep_pipi) % Lt ) * Bp1_src( tsrc );
    }
}

bubbleDataDoubleJackAllMomenta doubleJackknifeResampleBubble(const bubbleDataAllMomenta &bubbles){
  int Lt = bubbles.getLt();
  bubbleDataDoubleJackAllMomenta out(Lt,bubbles.getNsample());
  for(auto it = bubbles.begin(); it != bubbles.end(); it++){
    const threeMomentum & mom = it->first;
    const typename bubbleDataAllMomenta::ContainerType & raw = it->second;
    for(int t=0;t<Lt;t++) out(mom)(t).resample(raw(t));
  }
  return out;
}


#endif
