#ifndef _PIPI_RESAMPLED_DATA_H_
#define _PIPI_RESAMPLED_DATA_H_

#include<tuple>

#include<config.h>
#include<utils/macros.h>

#include "mom_data_containers.h"
#include "mom_project.h"

CPSFIT_START_NAMESPACE

//Combine the computation of the V diagram with A2 projection and source average to avoid large intermediate data storage
template<typename BubbleDataType>
auto computePiPi2ptFigureVprojectSourceAvg(const BubbleDataType &raw_bubble_data, const int tsep_pipi, const PiPiCorrelatorSelector &corr_select, const std::vector<threeMomentum> &pion_momenta)
  ->correlationFunction<double,typename std::decay<decltype(raw_bubble_data(Source,*((threeMomentum*)NULL))(0))>::type>{

  (std::cout << "Computing projected, src-averaged V diagrams with BubbleDataType = " << printType<BubbleDataType>() << " and " << omp_get_max_threads() << " threads\n").flush(); 
  boost::timer::auto_cpu_timer t(std::string("Report: Computed projected, src-averaged V diagrams with BubbleType = ") + printType<BubbleDataType>() + " in %w s\n");

  const int nmom = pion_momenta.size();
  
  const int Lt = raw_bubble_data.getLt();
  const int Nsample = raw_bubble_data.getNsample();

  typedef typename std::decay<decltype(raw_bubble_data(Source,*((threeMomentum*)NULL))(0))>::type  DistributionType;
  correlationFunction<double,DistributionType> out(Lt,
					    [&](const int t)
					    {
					      return typename correlationFunction<double,DistributionType>::ElementType(double(t), DistributionType(Nsample,0.));
					    }
					    );
  int nthr = omp_get_max_threads();
  std::vector<correlationFunction<double,DistributionType> > thr_sum(nthr, out);
  
  std::vector<std::tuple<int,int,double> > todo;
  double coeff;
  for(int psnk=0;psnk<nmom;psnk++)
    for(int psrc=0;psrc<nmom;psrc++)
      if(corr_select(coeff,pion_momenta[psrc],pion_momenta[psnk])){
	todo.push_back(std::make_tuple(psrc,psnk,coeff));
      }

#pragma omp parallel for
  for(int pp=0;pp<todo.size();pp++){
    int me = omp_get_thread_num();
    int psrc = std::get<0>(todo[pp]);
    int psnk = std::get<1>(todo[pp]);
    double coeff = std::get<2>(todo[pp]);

    const auto &Bp1_snk = raw_bubble_data(Sink, pion_momenta[psnk] ); 
    const auto &Bp1_src  = raw_bubble_data(Source,  pion_momenta[psrc] );

    for(int tsrc=0;tsrc<Lt;tsrc++)
      for(int tsep=0;tsep<Lt;tsep++)
	thr_sum[me].value(tsep) = thr_sum[me].value(tsep) + coeff * Bp1_snk( (tsrc + tsep) % Lt ) * Bp1_src( tsrc );
  }
  for(int tsep=0;tsep<Lt;tsep++){
    for(int thr=0;thr<nthr;thr++)
      out.value(tsep) = out.value(tsep) + thr_sum[thr].value(tsep);
    out.value(tsep) = out.value(tsep)/double(Lt);
  }
  return out;  
}

bubbleDataDoubleJackAllMomenta binDoubleJackknifeResampleBubble(const bubbleDataAllMomenta &bubbles, const int bin_size){
  int Lt = bubbles.getLt();
  bubbleDataDoubleJackAllMomenta out(Lt, bubbles.getTsepPiPi(), bubbles.getNsample()/bin_size);

  for(auto it = bubbles.begin(); it != bubbles.end(); it++){
    const auto & key = it->first;
    const typename bubbleDataAllMomenta::ContainerType & raw = it->second;
    for(int t=0;t<Lt;t++)
      out(key)(t).resample(raw(t).bin(bin_size));
  }
  return out;
}

/*
Generic fold routine. 
The data typically is symmetric as   C(Lt - fold_offset - t) ~ C(t) 
fold_offset should be  2*tsep_pipi  for pipi2pt,   tsep_pipi  for pipi->sigma and 0 for sigma 2pt
*/
template<typename T>
inline correlationFunction<double,T> fold(const correlationFunction<double,T> &f, const int fold_offset){
  const int Lt = f.size();
  correlationFunction<double,T> out(Lt);
  const int Tref = Lt-fold_offset;
  for(int t=0;t<Lt;t++){
    out.coord(t) = f.coord(t);
    out.value(t) = ( f.value(t) + f.value( (Tref-t+Lt) % Lt ) )/2.;
  }
  return out;
}


template<typename T>
inline correlationFunction<double,T> foldPiPi2pt(const correlationFunction<double,T> &f, const int tsep_pipi){
  return fold(f, 2*tsep_pipi);
}




CPSFIT_END_NAMESPACE

#endif
