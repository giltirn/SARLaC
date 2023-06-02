#ifndef _PIPI_RESAMPLED_DATA_H_
#define _PIPI_RESAMPLED_DATA_H_

#include<tuple>

#include<config.h>
#include<utils/macros.h>

#include "../mom_data_containers.h"
#include "../correlator_utils/mom_project.h"

CPSFIT_START_NAMESPACE

//Combine the computation of the V diagram with A2 projection and source average to avoid large intermediate data storage
template<typename BubbleDataType>
auto computePiPi2ptFigureVprojectSourceAvg(const BubbleDataType &raw_bubble_data, const int tsep_pipi, const PiPiProjectorBase &proj_src, const PiPiProjectorBase &proj_snk)
  ->correlationFunction<double,typename BubbleDataType::DistributionType>{

  (std::cout << "Computing projected, src-averaged V diagrams with BubbleDataType = " << printType<BubbleDataType>() << " and " << omp_get_max_threads() << " threads\n").flush(); 
  boost::timer::auto_cpu_timer t(std::string("Report: Computed projected, src-averaged V diagrams with BubbleType = ") + printType<BubbleDataType>() + " in %w s\n");

  const int Lt = raw_bubble_data.getLt();

  typedef typename BubbleDataType::DistributionType DistributionType;

  DistributionType zero(raw_bubble_data.begin()->second(0)); zeroit(zero);

  correlationFunction<double,DistributionType> out(Lt,
					    [&](const int t)
					    {
					      return typename correlationFunction<double,DistributionType>::ElementType(double(t), zero);
					    }
					    );
  int nthr = omp_get_max_threads();
  std::vector<correlationFunction<double,DistributionType> > thr_sum(nthr, out);
  std::vector<DistributionType> thr_tmp(nthr,zero);

  int npsrc = proj_src.nMomenta();
  int npsnk = proj_snk.nMomenta();

#pragma omp parallel for
  for(int pp=0;pp<npsrc*npsnk;pp++){ //psnk + npsnk*psrc
    int me = omp_get_thread_num();
    threeMomentum psrc = proj_src.momentum(pp/npsnk);
    threeMomentum psnk = proj_snk.momentum(pp%npsnk);
    double coeff = std::real(proj_src.coefficient(pp/npsnk) * proj_snk.coefficient(pp%npsnk));

    const auto &Bp1_snk = raw_bubble_data(Sink, psnk); 
    const auto &Bp1_src  = raw_bubble_data(Source,  psrc);

    for(int tsrc=0;tsrc<Lt;tsrc++){
      thr_tmp[me] = coeff * Bp1_src( tsrc );
      for(int tsep=0;tsep<Lt;tsep++)
	thr_sum[me].value(tsep) = thr_sum[me].value(tsep) + Bp1_snk( (tsrc + tsep) % Lt ) * thr_tmp[me];
    }
  }
#pragma omp parallel for
  for(int tsep=0;tsep<Lt;tsep++){
    for(int thr=0;thr<nthr;thr++)
      out.value(tsep) = out.value(tsep) + thr_sum[thr].value(tsep);
    out.value(tsep) = out.value(tsep)/double(Lt);
  }
  return out;  
}

template<typename resampledCorrelationFunction, typename binResampler>
resampledCorrelationFunction binResample(const rawDataCorrelationFunctionD &raw, const binResampler &resampler){
  typedef typename resampledCorrelationFunction::DataType DistributionType;
  return resampledCorrelationFunction(raw.size(),
				       [&](const int t){ 
					return typename resampledCorrelationFunction::ElementType(t, binResample<DistributionType>(raw.value(t),resampler));
				      }
				      );
}


template<typename resampledBubbleDataAllMomentumType, typename binResampler>
resampledBubbleDataAllMomentumType binResample(const bubbleDataAllMomenta &bubbles, const binResampler &resampler){
  int Lt = bubbles.getLt();

  resampledBubbleDataAllMomentumType out(Lt, bubbles.getTsepPiPi());
  typedef typename resampledBubbleDataAllMomentumType::DistributionType DistributionType;
  
  for(auto it = bubbles.begin(); it != bubbles.end(); it++){
    for(int t=0;t<Lt;t++)
      out(it->first)(t) = binResample<DistributionType>(it->second(t), resampler);
  }
  return out;
}

template<typename resampledCorrelationFunction, typename binResampler>
resampledCorrelationFunction computePiPi2ptVacSub(const bubbleDataAllMomenta &raw, const binResampler &resampler, const int tsep_pipi, 
						  const PiPiProjectorBase &proj_src, const PiPiProjectorBase &proj_snk){
  typedef typename resampledCorrelationFunction::DataType ResampledDistributionType;
  typedef bubbleDataAllMomentaSelect<ResampledDistributionType> bubbleDataAllMomentaType;

  auto dj_bubble_data = binResample<bubbleDataAllMomentaType>(raw, resampler);
  auto A2_realavg_V_dj = computePiPi2ptFigureVprojectSourceAvg(dj_bubble_data,tsep_pipi,proj_src,proj_snk);
  return resampledCorrelationFunction(3. * A2_realavg_V_dj);
}


//Vacuum subtraction for p_tot = (0,0,0)
template<typename resampledCorrelationFunction, typename binResampler>
inline resampledCorrelationFunction computePiPi2ptVacSub(const bubbleDataAllMomenta &raw, const binResampler &resampler, const int tsep_pipi, 
						   const PiPiProjector proj_src_t = PiPiProjector::A1momSet111, const PiPiProjector proj_snk_t = PiPiProjector::A1momSet111){
  std::unique_ptr<PiPiProjectorBase> proj_src( getProjector(proj_src_t, {0,0,0}) );
  std::unique_ptr<PiPiProjectorBase> proj_snk( getProjector(proj_snk_t, {0,0,0}) );
  return computePiPi2ptVacSub<resampledCorrelationFunction>(raw, resampler, tsep_pipi, *proj_src, *proj_snk);
}

//Non-zero p_tot
template<typename resampledCorrelationFunction, typename binResampler>
inline resampledCorrelationFunction computePiPi2ptVacSub(const bubbleDataAllMomenta &raw, const binResampler &resampler, const int tsep_pipi, const threeMomentum &p_tot,
						   const PiPiProjector proj_src_t = PiPiProjector::MovingSwaveGround, const PiPiProjector proj_snk_t = PiPiProjector::MovingSwaveGround){
  std::unique_ptr<PiPiProjectorBase> proj_src( getProjector(proj_src_t, p_tot) );
  std::unique_ptr<PiPiProjectorBase> proj_snk( getProjector(proj_snk_t, -p_tot) );
  return computePiPi2ptVacSub<resampledCorrelationFunction>(raw, resampler, tsep_pipi, *proj_src, *proj_snk);
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

struct getResampledPiPi2ptDataOpts{
  bool load_hdf5_data_checkpoint;
  std::string load_hdf5_data_checkpoint_stub;
  bool save_hdf5_data_checkpoint;
  std::string save_hdf5_data_checkpoint_stub;
  getResampledPiPi2ptDataOpts(): load_hdf5_data_checkpoint(false), save_hdf5_data_checkpoint(false){}
};


//Read and combine/double-jack resample data from original files or a checkpoint of the entire data set
template<typename FigureFilenamePolicy, typename BubbleFilenamePolicy>
void getResampledPiPi2ptData(jackknifeCorrelationFunctionD &pipi_j, doubleJackknifeCorrelationFunctionD &pipi_dj,
			     const PiPiProjectorBase &proj_src, const PiPiProjectorBase &proj_snk,
			     const int isospin, const int Lt, const int tsep_pipi, const int tstep_pipi, bool do_vacuum_subtraction,
			     const std::string &data_dir, const int traj_start, const int traj_inc, const int traj_lessthan, const int bin_size,
			     const FigureFilenamePolicy &ffn, const BubbleFilenamePolicy &bfn_src, const BubbleFilenamePolicy &bfn_snk, const std::string &extra_descr,
			     const getResampledPiPi2ptDataOpts &opts = getResampledPiPi2ptDataOpts()){

  //Read the data
  bubbleDataAllMomenta raw_bubble_data;
  rawDataCorrelationFunctionD pipi_raw;
  {
    figureDataAllMomenta raw_data;
    if(opts.load_hdf5_data_checkpoint){
      loadRawDataCheckpoint(raw_data, raw_bubble_data, opts.load_hdf5_data_checkpoint_stub, extra_descr);
    }else{
      std::cout << "Reading raw data " << extra_descr << std::endl;
      
      readRawPiPi2ptData(raw_data, raw_bubble_data, ffn, 
			 bfn_src, bfn_snk, data_dir, traj_start, traj_inc, traj_lessthan, 
			 Lt, tstep_pipi, tsep_pipi, proj_src, proj_snk);
    }
    if(opts.save_hdf5_data_checkpoint) saveRawDataCheckpoint(raw_data, raw_bubble_data, opts.save_hdf5_data_checkpoint_stub, extra_descr);
    combineRawPiPiContractions(pipi_raw, raw_data, proj_src, proj_snk, isospin, bin_size, extra_descr);
  }

  const int nsample = (traj_lessthan - traj_start)/traj_inc/bin_size;
  
  basicBinResampler resampler(bin_size);

  pipi_j = binResample<jackknifeCorrelationFunctionD>(pipi_raw, resampler);
  pipi_dj = binResample<doubleJackknifeCorrelationFunctionD>(pipi_raw, resampler);

  if(isospin == 0 && do_vacuum_subtraction){
    pipi_j = pipi_j - computePiPi2ptVacSub<jackknifeCorrelationFunctionD>(raw_bubble_data, resampler, tsep_pipi, proj_src, proj_snk);
    pipi_dj = pipi_dj - computePiPi2ptVacSub<doubleJackknifeCorrelationFunctionD>(raw_bubble_data, resampler, tsep_pipi, proj_src, proj_snk);
  }
  
  pipi_j = fold(pipi_j, 2*tsep_pipi);
  pipi_dj = fold(pipi_dj, 2*tsep_pipi);
}


CPSFIT_END_NAMESPACE

#endif
