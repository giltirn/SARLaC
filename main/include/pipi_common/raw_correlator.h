#ifndef _PIPI_RAW_CORRELATOR_H_
#define _PIPI_RAW_CORRELATOR_H_

#include<config.h>
#include<utils/macros.h>

#include "mom_data_containers.h"
#include "mom_project.h"
#include "raw_data.h"

CPSFIT_START_NAMESPACE

//Generate a correlation function by source-averaging the raw data
template<typename FigureDataType>
auto sourceAverage(const FigureDataType & data)->correlationFunction<double,typename std::decay<decltype(data(0,0))>::type>{
  typedef typename std::decay<decltype(data(0,0))>::type DistributionType;

  int Lt = data.getLt();
  int nsample = data.getNsample();
  correlationFunction<double,DistributionType> into(data.getLt(),
					     [nsample](int i) {  return typename correlationFunction<double,DistributionType>::ElementType(i, DistributionType(nsample,0.)); }
					     );
  std::vector<int> tsrc_include;
  for(int tsrc=0;tsrc<Lt;tsrc++){
    bool is_nonzero = !data.isZero(tsrc);
    if(is_nonzero)
      tsrc_include.push_back(tsrc);
  }
  double N(tsrc_include.size());

  std::cout << "sourceAverage detected " << N << " non-zero timeslices\n";

  for(int tsep=0;tsep<Lt;tsep++){
    auto & v = into.value(tsep);
    v.zero();
    for(int i=0;i<tsrc_include.size();i++)
      v = v + data(tsrc_include[i],tsep);
    v = v/N;
  }
  return into;
}

void outputRawCorrelator(const std::string &filename, const correlationFunction<double,rawDataDistributionD> &data, const double coeff){
  std::ofstream of(filename.c_str());
  of << std::setprecision(11) << std::scientific;
  int Lt = data.size();
  int nsample = data.value(0).size();

  for(int t=0;t<data.size();t++)
    for(int s=0;s<nsample;s++)
      of << t << " " << s << " " << coeff * data.value(t).sample(s) << " " << 0. << std::endl;
  of.close();
}

typedef correlationFunction<double,rawDataDistributionD> rawCorrelationFunction;

inline void bin(rawCorrelationFunction &raw, const int bin_size){
  for(int i=0;i<raw.size();i++) raw.value(i) = raw.value(i).bin(bin_size);
}

//Given the parsed, raw data, compute the raw , unbinned, unresampled pipi correlation function from the underlying contraction data. This includes projecting the pipi states onto
//a user-selected linear combination (for example projecting onto the A1 cubic representation)
void getRawPiPiCorrFunc(rawCorrelationFunction &pipi_raw, const figureDataAllMomenta &raw_data,
			const PiPiProjectorBase &proj_src, const PiPiProjectorBase &proj_snk, const int isospin, 
			const int bin_size, const std::string &extra_descr = "", bool output_raw_data = true){
  if(isospin != 0 && isospin != 2) error_exit(std::cout << "getRawPiPiCorrFunc only supports isospin 0,2\n");
  const char figs[4] = {'C','D','R','V'};
  const std::vector<double> coeffs = isospin == 0 ? std::vector<double>({1., 2., -6., 3.}) : std::vector<double>({-2., 2., 0., 0.});
  std::string ee = extra_descr != "" ? "_" + extra_descr : "";  

  rawCorrelationFunction fig_corr[4];  

  for(int f=0;f<4;f++){
    figureData proj_data = project(figs[f], raw_data, proj_src, proj_snk);
    fig_corr[f] = sourceAverage(proj_data);

    if(output_raw_data){
      //These data are saved after binning
      rawCorrelationFunction realavg_b(fig_corr[f]);
      bin(realavg_b, bin_size);
      outputRawCorrelator(stringize("raw_data_%cpart%s.dat",figs[f],ee.c_str()), realavg_b, coeffs[f]);
    }
  }    
  
  pipi_raw = coeffs[0]*fig_corr[0] + coeffs[1]*fig_corr[1] + coeffs[2]*fig_corr[2] + coeffs[3]*fig_corr[3];
  std::cout << "Raw data " << extra_descr << ":\n" << pipi_raw << std::endl;
}




CPSFIT_END_NAMESPACE

#endif
