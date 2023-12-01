#ifndef _PIPI_RAW_CORRELATOR_H_
#define _PIPI_RAW_CORRELATOR_H_

#include<config.h>
#include<utils/macros.h>

#include <data_series.h>
#include <common.h>

SARLAC_START_NAMESPACE

//Generate a correlation function by source-averaging the raw data
template<typename FigureDataType>
auto sourceAverage(const FigureDataType & data)->correlationFunction<double,typename FigureDataType::DistributionType>{
  const int Lt = data.getLt();

  std::vector<int> tsrc_include;
  for(int tsrc=0;tsrc<Lt;tsrc++){
    bool is_nonzero = !data.isZero(tsrc);
    if(is_nonzero)
      tsrc_include.push_back(tsrc);
  }
  const double N(tsrc_include.size());

  std::cout << "sourceAverage detected " << N << " non-zero timeslices\n";

  correlationFunction<double, typename FigureDataType::DistributionType> into(Lt);

  for(int tsep=0;tsep<Lt;tsep++){
    into.coord(tsep) = tsep;
    auto & v = into.value(tsep);
    v = data(tsrc_include[0],tsep);
    for(int i=1;i<tsrc_include.size();i++)
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

inline void bin(rawDataCorrelationFunctionD &raw, const int bin_size){
  for(int i=0;i<raw.size();i++) raw.value(i) = raw.value(i).bin(bin_size);
}


SARLAC_END_NAMESPACE

#endif
