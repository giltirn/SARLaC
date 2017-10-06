#ifndef _KTOPIPI_PLOT_H_
#define _KTOPIPI_PLOT_H_

#include <plot.h>

//Perform the error-weighted average of all data with fixed  tsep_op_pi  subject to a cut on the minimum tsep_k_op separation. Output coord is tsep_op_pi
void errorWeightedAverage(std::vector<correlationFunction<double, jackknifeDistributionD> > &out,
			  const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &in,
			  const int tmin_k_op){
  assert(in[0].size() > 0);
  const int nsample = in[0].value(0).size();
  
  for(int q=0;q<10;q++){
    std::map<int, std::vector<int> > dmap;
    for(int d=0;d<in[q].size();d++){
      const int tsep_k_op = int(in[q].coord(d).t);
      const int tsep_k_pi = in[q].coord(d).tsep_k_pi;
      const int tsep_op_pi = tsep_k_pi - tsep_k_op;
      if(tsep_op_pi >= 0 && tsep_k_op >= tmin_k_op) dmap[tsep_op_pi].push_back(d);
    }
    for(auto it = dmap.begin(); it != dmap.end(); it++){
      const int tsep_op_pi = it->first;
      const std::vector<int> &include = it->second;

      double wsum = 0.;
      jackknifeDistributionD wavg(nsample,0.);

      for(int dd=0; dd<include.size();dd++){
	const int d = include[dd];
	const double stderr = in[q].value(d).standardError();
	const double w = 1/stderr/stderr;
	wsum = wsum + w;
	wavg = wavg + w*in[q].value(d);
      }
      wavg = wavg / wsum;

      out[q].push_back(tsep_op_pi, wavg);
    }
  }
  
}

template<typename FitFunc>
struct extractMdata{};

//Separate Q fit
template<>
struct extractMdata<FitKtoPiPi>{
  const std::vector<jackknifeDistribution<FitKtoPiPi::Params> > &fit_params;

  extractMdata(const std::vector<jackknifeDistribution<FitKtoPiPi::Params> > &_fit_params): fit_params(_fit_params){}
  
  //Extract the matrix element from the data for a given coordinate assuming we know the remaining fit parameters
  double getMdata(const amplitudeDataCoord &x, const double y, const int q, const int sample) const{
    FitKtoPiPi::Params p1(fit_params[q].sample(sample)); p1.M = 1.;
    return y/FitKtoPiPi::value(x,p1);
  }

  inline double getMfit(const int q, const int s) const{ return fit_params[q].sample(s).M; }
};

//Simultaneous Q fit
template<>
struct extractMdata<FitKtoPiPiSim>{
  const std::vector<jackknifeDistribution<FitKtoPiPiSim::Params> > &fit_params;

  extractMdata(const std::vector<jackknifeDistribution<FitKtoPiPiSim::Params> > &_fit_params): fit_params(_fit_params){}
  
  //Extract the matrix element from the data for a given coordinate assuming we know the remaining fit parameters
  double getMdata(const amplitudeDataCoord &x, const double y, const int q, const int sample) const{
    FitKtoPiPiSim::Params p1(fit_params[0].sample(sample)); p1.M[q] = 1.;
    amplitudeDataCoordSim xx(x,q);
    return y/FitKtoPiPiSim::value(xx,p1);
  }

  inline double getMfit(const int q, const int s) const{ return fit_params[0].sample(s).M[q]; }
};


template<typename MdataExtractor>
void plotErrorWeightedData(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &data, const MdataExtractor &extractor, const Args &args){
  //Extract the matrix element from the data
  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > Mdata(10);
  assert(data[0].size() > 0);
  const int nsample = data[0].value(0).size();
    
  for(int q=0;q<10;q++)
    for(int x=0;x<data[q].size();x++){
      jackknifeDistributionD tmp = data[q].value(x);
      for(int s=0;s<nsample;s++) tmp.sample(s) = extractor.getMdata(data[q].coord(x), tmp.sample(s), q, s);
      Mdata[q].push_back(data[q].coord(x),tmp);
    }
  
  //Compute the weighted average
  std::vector<correlationFunction<double, jackknifeDistributionD> > wavg(10);
  errorWeightedAverage(wavg, Mdata, args.tmin_k_op);
  
  typedef DataSeriesAccessor<correlationFunction<double, jackknifeDistributionD>, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistributionD> > accessor;
  
  for(int q=0;q<10;q++){
    //Largest tsep_k_pi
    double largest_tsep_op_pi = -1;
    for(int i=0;i<wavg[q].size();i++) if(wavg[q].coord(i) > largest_tsep_op_pi) largest_tsep_op_pi = wavg[q].coord(i);

    //Setup the fit curve
    jackknifeDistributionD Mfit(nsample, [&](const int s){ return extractor.getMfit(q,s); });
    correlationFunction<double, jackknifeDistributionD> curve(2);
    curve.value(0) = Mfit;
    curve.coord(0) = args.tmin_op_pi;
    curve.value(1) = Mfit;
    curve.coord(1) = largest_tsep_op_pi;
        
    //Plot
    MatPlotLibScriptGenerate plotter;
    accessor dset_accessor(wavg[q]);
    accessor fitcurve_accessor(curve);
    MatPlotLibScriptGenerate::handleType dset_handle = plotter.plotData(dset_accessor);
    MatPlotLibScriptGenerate::handleType fitcurve_handle = plotter.errorBand(fitcurve_accessor);
    
    plotter.setXlabel("$t$");
    std::ostringstream ylabel; ylabel << "$M^{1/2,\\ \\rm{lat}}_" << q+1 << "$";
    plotter.setYlabel(ylabel.str());

    std::ostringstream filename_stub; filename_stub << "plot_errw_Q" << q+1;
    plotter.write( filename_stub.str()+".py", filename_stub.str()+".pdf");
  }
}


#endif
