#ifndef _KTOPIPI_PLOT_H_
#define _KTOPIPI_PLOT_H_

#include<set>

#include<config.h>
#include<utils/macros.h>
#include<utils/utils/error.h>

#include<common.h>
#include<plot.h>

#include "fitfunc.h"

CPSFIT_START_NAMESPACE

//Perform the error-weighted average of all data with fixed  tsep_op_pi  subject to a cut on the minimum tsep_k_op separation. Output coord is tsep_op_pi
template<typename DistributionType>
void errorWeightedAverage(std::vector<correlationFunction<double, DistributionType> > &out,
			  const std::vector<correlationFunction<amplitudeDataCoord, DistributionType> > &in,
			  const int tmin_k_op){
  assert(in[0].size() > 0);
  
  for(int q=0;q<in.size();q++){
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

      std::vector<DistributionType const*> ptrs(include.size());
      for(int dd=0;dd<include.size();dd++) ptrs[dd] = &in[q].value(include[dd]);
      
      DistributionType wavg = weightedAvg(ptrs);

      out[q].push_back(tsep_op_pi, wavg);
    }
  }
  
}

#define DIST template<typename, template<typename> class> class DistributionType 
#define VEC template<typename> class V = basic_vector


template<typename FitFunc, typename ParamsDist>
struct extractMdata{};

//Separate Q fit
template<typename ParamsDist>
struct extractMdata<FitKtoPiPi, ParamsDist>{
  typedef iterate<ParamsDist> iter;

  const std::vector<ParamsDist> &fit_params;

  extractMdata(const std::vector<ParamsDist> &_fit_params): fit_params(_fit_params){}
  
  //Extract the matrix element from the data for a given coordinate assuming we know the remaining fit parameters
  double getMdata(const amplitudeDataCoord &x, const double y, const int q, const int sample) const{
    FitKtoPiPi::Params p1(iter::at(sample, fit_params[q])); p1.M = 1.;
    return y/FitKtoPiPi::value(x,p1);
  }

  inline double getMfit(const int q, const int s) const{ return iter::at(s, fit_params[q]).M; }
};

//Simultaneous Q fit
template<int N, typename ParamsDist>
struct extractMdata<FitKtoPiPiSim<N>, ParamsDist>{
  typedef iterate<ParamsDist> iter;

  const ParamsDist &fit_params;

  extractMdata(const ParamsDist &_fit_params): fit_params(_fit_params){}
  
  //Extract the matrix element from the data for a given coordinate assuming we know the remaining fit parameters
  double getMdata(const amplitudeDataCoord &x, const double y, const int q, const int sample) const{
    typename FitKtoPiPiSim<N>::Params p1(iter::at(sample, fit_params)); p1.M[q] = 1.;
    amplitudeDataCoordSim xx(x,q);
    return y/FitKtoPiPiSim<N>::value(xx,p1);
  }

  inline double getMfit(const int q, const int s) const{ return iter::at(s, fit_params).M[q]; }
};
template<typename ParamsDist>
struct extractMdata<FitKtoPiPiSim<7>, ParamsDist>{ //These data are converted to the 10 basis prior to plotting
  typedef iterate<ParamsDist> iter;

  const ParamsDist &fit_params;

  extractMdata(const ParamsDist &_fit_params): fit_params(_fit_params){}
  
  //Extract the matrix element from the data for a given coordinate assuming we know the remaining fit parameters
  double getMdata(const amplitudeDataCoord &x, const double y, const int q, const int sample) const{
    typename FitKtoPiPiSim<10>::Params p1(iter::at(sample, fit_params)); p1.M[q] = 1.;
    amplitudeDataCoordSim xx(x,q);
    return y/FitKtoPiPiSim<10>::value(xx,p1);
  }

  inline double getMfit(const int q, const int s) const{ return iter::at(s, fit_params).M[q]; }
};


//Separate Q fit with constant
template<typename ParamsDist>
struct extractMdata<FitKtoPiPiWithConstant, ParamsDist>{
  typedef iterate<ParamsDist> iter;

  const std::vector<ParamsDist> &fit_params;

  extractMdata(const std::vector<ParamsDist> &_fit_params): fit_params(_fit_params){}
  
  //Extract the matrix element from the data for a given coordinate assuming we know the remaining fit parameters
  double getMdata(const amplitudeDataCoord &x, const double y, const int q, const int sample) const{
    FitKtoPiPiWithConstant::Params p1(iter::at(sample, fit_params[q])); p1.M = 1.;
    return y/FitKtoPiPiWithConstant::value(x,p1);
  }

  inline double getMfit(const int q, const int s) const{ return iter::at(s, fit_params[q]).M; }
};


template<typename ParamsDist>
struct extractMdata<FitKtoPiPiTwoExp, ParamsDist>{
  typedef iterate<ParamsDist> iter;

  const std::vector<ParamsDist> &fit_params;

  extractMdata(const std::vector<ParamsDist> &_fit_params): fit_params(_fit_params){}
  
  //Extract the matrix element from the data for a given coordinate assuming we know the remaining fit parameters
  double getMdata(const amplitudeDataCoord &x, const double y, const int q, const int sample) const{
    return FitKtoPiPiTwoExp::getM0data(y, x, iter::at(sample, fit_params[q]));
  }

  inline double getMfit(const int q, const int s) const{ return iter::at(s, fit_params[q]).M0; }
};


template<typename MdataExtractor, typename DistributionType>
void plotErrorWeightedData(const std::vector<correlationFunction<amplitudeDataCoord, DistributionType> > &data, const MdataExtractor &extractor, const int tmin_k_op, const int tmin_op_pi){
  typedef iterate<DistributionType> iter;

  //Extract the matrix element from the data
  const int nsample = iter::size(data[0].value(0));
  const int nq = data.size();

  std::vector<correlationFunction<amplitudeDataCoord, DistributionType> > Mdata(nq);
  assert(data[0].size() > 0);
  
  for(int q=0;q<nq;q++)
    for(int x=0;x<data[q].size();x++){
      DistributionType tmp = data[q].value(x);
      for(int s=0;s<nsample;s++) iter::at(s, tmp) = extractor.getMdata(data[q].coord(x), iter::at(s, tmp), q, s);
      Mdata[q].push_back(data[q].coord(x),tmp);
    }
  
  //Compute the weighted average
  std::vector<correlationFunction<double, DistributionType> > wavg(nq);
  errorWeightedAverage(wavg, Mdata, tmin_k_op);
 
  typedef DataSeriesAccessor<correlationFunction<double, DistributionType>, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<DistributionType> > accessor;
  
  for(int q=0;q<nq;q++){
    //Break up data by tsep_k_pi
    std::vector<int> tsep_k_pi_all;
    std::vector< correlationFunction<double, DistributionType> > tsep_k_pi_data;

    for(int i=0;i<Mdata[q].size();i++){
      int tsep_k_pi = Mdata[q].coord(i).tsep_k_pi;
      int top_pi = tsep_k_pi - (int)Mdata[q].coord(i).t;

      int idx = -1;
      for(int t=0;t<tsep_k_pi_all.size();t++) if(tsep_k_pi_all[t] == tsep_k_pi){ idx = t; break; }

      if(idx == -1){
	tsep_k_pi_all.push_back(tsep_k_pi);
	tsep_k_pi_data.push_back(correlationFunction<double, DistributionType>());
	idx = tsep_k_pi_all.size()-1;
      }

      tsep_k_pi_data[idx].push_back(top_pi, Mdata[q].value(i));
    }

    //Largest tsep_k_pi
    double largest_tsep_op_pi = -1;
    for(int i=0;i<wavg[q].size();i++) if(wavg[q].coord(i) > largest_tsep_op_pi) largest_tsep_op_pi = wavg[q].coord(i);

    //Setup the fit curve
    DistributionType Mfit(data[0].value(0).getInitializer());
    for(int s=0;s<iter::size(Mfit);s++) iter::at(s, Mfit) = extractor.getMfit(q,s);

    correlationFunction<double, DistributionType> curve(2);
    curve.value(0) = Mfit;
    curve.coord(0) = tmin_op_pi;
    curve.value(1) = Mfit;
    curve.coord(1) = largest_tsep_op_pi;
        
    //Plot
    MatPlotLibScriptGenerate plotter;

    for(int t=0;t<tsep_k_pi_all.size();t++)
      plotter.plotData(accessor(tsep_k_pi_data[t]),stringize("tsep_k_pi_%d", tsep_k_pi_all[t]));

    plotter.plotData(accessor(wavg[q]),"wavg");

    MatPlotLibScriptGenerate::kwargsType kwargs; kwargs["alpha"] = 0.6;
    plotter.errorBand(accessor(curve), kwargs, "fit");
    
    plotter.setXlabel("$t$");
    std::ostringstream ylabel; ylabel << "$M^{1/2,\\ \\rm{lat}}_" << q+1 << "$";
    plotter.setYlabel(ylabel.str());

    std::ostringstream filename_stub; filename_stub << "plot_errw_Q" << q+1;
    plotter.write( filename_stub.str()+".py", filename_stub.str()+".pdf");
  }
}

void plotErrorWeightedData2exp(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &data, 
			       const std::vector<jackknifeDistribution<FitKtoPiPiTwoExp::Params> > &fit_params,
			       const int Lt, const int tmin_k_op, const int tmin_op_pi){
  auto mK = distributionStructPeek(fit_params[0], &FitKtoPiPiTwoExp::Params::mK);
  FitKtoPiPiTwoExp fitfunc;
  auto zero = mK; zeroit(zero);

  for(int q=0;q<10;q++){
#define GETP(NM) auto NM = distributionStructPeek(fit_params[q], &FitKtoPiPiTwoExp::Params:: NM)
    GETP(AK);
    GETP(A0);
    GETP(M0);
    GETP(E0);
    GETP(A1);
    GETP(M1);
    GETP(E1);
#undef GETP

    //Time dep is currently
    //p.AK*p.A0*p.M0*exp(-p.E0 * c.tsep_k_pi)*exp( -(p.mK - p.E0)*c.t )/sqrt(2.) + 
    //  p.AK*p.A1*p.M1*exp(-p.E1 * c.tsep_k_pi)*exp( -(p.mK - p.E1)*c.t )/sqrt(2.);
    
    //Remove kaon time dependence from data

    typedef correlationFunction<amplitudeDataCoord, jackknifeDistributionD>::ElementType Etype;

    correlationFunction<amplitudeDataCoord, jackknifeDistributionD> data_tkrem = 
      correlationFunction<amplitudeDataCoord, jackknifeDistributionD>(data[q].size(),
								      [&](const int i){
									return Etype(data[q].coord(i),
										     data[q].value(i) * exp(mK * data[q].coord(i).t) );
								      }
								      );
    //Time dep now:
    //p.AK*p.A0*p.M0*exp(-p.E0 * (c.tsep_k_pi - c.t)/sqrt(2.) + 
    //  p.AK*p.A1*p.M1*exp(-p.E1 * (c.tsep_k_pi - c.t)/sqrt(2.);
    
    //tsep_k_pi - t = t_op_pi
    
    //Error weighted average over data with same t_op_pi
    std::vector< std::vector<int> > dmap(Lt);
    std::vector< std::vector<double> > wmap(Lt);
    std::vector<double> wsum(Lt,0.);

    int tsep_k_pi_max = -1;
    for(int i=0;i<data_tkrem.size();i++){
      if(data_tkrem.coord(i).t < tmin_k_op) continue;

      if((int)data_tkrem.coord(i).t > data_tkrem.coord(i).tsep_k_pi) continue; 

      if(data_tkrem.coord(i).tsep_k_pi > tsep_k_pi_max) tsep_k_pi_max = data_tkrem.coord(i).tsep_k_pi;

      int t_op_pi = data_tkrem.coord(i).tsep_k_pi - (int)data_tkrem.coord(i).t;
      double w = data_tkrem.value(i).standardError();
      w = 1./w/w;

      std::cout << "Q" << q+1 << " including data point " << i << " with t_op_pi " << t_op_pi << " and value " << data_tkrem.value(i) << " and weight " << w << ": running sum of weights for this t " << wsum[t_op_pi] << std::endl;

      dmap[t_op_pi].push_back(i);
      wmap[t_op_pi].push_back(w);
      wsum[t_op_pi] += w;
    }
    
    correlationFunction<double, jackknifeDistributionD> data_wavg;
    for(int t=0;t<tsep_k_pi_max;t++){
      if(dmap[t].size() > 0){
	jackknifeDistributionD v = wmap[t][0] * data_tkrem.value( dmap[t][0] );

	for(int i=1;i<dmap[t].size();i++)
	  v = v + wmap[t][i] * data_tkrem.value( dmap[t][i] );

	v = v / wsum[t];

	std::cout << "t_op_pi " << t << " " << " weighted avg of " << dmap[t].size() << " data (";
	for(int i=0;i<dmap[t].size();i++) std::cout << "[" << data_tkrem.coord( dmap[t][i] ) << "]";
	std::cout << ") to give result " << v << std::endl;

	data_wavg.push_back(t, v);
      }
    }
    
    //Construct fit curves
    correlationFunction<double, jackknifeDistributionD> curve_ground;
    correlationFunction<double, jackknifeDistributionD> curve_excited;
    correlationFunction<double, jackknifeDistributionD> curve_sum;

    int npoint = 60;
    double delta = double(tsep_k_pi_max-tmin_op_pi)/(npoint - 1);
    for(int i=0;i<npoint;i++){
      double t_op_pi = tmin_op_pi + i*delta;
      jackknifeDistributionD gnd = AK*A0*M0*exp(-E0*t_op_pi)/sqrt(2.);
      jackknifeDistributionD exc = AK*A1*M1*exp(-E1*t_op_pi)/sqrt(2.);
      jackknifeDistributionD sum = gnd+exc;
      curve_ground.push_back(t_op_pi, gnd);
      curve_excited.push_back(t_op_pi, exc);
      curve_sum.push_back(t_op_pi, sum);
    }    
  
    //Plot
    MatPlotLibScriptGenerate plotter;
    typedef DataSeriesAccessor<correlationFunction<double, jackknifeDistributionD>, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistributionD> > accessor;
    accessor dset_accessor(data_wavg);
    accessor fitcurve_gnd_accessor(curve_ground);
    accessor fitcurve_exc_accessor(curve_excited);
    accessor fitcurve_sum_accessor(curve_sum);

    MatPlotLibScriptGenerate::handleType handle;
    MatPlotLibScriptGenerate::kwargsType kwargs;

    handle = plotter.plotData(dset_accessor);
    kwargs["alpha"] = 0.5;
    kwargs["color"] = 'r';
    handle = plotter.errorBand(fitcurve_gnd_accessor,kwargs);
    kwargs["color"] = 'g';
    handle = plotter.errorBand(fitcurve_exc_accessor,kwargs);
    kwargs["color"] = 'b';
    handle = plotter.errorBand(fitcurve_sum_accessor,kwargs);
    
    plotter.setXlabel("$t$");
    std::ostringstream ylabel; ylabel << "$M^{1/2,\\ \\rm{lat}}_" << q+1 << "$";
    plotter.setYlabel(ylabel.str());

    std::ostringstream filename_stub; filename_stub << "plot_2exp_errw_Q" << q+1;
    plotter.write( filename_stub.str()+".py", filename_stub.str()+".pdf");
  }
}



void plotErrorWeightedData2expFlat(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &data, 
			       const std::vector<jackknifeDistribution<FitKtoPiPiTwoExp::Params> > &fit_params,
			       const int Lt, const int tmin_k_op, const int tmin_op_pi){
  auto mK = distributionStructPeek(fit_params[0], &FitKtoPiPiTwoExp::Params::mK);
  FitKtoPiPiTwoExp fitfunc;
  auto zero = mK; zeroit(zero);

  for(int q=0;q<10;q++){
#define GETP(NM) auto NM = distributionStructPeek(fit_params[q], &FitKtoPiPiTwoExp::Params:: NM)
    GETP(AK);
    GETP(A0);
    GETP(M0);
    GETP(E0);
    GETP(A1);
    GETP(M1);
    GETP(E1);
#undef GETP

    //Time dep is currently
    //p.AK*p.A0*p.M0*exp(-p.E0 * c.tsep_k_pi)*exp( -(p.mK - p.E0)*c.t )/sqrt(2.) + 
    //  p.AK*p.A1*p.M1*exp(-p.E1 * c.tsep_k_pi)*exp( -(p.mK - p.E1)*c.t )/sqrt(2.);
    
    //Remove kaon and ground state pipi time dependence from data

    typedef correlationFunction<amplitudeDataCoord, jackknifeDistributionD>::ElementType Etype;

    correlationFunction<amplitudeDataCoord, jackknifeDistributionD> data_tkrem = 
      correlationFunction<amplitudeDataCoord, jackknifeDistributionD>(data[q].size(),
								      [&](const int i){
									const amplitudeDataCoord &coord = data[q].coord(i); 
									int t_op_pi = coord.tsep_k_pi - (int)coord.t;
									return Etype(coord,
										     data[q].value(i) * exp(mK * coord.t) * exp(E0*t_op_pi) );
								      }
								      );
    //Time dep now:
    //p.AK*p.A0*p.M0/sqrt(2.) + 
    //  p.AK*p.A1*p.M1*exp(-(p.E1-p.E0) * (c.tsep_k_pi - c.t))/sqrt(2.);
    
    //tsep_k_pi - t = t_op_pi
    
    //Error weighted average over data with same t_op_pi
    std::vector< std::vector<int> > dmap(Lt);
    std::vector< std::vector<double> > wmap(Lt);
    std::vector<double> wsum(Lt,0.);

    std::set<int> tsep_k_pi_set;

    int tsep_k_pi_max = -1;
    for(int i=0;i<data_tkrem.size();i++){
      if(data_tkrem.coord(i).t < tmin_k_op) continue;

      if((int)data_tkrem.coord(i).t > data_tkrem.coord(i).tsep_k_pi) continue; 

      if(data_tkrem.coord(i).tsep_k_pi > tsep_k_pi_max) tsep_k_pi_max = data_tkrem.coord(i).tsep_k_pi;

      tsep_k_pi_set.insert(data_tkrem.coord(i).tsep_k_pi);

      int t_op_pi = data_tkrem.coord(i).tsep_k_pi - (int)data_tkrem.coord(i).t;
      double w = data_tkrem.value(i).standardError();
      w = 1./w/w;

      std::cout << "Q" << q+1 << " including data point " << i << " with t_op_pi " << t_op_pi << " and value " << data_tkrem.value(i) << " and weight " << w << ": running sum of weights for this t " << wsum[t_op_pi] << std::endl;

      dmap[t_op_pi].push_back(i);
      wmap[t_op_pi].push_back(w);
      wsum[t_op_pi] += w;
    }

    
    std::map<int, int> tsep_k_pi_idx_map;
    std::map<int, int> idx_tsep_k_pi_idx_map;
    int idx = 0;
    for(auto it = tsep_k_pi_set.begin(); it != tsep_k_pi_set.end(); ++it){
      idx_tsep_k_pi_idx_map[idx] = *it;
      tsep_k_pi_idx_map[*it] = idx++;
    }
    
    correlationFunction<double, jackknifeDistributionD> data_wavg;
    std::vector<correlationFunction<double, jackknifeDistributionD> > data_tsepkpi(tsep_k_pi_set.size());

    for(int t=0;t<tsep_k_pi_max;t++){
      if(dmap[t].size() > 0){

	//Add to separate tsep_k_pi sets
	for(int i=0;i<dmap[t].size();i++){
	  int tsep_k_pi = data_tkrem.coord( dmap[t][i] ).tsep_k_pi;
	  data_tsepkpi[ tsep_k_pi_idx_map[tsep_k_pi] ].push_back( t, data_tkrem.value( dmap[t][i] ) );
	}
	
	//Construct weighted avg
	jackknifeDistributionD v = wmap[t][0] * data_tkrem.value( dmap[t][0] );

	for(int i=1;i<dmap[t].size();i++)
	  v = v + wmap[t][i] * data_tkrem.value( dmap[t][i] );

	v = v / wsum[t];

	std::cout << "t_op_pi " << t << " " << " weighted avg of " << dmap[t].size() << " data (";
	for(int i=0;i<dmap[t].size();i++) std::cout << "[" << data_tkrem.coord( dmap[t][i] ) << ":" << data_tkrem.value(dmap[t][i]) << "]";
	std::cout << ") to give result " << v << std::endl;

	data_wavg.push_back(t, v);
      }
    }
    
    //Construct fit curves
    correlationFunction<double, jackknifeDistributionD> curve_ground;
    correlationFunction<double, jackknifeDistributionD> curve_excited;
    correlationFunction<double, jackknifeDistributionD> curve_sum;

    int npoint = 60;
    double delta = double(tsep_k_pi_max-tmin_op_pi)/(npoint - 1);
    for(int i=0;i<npoint;i++){
      double t_op_pi = tmin_op_pi + i*delta;
      jackknifeDistributionD gnd = AK*A0*M0/sqrt(2.);
      jackknifeDistributionD exc = AK*A1*M1*exp(-(E1-E0)*t_op_pi)/sqrt(2.);
      jackknifeDistributionD sum = gnd+exc;
      curve_ground.push_back(t_op_pi, gnd);
      curve_excited.push_back(t_op_pi, exc);
      curve_sum.push_back(t_op_pi, sum);
    }    
  
    {
      //Plot
      MatPlotLibScriptGenerate plotter;
      typedef DataSeriesAccessor<correlationFunction<double, jackknifeDistributionD>, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistributionD> > accessor;
     
      MatPlotLibScriptGenerate::handleType handle;
      MatPlotLibScriptGenerate::kwargsType kwargs;

      handle = plotter.plotData(accessor(data_wavg),"wavg");
      plotter.setLegend(handle,"Weighted avg");
      
      for(int i=0;i<data_tsepkpi.size();i++){
	handle = plotter.plotData(accessor(data_tsepkpi[i]), stringize("tsep_k_pi%d",idx_tsep_k_pi_idx_map[i]));
	plotter.setLegend(handle,stringize("$t_{\\rm sep}^{K\\to\\pi}=%d$",idx_tsep_k_pi_idx_map[i]));
      }      

      kwargs["alpha"] = 0.5;
      kwargs["color"] = 'r';
      handle = plotter.errorBand(accessor(curve_ground),kwargs,"fit_gnd");
      plotter.setLegend(handle,"Gnd");
      kwargs["color"] = 'g';
      handle = plotter.errorBand(accessor(curve_excited),kwargs,"fit_exc");
      plotter.setLegend(handle,"Exc");
      kwargs["color"] = 'b';
      handle = plotter.errorBand(accessor(curve_sum),kwargs,"fit_sum");
      plotter.setLegend(handle,"Tot");
      
      plotter.setXlabel("$t$");
      std::ostringstream ylabel; ylabel << "$M^{1/2,\\ \\rm{lat}}_" << q+1 << "$";
      plotter.setYlabel(ylabel.str());
      
      plotter.createLegend();
      std::ostringstream filename_stub; filename_stub << "plot_2exp_errw_Q" << q+1 << "_flat";
      plotter.write( filename_stub.str()+".py", filename_stub.str()+".pdf");
    }

    //Histogram plot of data
    {
      MatPlotLibScriptGenerate plotter;
      std::ostream &preinv = plotter.preinvoke();
      std::ostream &inv = plotter.invoke();

      MatPlotLibScriptGenerate::kwargsType kwargs;
      kwargs["nbins"] = 50;

      for(int i=0;i<data_wavg.size();i++){
	int t = (int)data_wavg.coord(i);
	DistributionSampleAccessor<jackknifeDistributionD> accessor(data_wavg.value(i));
	plotter.histogram(accessor, kwargs, stringize("t_op_pi%d",t));

	double err = data_wavg.value(i).standardError()/sqrt(data_wavg.value(i).size()-1); //remove jackknife scaling factor
	double cen = data_wavg.value(i).best();
	double lo = cen - err;
	double hi = cen + err;
	
	std::string dist_info_var = stringize("distinfo_t_op_pi%d",t);
	
	preinv << dist_info_var << "=[" << hi << ", " << cen << ", " << lo << "]\n"; 
	inv << '\t' << "pyplot.plt.axvline(x=" << dist_info_var << "[0])\n";
	inv << '\t' << "pyplot.plt.axvline(x=" << dist_info_var << "[1])\n";
	inv << '\t' << "pyplot.plt.axvline(x=" << dist_info_var << "[2])\n";
      }
      plotter.setXlabel("sample");
      plotter.setYlabel("$N$");
      
      std::ostringstream filename_stub; filename_stub << "plot_2exp_errw_Q" << q+1 << "_flat_jackhist";
      plotter.write( filename_stub.str()+".py", filename_stub.str()+".pdf");
    }
  }
}




template<typename FitFunc, typename DistributionD>
struct plotFF{
  inline static void plot(const std::vector<correlationFunction<amplitudeDataCoord, DistributionD> > &A0_all_j,
			  const typename fitReturnType<FitFunc, DistributionD>::type &fit_params,
			  const int Lt, const int tmin_k_op, const int tmin_op_pi){
    extractMdata<FitFunc, typename fitReturnType<FitFunc, DistributionD>::ParamsDist> extractor(fit_params);
    plotErrorWeightedData(A0_all_j,extractor,tmin_k_op,tmin_op_pi);
  }
};
template<>
struct plotFF<FitKtoPiPiTwoExp, jackknifeDistributionD>{
  inline static void plot(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_all_j,
			  const std::vector<jackknifeDistribution<FitKtoPiPiTwoExp::Params> > &fit_params,
			  const int Lt, const int tmin_k_op, const int tmin_op_pi){
    int nsample = A0_all_j[0].value(0).size();
    plotErrorWeightedData2exp(A0_all_j, fit_params, Lt, tmin_k_op, tmin_op_pi);
    plotErrorWeightedData2expFlat(A0_all_j, fit_params, Lt, tmin_k_op, tmin_op_pi);
  }
};

struct PlotOnlyInputs{
  const int Lt;
  const int tmin_k_op;
  const int tmin_op_pi;

  std::string load_amplitude_data_file;

  PlotOnlyInputs(const int Lt, const int tmin_k_op, const int tmin_op_pi, const std::string &load_amplitude_data_file): Lt(Lt), tmin_k_op(tmin_k_op), tmin_op_pi(tmin_op_pi), load_amplitude_data_file(load_amplitude_data_file){}
};
template<typename FitFunc>
struct PlotOnlyCall{
  static void call(const PlotOnlyInputs &inputs){
    if(!fileExists(inputs.load_amplitude_data_file)) error_exit(std::cout << "PlotOnlyCall::call Must provide a valid checkpoint file, got \"" << inputs.load_amplitude_data_file << "\"\n");
    std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > A0_all_j(10);
#ifdef HAVE_HDF5
    {
      HDF5reader reader(inputs.load_amplitude_data_file);
      read(reader, A0_all_j, "A0_all_j");
    }
#else
    assert(0);
#endif

    typename fitReturnType<FitFunc, jackknifeDistributionD>::type fit_params;
#ifdef HAVE_HDF5
    readParamsStandard(fit_params, "params.hdf5");
#else
    assert(0);
#endif
    plotFF<FitFunc, jackknifeDistributionD>::plot(A0_all_j, fit_params, inputs.Lt, inputs.tmin_k_op, inputs.tmin_op_pi);
  }
    
};

CPSFIT_END_NAMESPACE

#endif
