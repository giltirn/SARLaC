#ifndef _FIT_KTOPIPI_KTOSIGMA_GPARITY_SIMFIT_PLOT_H
#define _FIT_KTOPIPI_KTOSIGMA_GPARITY_SIMFIT_PLOT_H

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

template< template<typename, template<typename> class> class DistributionType > 
inline DistributionType<double, basic_vector> peek(const std::string &tag, 
						   const DistributionType<taggedValueContainer<double,std::string>, basic_vector > &params){
  typedef iterate< DistributionType<taggedValueContainer<double,std::string>, basic_vector > > iter_p;
  typedef iterate< DistributionType<double, basic_vector> > iter_v;
  int index = params.sample(0).index(tag);
  DistributionType<double, basic_vector> out(params.getInitializer());
  for(int s=0;s<iter_p::size(params);s++) iter_v::at(s, out) = iter_p::at(s, params)(index);
  return out;
}

template<typename DistributionType>
void filterAndErrWeightedAvgData(correlationFunction<double, DistributionType> &data_wavg, //data that has been weighted avg
				 std::vector<correlationFunction<double, DistributionType> > &data_tsepkpi, //for each tsep_k_pi index, the data
				 std::map<int, int> &tsep_k_pi_idx_to_tsep_k_pi_map, //map of tsep_k_pi index to tsep_k_pi
				 const correlationFunction<amplitudeDataCoord, DistributionType> &data_in,
				 const int tmin_k_op, const int Lt, const std::string &descr){
  //Pick out data with t >= tmin_k_op and compute weights
  std::vector< std::vector<int> > dmap(Lt); //indices  [t_op_pi]
  std::vector< std::vector<double> > wmap(Lt); //weights
  std::vector<double> wsum(Lt,0.); //sum of weights

  std::set<int> tsep_k_pi_set; //set of tsep_k_op available

  int tsep_k_pi_max = -1;
  for(int i=0;i<data_in.size();i++){
    if(data_in.coord(i).t < tmin_k_op) continue;

    if((int)data_in.coord(i).t > data_in.coord(i).tsep_k_pi) continue; 

    if(data_in.coord(i).tsep_k_pi > tsep_k_pi_max) tsep_k_pi_max = data_in.coord(i).tsep_k_pi;

    tsep_k_pi_set.insert(data_in.coord(i).tsep_k_pi);

    int t_op_pi = data_in.coord(i).tsep_k_pi - (int)data_in.coord(i).t;
    double w = data_in.value(i).standardError();
    w = 1./w/w;

    std::cout << descr << " including data point " << i << " with t_op_pi " << t_op_pi << "(tk_pi " << data_in.coord(i).tsep_k_pi << ") and value " << data_in.value(i) << " and weight " << w << ": running sum of weights for this t " << wsum[t_op_pi] << std::endl;

    dmap[t_op_pi].push_back(i);
    wmap[t_op_pi].push_back(w);
    wsum[t_op_pi] += w;
  }

  //Map available tsep_k_pi to an index from 0 to the number of available tsep_k_pi
  std::map<int, int> tsep_k_pi_idx_map;
  int idx = 0;
  std::cout << "t_k_pi index mapping:" << std::endl;
  for(auto it = tsep_k_pi_set.begin(); it != tsep_k_pi_set.end(); ++it){
    std::cout << idx << " <-> " << *it << std::endl;
    tsep_k_pi_idx_to_tsep_k_pi_map[idx] = *it; //idx -> tsep_k_pi
    tsep_k_pi_idx_map[*it] = idx++; //tsep_k_pi -> idx
  }
    
  //Sort data and construct error weighted averages
  data_tsepkpi.resize(tsep_k_pi_set.size());

  for(int t=0;t<tsep_k_pi_max;t++){
    if(dmap[t].size() > 0){
      
      //Add to separate tsep_k_pi sets
      for(int i=0;i<dmap[t].size();i++){
	int tsep_k_pi = data_in.coord( dmap[t][i] ).tsep_k_pi;
	data_tsepkpi[ tsep_k_pi_idx_map[tsep_k_pi] ].push_back( t, data_in.value( dmap[t][i] ) );
      }
      
      //Construct weighted avg
      DistributionType v = wmap[t][0] * data_in.value( dmap[t][0] );
      
      for(int i=1;i<dmap[t].size();i++)
	v = v + wmap[t][i] * data_in.value( dmap[t][i] );
      
      v = v / wsum[t];
      
      std::cout << "t_op_pi " << t << " " << " weighted avg of " << dmap[t].size() << " data (";
      for(int i=0;i<dmap[t].size();i++) std::cout << "[" << data_in.coord( dmap[t][i] ) << ":" << data_in.value(dmap[t][i]) << "]";
      std::cout << ") to give result " << v << std::endl;
      
      data_wavg.push_back(t, v);
    }
  }
}

template<typename DistributionType>
inline correlationFunction<amplitudeDataCoord, DistributionType> removeGroundStateDependence(const correlationFunction<amplitudeDataCoord, DistributionType> &in_data,
											     const DistributionType &AK, const DistributionType &mK,
											     const DistributionType &Aop0, const DistributionType &E0){
    typedef correlationFunction<amplitudeDataCoord, DistributionType> CorrFuncType;
    typedef typename CorrFuncType::ElementType Etype;

    return CorrFuncType(in_data.size(),
			[&](const int i){
			  const amplitudeDataCoord &coord = in_data.coord(i); 
			  int t_op_pi = coord.tsep_k_pi - (int)coord.t;
			  return Etype(coord,
				       in_data.value(i) * exp(mK * coord.t) * exp(E0*t_op_pi) * sqrt(2.) / AK / Aop0 );
			}
			);
}



void plotErrorWeightedData2expFlat(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &ktopipi_data,
				   const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &ktosigma_data,
				   const std::vector<jackknifeDistribution<FitSimGenTwoState::Params> > &fit_params,
				   const int Lt, const int tmin_k_op, const int tmin_op_snk,
				   const FitSimGenTwoState &fitfunc){
  for(int q=0;q<10;q++){
#define GETP(NM) auto NM = peek(#NM, fit_params[q])
    GETP(AK);
    GETP(mK);
    GETP(Apipi0);
    GETP(Apipi1);
    GETP(Asigma0);
    GETP(Asigma1);
    GETP(E0);
    GETP(E1);
    GETP(M0);
    GETP(M1);
#undef GETP

    auto zero = mK; zeroit(zero);

    //Time dep is currently
    //AK*A0*M0*exp(-E0 * tsep_k_pi)*exp( -(mK - E0)*t )/sqrt(2.) + 
    //  AK*A1*M1*exp(-E1 * tsep_k_pi)*exp( -(mK - E1)*t )/sqrt(2.);
    
    //Remove kaon and ground state pipi time dependence from data

    typedef correlationFunction<amplitudeDataCoord, jackknifeDistributionD> CorrFuncType;

    CorrFuncType ktopipi_flat = removeGroundStateDependence(ktopipi_data[q], AK, mK, Apipi0, E0);
    CorrFuncType ktosigma_flat = removeGroundStateDependence(ktosigma_data[q], AK, mK, Asigma0, E0);

    //Time dep now:
    //M0 + (A1/A0)*M1*exp(-(E1-E0) * (tsep_k_pi - t));
    
    //tsep_k_pi - t = t_op_pi
    
    //Error weighted average over data with same t_op_pi
    typedef correlationFunction<double, jackknifeDistributionD> PlotCorrFuncType;
    PlotCorrFuncType ktopipi_wavg, ktosigma_wavg;
    std::vector<PlotCorrFuncType> ktopipi_tsepkop, ktosigma_tsepkop;
    std::map<int, int> ktopipi_idx_tsep_k_op_map, ktosigma_idx_tsep_k_op_map;
  

    filterAndErrWeightedAvgData(ktopipi_wavg, ktopipi_tsepkop, ktopipi_idx_tsep_k_op_map, ktopipi_flat, tmin_k_op, Lt, stringize("K->pipi Q%d",q+1) );
    filterAndErrWeightedAvgData(ktosigma_wavg, ktosigma_tsepkop, ktosigma_idx_tsep_k_op_map, ktosigma_flat, tmin_k_op, Lt, stringize("K->sigma Q%d",q+1));

    int tsep_k_op_max = -1;
    for(auto it=ktopipi_idx_tsep_k_op_map.begin(); it != ktopipi_idx_tsep_k_op_map.end(); ++it) tsep_k_op_max = std::max(tsep_k_op_max, it->second);
    for(auto it=ktosigma_idx_tsep_k_op_map.begin(); it != ktosigma_idx_tsep_k_op_map.end(); ++it) tsep_k_op_max = std::max(tsep_k_op_max, it->second);
    
    //Construct fit curves
    PlotCorrFuncType curve_ground, curve_excited_ktopipi, curve_sum_ktopipi, curve_excited_ktosigma, curve_sum_ktosigma;

    int npoint = 60;
    double delta = double(tsep_k_op_max-tmin_op_snk)/(npoint - 1);
    for(int i=0;i<npoint;i++){
      double t_op_snk = tmin_op_snk + i*delta;
      jackknifeDistributionD gnd = M0;
      curve_ground.push_back(t_op_snk, gnd);

      jackknifeDistributionD exc = Apipi1*M1*exp(-(E1-E0)*t_op_snk)/Apipi0;
      jackknifeDistributionD sum = gnd+exc;

      curve_excited_ktopipi.push_back(t_op_snk, exc);
      curve_sum_ktopipi.push_back(t_op_snk, sum);

      exc = Asigma1*M1*exp(-(E1-E0)*t_op_snk)/Asigma0;
      sum = gnd+exc;

      curve_excited_ktosigma.push_back(t_op_snk, exc);
      curve_sum_ktosigma.push_back(t_op_snk, sum);
    }
  
    {
      //Plot
      MatPlotLibScriptGenerate plotter;
      typedef DataSeriesAccessor<PlotCorrFuncType, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistributionD> > accessor;
     
      //Weighted avgs
      plotter.setLegend(plotter.plotData(accessor(ktopipi_wavg),{{"color",'r'}},"ktopipi_wavg"),
			"$K\\to\\pi\\pi$ weighted avg");

      plotter.setLegend(plotter.plotData(accessor(ktosigma_wavg),{{"color",'b'}},"ktosigma_wavg"),
			"$K\\to\\sigma$ weighted avg");

      //Individual tseps
      for(int i=0;i<ktopipi_tsepkop.size();i++){
	int tsep_k_op =  ktopipi_idx_tsep_k_op_map[i];
	plotter.setLegend(plotter.plotData(accessor(ktopipi_tsepkop[i]), {{"color",'r'}}, stringize("ktopipi_tsep_k_pi%d",tsep_k_op)),
			  stringize("$K\\to\\pi\\pi$ $t_{\\rm sep}^{K\\to{\\rm op}}=%d$",tsep_k_op));
      }      
      for(int i=0;i<ktosigma_tsepkop.size();i++){
	int tsep_k_op =  ktosigma_idx_tsep_k_op_map[i];
	plotter.setLegend(plotter.plotData(accessor(ktosigma_tsepkop[i]), {{"color",'b'}}, stringize("ktosigma_tsep_k_pi%d",tsep_k_op)),
			  stringize("$K\\to\\sigma$ $t_{\\rm sep}^{K\\to{\\rm op}}=%d$",tsep_k_op));
      }      

      plotter.setLegend(plotter.errorBand(accessor(curve_ground),{{"color",'r'}, {"alpha",0.5}},"fit_gnd"),
			"Gnd");

      plotter.setLegend(plotter.errorBand(accessor(curve_excited_ktopipi),{{"color",'g'}, {"alpha",0.5}},"fit_exc_ktopipi"),
			"$K\\to\\pi\\pi$ Exc");

      plotter.setLegend(plotter.errorBand(accessor(curve_sum_ktopipi),{{"color",'b'}, {"alpha",0.5}},"fit_sum_ktopipi"),
			"$K\\to\\pi\\pi$ Tot");

      plotter.setLegend(plotter.errorBand(accessor(curve_excited_ktosigma),{{"color",'g'}, {"alpha",0.5}},"fit_exc_ktosigma"),
			"$K\\to\\sigma$ Exc");

      plotter.setLegend(plotter.errorBand(accessor(curve_sum_ktosigma),{{"color",'b'}, {"alpha",0.5}},"fit_sum_ktosigma"),
			"$K\\to\\sigma$ Tot");
      
      plotter.setXlabel("$t$");
      std::ostringstream ylabel; ylabel << "$M^{1/2,\\ \\rm{lat}}_" << q+1 << "$";
      plotter.setYlabel(ylabel.str());
      
      plotter.createLegend();
      std::ostringstream filename_stub; filename_stub << "plot_2exp_errw_Q" << q+1 << "_flat";
      plotter.write( filename_stub.str()+".py", filename_stub.str()+".pdf");
    }
  }
}



void plotErrorWeightedData2expFlatDiffM1(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &ktopipi_data,
					 const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &ktosigma_data,
					 const std::vector<jackknifeDistribution<FitSimGenTwoState::Params> > &fit_params,
					 const int Lt, const int tmin_k_op, const int tmin_op_snk,
					 const FitSimGenTwoState &fitfunc){
  for(int q=0;q<10;q++){
#define GETP(NM) auto NM = peek(#NM, fit_params[q])
    GETP(AK);
    GETP(mK);
    GETP(Apipi0);
    GETP(Apipi1);
    GETP(Asigma0);
    GETP(Asigma1);
    GETP(E0);
    GETP(E1);
    GETP(M0);
    GETP(M1pipi);
    GETP(M1sigma);
#undef GETP

    auto zero = mK; zeroit(zero);

    //Time dep is currently
    //AK*A0*M0*exp(-E0 * tsep_k_pi)*exp( -(mK - E0)*t )/sqrt(2.) + 
    //  AK*A1*M1*exp(-E1 * tsep_k_pi)*exp( -(mK - E1)*t )/sqrt(2.);
    
    //Remove kaon and ground state pipi time dependence from data

    typedef correlationFunction<amplitudeDataCoord, jackknifeDistributionD> CorrFuncType;

    CorrFuncType ktopipi_flat = removeGroundStateDependence(ktopipi_data[q], AK, mK, Apipi0, E0);
    CorrFuncType ktosigma_flat = removeGroundStateDependence(ktosigma_data[q], AK, mK, Asigma0, E0);

    //Time dep now:
    //M0 + (A1/A0)*M1*exp(-(E1-E0) * (tsep_k_pi - t));
    
    //tsep_k_pi - t = t_op_pi
    
    //Error weighted average over data with same t_op_pi
    typedef correlationFunction<double, jackknifeDistributionD> PlotCorrFuncType;
    PlotCorrFuncType ktopipi_wavg, ktosigma_wavg;
    std::vector<PlotCorrFuncType> ktopipi_tsepkop, ktosigma_tsepkop;
    std::map<int, int> ktopipi_idx_tsep_k_op_map, ktosigma_idx_tsep_k_op_map;
  

    filterAndErrWeightedAvgData(ktopipi_wavg, ktopipi_tsepkop, ktopipi_idx_tsep_k_op_map, ktopipi_flat, tmin_k_op, Lt, stringize("K->pipi Q%d",q+1) );
    filterAndErrWeightedAvgData(ktosigma_wavg, ktosigma_tsepkop, ktosigma_idx_tsep_k_op_map, ktosigma_flat, tmin_k_op, Lt, stringize("K->sigma Q%d",q+1));

    int tsep_k_op_max = -1;
    for(auto it=ktopipi_idx_tsep_k_op_map.begin(); it != ktopipi_idx_tsep_k_op_map.end(); ++it) tsep_k_op_max = std::max(tsep_k_op_max, it->second);
    for(auto it=ktosigma_idx_tsep_k_op_map.begin(); it != ktosigma_idx_tsep_k_op_map.end(); ++it) tsep_k_op_max = std::max(tsep_k_op_max, it->second);
    
    //Construct fit curves
    PlotCorrFuncType curve_ground, curve_excited_ktopipi, curve_sum_ktopipi, curve_excited_ktosigma, curve_sum_ktosigma;

    int npoint = 60;
    double delta = double(tsep_k_op_max-tmin_op_snk)/(npoint - 1);
    for(int i=0;i<npoint;i++){
      double t_op_snk = tmin_op_snk + i*delta;
      jackknifeDistributionD gnd = M0;
      curve_ground.push_back(t_op_snk, gnd);

      jackknifeDistributionD exc = Apipi1*M1pipi*exp(-(E1-E0)*t_op_snk)/Apipi0;
      jackknifeDistributionD sum = gnd+exc;

      curve_excited_ktopipi.push_back(t_op_snk, exc);
      curve_sum_ktopipi.push_back(t_op_snk, sum);

      exc = Asigma1*M1sigma*exp(-(E1-E0)*t_op_snk)/Asigma0;
      sum = gnd+exc;

      curve_excited_ktosigma.push_back(t_op_snk, exc);
      curve_sum_ktosigma.push_back(t_op_snk, sum);
    }
  
    {
      //Plot
      MatPlotLibScriptGenerate plotter;
      typedef DataSeriesAccessor<PlotCorrFuncType, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistributionD> > accessor;
     
      //Weighted avgs
      plotter.setLegend(plotter.plotData(accessor(ktopipi_wavg),{{"color",'r'}},"ktopipi_wavg"),
			"$K\\to\\pi\\pi$ weighted avg");

      plotter.setLegend(plotter.plotData(accessor(ktosigma_wavg),{{"color",'b'}},"ktosigma_wavg"),
			"$K\\to\\sigma$ weighted avg");

      //Individual tseps
      for(int i=0;i<ktopipi_tsepkop.size();i++){
	int tsep_k_op =  ktopipi_idx_tsep_k_op_map[i];
	plotter.setLegend(plotter.plotData(accessor(ktopipi_tsepkop[i]), {{"color",'r'}}, stringize("ktopipi_tsep_k_pi%d",tsep_k_op)),
			  stringize("$K\\to\\pi\\pi$ $t_{\\rm sep}^{K\\to{\\rm op}}=%d$",tsep_k_op));
      }      
      for(int i=0;i<ktosigma_tsepkop.size();i++){
	int tsep_k_op =  ktosigma_idx_tsep_k_op_map[i];
	plotter.setLegend(plotter.plotData(accessor(ktosigma_tsepkop[i]), {{"color",'b'}}, stringize("ktosigma_tsep_k_pi%d",tsep_k_op)),
			  stringize("$K\\to\\sigma$ $t_{\\rm sep}^{K\\to{\\rm op}}=%d$",tsep_k_op));
      }      

      plotter.setLegend(plotter.errorBand(accessor(curve_ground),{{"color",'r'}, {"alpha",0.5}},"fit_gnd"),
			"Gnd");

      plotter.setLegend(plotter.errorBand(accessor(curve_excited_ktopipi),{{"color",'g'}, {"alpha",0.5}},"fit_exc_ktopipi"),
			"$K\\to\\pi\\pi$ Exc");

      plotter.setLegend(plotter.errorBand(accessor(curve_sum_ktopipi),{{"color",'b'}, {"alpha",0.5}},"fit_sum_ktopipi"),
			"$K\\to\\pi\\pi$ Tot");

      plotter.setLegend(plotter.errorBand(accessor(curve_excited_ktosigma),{{"color",'g'}, {"alpha",0.5}},"fit_exc_ktosigma"),
			"$K\\to\\sigma$ Exc");

      plotter.setLegend(plotter.errorBand(accessor(curve_sum_ktosigma),{{"color",'b'}, {"alpha",0.5}},"fit_sum_ktosigma"),
			"$K\\to\\sigma$ Tot");
      
      plotter.setXlabel("$t$");
      std::ostringstream ylabel; ylabel << "$M^{1/2,\\ \\rm{lat}}_" << q+1 << "$";
      plotter.setYlabel(ylabel.str());
      
      plotter.createLegend();
      std::ostringstream filename_stub; filename_stub << "plot_2exp_errw_Q" << q+1 << "_flat";
      plotter.write( filename_stub.str()+".py", filename_stub.str()+".pdf");
    }
  }
}




CPSFIT_END_NAMESPACE

#endif
