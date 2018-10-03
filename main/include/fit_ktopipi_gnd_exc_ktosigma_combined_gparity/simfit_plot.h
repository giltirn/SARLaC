#ifndef _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_SIMFIT_PLOT_H
#define _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_SIMFIT_PLOT_H

#include<config.h>
#include<utils/macros.h>

#include<fit_ktopipi_ktosigma_combined_gparity/simfit_plot.h>

CPSFIT_START_NAMESPACE



void plotErrorWeightedData2expFlat(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &ktopipi_data,
				   const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &ktopipi_exc_data,
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
    GETP(Apipi_exc_0);
    GETP(Apipi_exc_1);
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
    //  AK*A1*M1*exp(-E1 * tsep_k_pi)*exp( -(mK - E1)*t )/sqrt(2.)
    
    //Remove kaon and ground state pipi time dependence from data

    typedef correlationFunction<amplitudeDataCoord, jackknifeDistributionD> CorrFuncType;

    CorrFuncType ktopipi_flat = removeGroundStateDependence(ktopipi_data[q], AK, mK, Apipi0, E0);
    CorrFuncType ktopipi_exc_flat = removeGroundStateDependence(ktopipi_exc_data[q], AK, mK, Apipi_exc_0, E0);
    CorrFuncType ktosigma_flat = removeGroundStateDependence(ktosigma_data[q], AK, mK, Asigma0, E0);

    //Time dep now:
    //M0 + (A1/A0)*M1*exp(-(E1-E0) * (tsep_k_pi - t))
    
    //tsep_k_pi - t = t_op_pi
    
    //Error weighted average over data with same t_op_pi
    typedef correlationFunction<double, jackknifeDistributionD> PlotCorrFuncType;
    PlotCorrFuncType ktopipi_wavg, ktopipi_exc_wavg, ktosigma_wavg;
    std::vector<PlotCorrFuncType> ktopipi_tsepkop, ktopipi_exc_tsepkop, ktosigma_tsepkop;
    std::map<int, int> ktopipi_idx_tsep_k_op_map, ktopipi_exc_idx_tsep_k_op_map, ktosigma_idx_tsep_k_op_map;
  

    filterAndErrWeightedAvgData(ktopipi_wavg, ktopipi_tsepkop, ktopipi_idx_tsep_k_op_map, ktopipi_flat, tmin_k_op, Lt, stringize("K->pipi (111) Q%d",q+1) );
    filterAndErrWeightedAvgData(ktopipi_exc_wavg, ktopipi_exc_tsepkop, ktopipi_exc_idx_tsep_k_op_map, ktopipi_exc_flat, tmin_k_op, Lt, stringize("K->pipi (311) Q%d",q+1) );
    filterAndErrWeightedAvgData(ktosigma_wavg, ktosigma_tsepkop, ktosigma_idx_tsep_k_op_map, ktosigma_flat, tmin_k_op, Lt, stringize("K->sigma Q%d",q+1));

    int tsep_k_op_max = -1;
    for(auto it=ktopipi_idx_tsep_k_op_map.begin(); it != ktopipi_idx_tsep_k_op_map.end(); ++it) tsep_k_op_max = std::max(tsep_k_op_max, it->second);
    for(auto it=ktopipi_exc_idx_tsep_k_op_map.begin(); it != ktopipi_exc_idx_tsep_k_op_map.end(); ++it) tsep_k_op_max = std::max(tsep_k_op_max, it->second);
    for(auto it=ktosigma_idx_tsep_k_op_map.begin(); it != ktosigma_idx_tsep_k_op_map.end(); ++it) tsep_k_op_max = std::max(tsep_k_op_max, it->second);
    
    //Construct fit curves
    PlotCorrFuncType curve_ground, 
      curve_excited_ktopipi, curve_sum_ktopipi, 
      curve_excited_ktopipi_exc, curve_sum_ktopipi_exc, 
      curve_excited_ktosigma, curve_sum_ktosigma;

    int npoint = 60;
    double delta = double(tsep_k_op_max-tmin_op_snk)/(npoint - 1);
    for(int i=0;i<npoint;i++){
      double t_op_snk = tmin_op_snk + i*delta;
      jackknifeDistributionD gnd = M0;
      curve_ground.push_back(t_op_snk, gnd);

      //K->pipi gnd
      jackknifeDistributionD exc = Apipi1*M1*exp(-(E1-E0)*t_op_snk)/Apipi0;
      jackknifeDistributionD sum = gnd+exc;

      curve_excited_ktopipi.push_back(t_op_snk, exc);
      curve_sum_ktopipi.push_back(t_op_snk, sum);


      //K->pipi exc
      exc = Apipi_exc_1*M1*exp(-(E1-E0)*t_op_snk)/Apipi_exc_0;
      sum = gnd+exc;

      curve_excited_ktopipi_exc.push_back(t_op_snk, exc);
      curve_sum_ktopipi_exc.push_back(t_op_snk, sum);


      //K->sigma
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
      plotter.setLegend(plotter.plotData(accessor(ktopipi_wavg),{{"color",'r'}},"ktopipi_111_wavg"),
			"$K\\to\\pi\\pi_{111}$ weighted avg");

      plotter.setLegend(plotter.plotData(accessor(ktopipi_exc_wavg),{{"color",'g'}},"ktopipi_311_wavg"),
			"$K\\to\\pi\\pi_{311}$ weighted avg");

      plotter.setLegend(plotter.plotData(accessor(ktosigma_wavg),{{"color",'b'}},"ktosigma_wavg"),
			"$K\\to\\sigma$ weighted avg");

      //Individual tseps
      for(int i=0;i<ktopipi_tsepkop.size();i++){
	int tsep_k_op =  ktopipi_idx_tsep_k_op_map[i];
	plotter.setLegend(plotter.plotData(accessor(ktopipi_tsepkop[i]), {{"color",'r'}}, stringize("ktopipi_111_tsep_k_pi%d",tsep_k_op)),
			  stringize("$K\\to\\pi\\pi_{111}$ $t_{\\rm sep}^{K\\to{\\rm op}}=%d$",tsep_k_op));
      }      
      for(int i=0;i<ktopipi_exc_tsepkop.size();i++){
	int tsep_k_op =  ktopipi_exc_idx_tsep_k_op_map[i];
	plotter.setLegend(plotter.plotData(accessor(ktopipi_exc_tsepkop[i]), {{"color",'g'}}, stringize("ktopipi_311_tsep_k_pi%d",tsep_k_op)),
			  stringize("$K\\to\\pi\\pi_{311}$ $t_{\\rm sep}^{K\\to{\\rm op}}=%d$",tsep_k_op));      }      

      for(int i=0;i<ktosigma_tsepkop.size();i++){
	int tsep_k_op =  ktosigma_idx_tsep_k_op_map[i];
	plotter.setLegend(plotter.plotData(accessor(ktosigma_tsepkop[i]), {{"color",'b'}}, stringize("ktosigma_tsep_k_pi%d",tsep_k_op)),
			  stringize("$K\\to\\sigma$ $t_{\\rm sep}^{K\\to{\\rm op}}=%d$",tsep_k_op));
      }      

      plotter.setLegend(plotter.errorBand(accessor(curve_ground),{{"color",'r'}, {"alpha",0.5}},"fit_gnd"),
			"Gnd");
      
      //Fit curves
      //K->pipi gnd
      plotter.setLegend(plotter.errorBand(accessor(curve_excited_ktopipi),{{"color",'g'}, {"alpha",0.5}},"fit_exc_ktopipi_111"),
			"$K\\to\\pi\\pi_{111}$ Exc");

      plotter.setLegend(plotter.errorBand(accessor(curve_sum_ktopipi),{{"color",'b'}, {"alpha",0.5}},"fit_sum_ktopipi_111"),
			"$K\\to\\pi\\pi_{111}$ Tot");

      //K->pipi exc
      plotter.setLegend(plotter.errorBand(accessor(curve_excited_ktopipi_exc),{{"color",'g'}, {"alpha",0.5}},"fit_exc_ktopipi_311"),
			"$K\\to\\pi\\pi_{311}$ Exc");

      plotter.setLegend(plotter.errorBand(accessor(curve_sum_ktopipi_exc),{{"color",'b'}, {"alpha",0.5}},"fit_sum_ktopipi_311"),
			"$K\\to\\pi\\pi_{311}$ Tot");

      //K->sigma
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



void plotErrorWeightedData3expFlat(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &ktopipi_data,
				   const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &ktopipi_exc_data,
				   const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &ktosigma_data,
				   const std::vector<jackknifeDistribution<FitSimGenTwoState::Params> > &fit_params,
				   const int Lt, const int tmin_k_op, const int tmin_op_snk,
				   const FitSimGenThreeState &fitfunc){
  for(int q=0;q<10;q++){
#define GETP(NM) auto NM = peek(#NM, fit_params[q])
    GETP(AK);
    GETP(mK);
    GETP(Apipi0);
    GETP(Apipi1);
    GETP(Apipi2);
    GETP(Apipi_exc_0);
    GETP(Apipi_exc_1);
    GETP(Apipi_exc_2);
    GETP(Asigma0);
    GETP(Asigma1);
    GETP(Asigma2);
    GETP(E0);
    GETP(E1);
    GETP(E2);
    GETP(M0);
    GETP(M1);
    GETP(M2);
#undef GETP

    auto zero = mK; zeroit(zero);

    //Time dep is currently
    //AK*A0*M0*exp(-E0 * tsep_k_pi)*exp( -(mK - E0)*t )/sqrt(2.) + 
    //  AK*A1*M1*exp(-E1 * tsep_k_pi)*exp( -(mK - E1)*t )/sqrt(2.) +
    //  AK*A2*M2*exp(-E2 * tsep_k_pi)*exp( -(mK - E2)*t )/sqrt(2.)
    
    //Remove kaon and ground state pipi time dependence from data

    typedef correlationFunction<amplitudeDataCoord, jackknifeDistributionD> CorrFuncType;

    CorrFuncType ktopipi_flat = removeGroundStateDependence(ktopipi_data[q], AK, mK, Apipi0, E0);
    CorrFuncType ktopipi_exc_flat = removeGroundStateDependence(ktopipi_exc_data[q], AK, mK, Apipi_exc_0, E0);
    CorrFuncType ktosigma_flat = removeGroundStateDependence(ktosigma_data[q], AK, mK, Asigma0, E0);


    //Time dep now:
    //M0 + (A1/A0)*M1*exp(-(E1-E0) * (tsep_k_pi - t)) + (A2/A0)*M2*exp(-(E2-E0) * (tsep_k_pi - t))
    
    //tsep_k_pi - t = t_op_pi
    
    //Error weighted average over data with same t_op_pi
    typedef correlationFunction<double, jackknifeDistributionD> PlotCorrFuncType;
    PlotCorrFuncType ktopipi_wavg, ktopipi_exc_wavg, ktosigma_wavg;
    std::vector<PlotCorrFuncType> ktopipi_tsepkop, ktopipi_exc_tsepkop, ktosigma_tsepkop;
    std::map<int, int> ktopipi_idx_tsep_k_op_map, ktopipi_exc_idx_tsep_k_op_map, ktosigma_idx_tsep_k_op_map;
  

    filterAndErrWeightedAvgData(ktopipi_wavg, ktopipi_tsepkop, ktopipi_idx_tsep_k_op_map, ktopipi_flat, tmin_k_op, Lt, stringize("K->pipi (111) Q%d",q+1) );
    filterAndErrWeightedAvgData(ktopipi_exc_wavg, ktopipi_exc_tsepkop, ktopipi_exc_idx_tsep_k_op_map, ktopipi_exc_flat, tmin_k_op, Lt, stringize("K->pipi (311) Q%d",q+1) );
    filterAndErrWeightedAvgData(ktosigma_wavg, ktosigma_tsepkop, ktosigma_idx_tsep_k_op_map, ktosigma_flat, tmin_k_op, Lt, stringize("K->sigma Q%d",q+1));

    int tsep_k_op_max = -1;
    for(auto it=ktopipi_idx_tsep_k_op_map.begin(); it != ktopipi_idx_tsep_k_op_map.end(); ++it) tsep_k_op_max = std::max(tsep_k_op_max, it->second);
    for(auto it=ktopipi_exc_idx_tsep_k_op_map.begin(); it != ktopipi_exc_idx_tsep_k_op_map.end(); ++it) tsep_k_op_max = std::max(tsep_k_op_max, it->second);
    for(auto it=ktosigma_idx_tsep_k_op_map.begin(); it != ktosigma_idx_tsep_k_op_map.end(); ++it) tsep_k_op_max = std::max(tsep_k_op_max, it->second);
    
    //Construct fit curves
    PlotCorrFuncType curve_ground, 
      curve_excited_ktopipi, curve_2excited_ktopipi, curve_sum_ktopipi, 
      curve_excited_ktopipi_exc, curve_2excited_ktopipi_exc, curve_sum_ktopipi_exc, 
      curve_excited_ktosigma, curve_2excited_ktosigma, curve_sum_ktosigma;

    int npoint = 60;
    double delta = double(tsep_k_op_max-tmin_op_snk)/(npoint - 1);
    for(int i=0;i<npoint;i++){
      double t_op_snk = tmin_op_snk + i*delta;
      jackknifeDistributionD gnd = M0;
      curve_ground.push_back(t_op_snk, gnd);

      //K->pipi gnd
      jackknifeDistributionD exc = Apipi1*M1*exp(-(E1-E0)*t_op_snk)/Apipi0;
      jackknifeDistributionD exc2 = Apipi2*M2*exp(-(E2-E0)*t_op_snk)/Apipi0;
      jackknifeDistributionD sum = gnd+exc+exc2;

      curve_excited_ktopipi.push_back(t_op_snk, exc);
      curve_2excited_ktopipi.push_back(t_op_snk, exc2);
      curve_sum_ktopipi.push_back(t_op_snk, sum);


      //K->pipi exc
      exc = Apipi_exc_1*M1*exp(-(E1-E0)*t_op_snk)/Apipi_exc_0;
      exc2 = Apipi_exc_2*M2*exp(-(E2-E0)*t_op_snk)/Apipi_exc_0;
      sum = gnd+exc+exc2;

      curve_excited_ktopipi_exc.push_back(t_op_snk, exc);
      curve_2excited_ktopipi_exc.push_back(t_op_snk, exc2);
      curve_sum_ktopipi_exc.push_back(t_op_snk, sum);


      //K->sigma
      exc = Asigma1*M1*exp(-(E1-E0)*t_op_snk)/Asigma0;
      exc2 = Asigma2*M2*exp(-(E2-E0)*t_op_snk)/Asigma0;
      sum = gnd+exc+exc2;

      curve_excited_ktosigma.push_back(t_op_snk, exc);
      curve_2excited_ktosigma.push_back(t_op_snk, exc2);
      curve_sum_ktosigma.push_back(t_op_snk, sum);
    }
  
    {
      //Plot
      MatPlotLibScriptGenerate plotter;
      typedef DataSeriesAccessor<PlotCorrFuncType, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistributionD> > accessor;
     
      //Weighted avgs
      plotter.setLegend(plotter.plotData(accessor(ktopipi_wavg),{{"color",'r'}},"ktopipi_111_wavg"),
			"$K\\to\\pi\\pi_{111}$ weighted avg");

      plotter.setLegend(plotter.plotData(accessor(ktopipi_exc_wavg),{{"color",'g'}},"ktopipi_311_wavg"),
			"$K\\to\\pi\\pi_{311}$ weighted avg");

      plotter.setLegend(plotter.plotData(accessor(ktosigma_wavg),{{"color",'b'}},"ktosigma_wavg"),
			"$K\\to\\sigma$ weighted avg");

      //Individual tseps
      for(int i=0;i<ktopipi_tsepkop.size();i++){
	int tsep_k_op =  ktopipi_idx_tsep_k_op_map[i];
	plotter.setLegend(plotter.plotData(accessor(ktopipi_tsepkop[i]), {{"color",'r'}}, stringize("ktopipi_111_tsep_k_pi%d",tsep_k_op)),
			  stringize("$K\\to\\pi\\pi_{111}$ $t_{\\rm sep}^{K\\to{\\rm op}}=%d$",tsep_k_op));
      }      
      for(int i=0;i<ktopipi_exc_tsepkop.size();i++){
	int tsep_k_op =  ktopipi_exc_idx_tsep_k_op_map[i];
	plotter.setLegend(plotter.plotData(accessor(ktopipi_exc_tsepkop[i]), {{"color",'g'}}, stringize("ktopipi_311_tsep_k_pi%d",tsep_k_op)),
			  stringize("$K\\to\\pi\\pi_{311}$ $t_{\\rm sep}^{K\\to{\\rm op}}=%d$",tsep_k_op));      }      

      for(int i=0;i<ktosigma_tsepkop.size();i++){
	int tsep_k_op =  ktosigma_idx_tsep_k_op_map[i];
	plotter.setLegend(plotter.plotData(accessor(ktosigma_tsepkop[i]), {{"color",'b'}}, stringize("ktosigma_tsep_k_pi%d",tsep_k_op)),
			  stringize("$K\\to\\sigma$ $t_{\\rm sep}^{K\\to{\\rm op}}=%d$",tsep_k_op));
      }      

      plotter.setLegend(plotter.errorBand(accessor(curve_ground),{{"color",'r'}, {"alpha",0.5}},"fit_gnd"),
			"Gnd");
      
      //Fit curves
      //K->pipi gnd
      plotter.setLegend(plotter.errorBand(accessor(curve_excited_ktopipi),{{"color",'g'}, {"alpha",0.5}},"fit_exc_ktopipi_111"),
			"$K\\to\\pi\\pi_{111}$ Exc");

      plotter.setLegend(plotter.errorBand(accessor(curve_2excited_ktopipi),{{"color",'c'}, {"alpha",0.5}},"fit_exc2_ktopipi_111"),
			"$K\\to\\pi\\pi_{111}$ Exc2");

      plotter.setLegend(plotter.errorBand(accessor(curve_sum_ktopipi),{{"color",'b'}, {"alpha",0.5}},"fit_sum_ktopipi_111"),
			"$K\\to\\pi\\pi_{111}$ Tot");

      //K->pipi exc
      plotter.setLegend(plotter.errorBand(accessor(curve_excited_ktopipi_exc),{{"color",'g'}, {"alpha",0.5}},"fit_exc_ktopipi_311"),
			"$K\\to\\pi\\pi_{311}$ Exc");

      plotter.setLegend(plotter.errorBand(accessor(curve_2excited_ktopipi_exc),{{"color",'c'}, {"alpha",0.5}},"fit_exc2_ktopipi_311"),
			"$K\\to\\pi\\pi_{311}$ Exc2");

      plotter.setLegend(plotter.errorBand(accessor(curve_sum_ktopipi_exc),{{"color",'b'}, {"alpha",0.5}},"fit_sum_ktopipi_311"),
			"$K\\to\\pi\\pi_{311}$ Tot");

      //K->sigma
      plotter.setLegend(plotter.errorBand(accessor(curve_excited_ktosigma),{{"color",'g'}, {"alpha",0.5}},"fit_exc_ktosigma"),
			"$K\\to\\sigma$ Exc");

      plotter.setLegend(plotter.errorBand(accessor(curve_2excited_ktosigma),{{"color",'c'}, {"alpha",0.5}},"fit_exc2_ktosigma"),
			"$K\\to\\sigma$ Exc2");

      plotter.setLegend(plotter.errorBand(accessor(curve_sum_ktosigma),{{"color",'b'}, {"alpha",0.5}},"fit_sum_ktosigma"),
			"$K\\to\\sigma$ Tot");
      
      plotter.setXlabel("$t$");
      std::ostringstream ylabel; ylabel << "$M^{1/2,\\ \\rm{lat}}_" << q+1 << "$";
      plotter.setYlabel(ylabel.str());
      
      plotter.createLegend();
      std::ostringstream filename_stub; filename_stub << "plot_3exp_errw_Q" << q+1 << "_flat";
      plotter.write( filename_stub.str()+".py", filename_stub.str()+".pdf");
    }
  }
}

CPSFIT_END_NAMESPACE

#endif
