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

//Fit params must contain tagged parameters:  E* M*  Apipi* Apipi_exc_*  Asigma* (* \in (0..nstate-1) ). The last 3 are optional depending on operators used


template< template<typename, template<typename> class> class DistributionType > 
void plotErrorWeightedDataNexpFlat(const ResampledData<DistributionType<double, basic_vector> > &data_j,
				   const std::vector<PiPiOperator> &operators,
				   const std::vector<DistributionType<taggedValueContainer<double,std::string>, basic_vector> > &fit_params,
				   const DistributionType<double, basic_vector> &mK,
				   const DistributionType<double, basic_vector> &AK,
				   const int nstate,
				   const int Lt, const int tmin_k_op, const int tmin_op_snk){
  std::cout << "Generating plots of error weighted data versus fit" << std::endl;
  typedef DistributionType<double, basic_vector> DistributionTypeD;

  const int nop = operators.size();
  const int nQ = fit_params.size();

  typedef taggedValueContainer<double,std::string> Params;
  typedef correlationFunction<amplitudeDataCoord, DistributionTypeD> CorrFunc;

  static const std::vector<PiPiOperator> all_ops = {PiPiOperator::PiPiGnd, PiPiOperator::PiPiExc, PiPiOperator::Sigma};
  static const std::vector<std::string> all_ops_pfmt = {"Apipi%d", "Apipi_exc_%d", "Asigma%d"};
  static const std::vector<std::string> all_ops_descr = {"K->pipi(111)", "K->pipi(311)", "K->sigma"};
  static const std::vector<std::string> all_ops_descr_short = {"kpipi_111", "kpipi_311", "ksigma"};
  static const std::vector<std::string> all_ops_latex = { "$K\\to\\pi\\pi_{111}$", "$K\\to\\pi\\pi_{311}$", "$K\\to\\sigma$" };

  //Which of the above are in our operator list
  std::vector<int> op_idx(nop);
  for(int o=0;o<nop;o++){
    int opidx=-1; for(int a=0;a<all_ops.size();a++) if(all_ops[a] == operators[o]){ opidx=a; break; }
    assert(opidx != -1);
    op_idx[o] = opidx;
  }

  for(int q=0;q<nQ;q++){
    std::cout << "Starting plot for Q" << q+1 << std::endl;

    //Get the fit parameters
    std::vector<std::vector<DistributionTypeD> > A(nop, std::vector<DistributionTypeD>(nstate));
    for(int o=0;o<nop;o++){
      int opidx=op_idx[o];
      for(int s=0;s<nstate;s++){
	//Get the amplitude parameter tag for this op and state
	std::string pnm = stringize(all_ops_pfmt[opidx].c_str(),s);	
	A[o][s] = peek(pnm, fit_params[q]);
      }
    }
    std::vector<DistributionTypeD> E(nstate);
    std::vector<DistributionTypeD> M(nstate);
    
    for(int s=0;s<nstate;s++){
      E[s] = peek(stringize("E%d",s),fit_params[q]);
      M[s] = peek(stringize("M%d",s),fit_params[q]);
    }
      
    //Get the data for this q
    std::vector<CorrFunc> qdata(nop);
    for(int o=0;o<nop;o++) qdata[o] = data_j(operators[o])[q];

    //Divide out all the kaon and sink operator dependence associated with the ground state
    for(int o=0;o<nop;o++)
      qdata[o] = removeGroundStateDependence(qdata[o], AK, mK, A[o][0], E[0]);

    //Time dep now:
    //M0 + (A1/A0)*M1*exp(-(E1-E0) * (tsep_k_pi - t)) + (A2/A0)*M2*exp(-(E2-E0) * (tsep_k_pi - t)) + ....  
    //tsep_k_pi - t = t_op_pi
    
    //Generate error weighted average over data with same t_op_pi and also split data out over different tsep_k_op  (op here refers to the sink operator)
    typedef correlationFunction<double, DistributionTypeD> PlotCorrFuncType;
    std::vector<PlotCorrFuncType> qdata_wavg(nop);
    std::vector<std::vector<PlotCorrFuncType> > qdata_tsepkop(nop); //[op][tsep_k_pi_idx] {data}   note tsep_k_pi_idx is an index for the available tsep_k_pi (mapping below)
    std::vector<std::map<int, int> > op_idx_tsepkop_map(nop); //for each op, a map of tsep_k_op index to tsep_k_op

    int tsep_k_op_max = -1;
    for(int o=0;o<nop;o++){ 
      filterAndErrWeightedAvgData(qdata_wavg[o], qdata_tsepkop[o], op_idx_tsepkop_map[o], qdata[o], tmin_k_op, Lt, stringize("%s Q%d", all_ops_descr[op_idx[o]].c_str(),q+1) );
      for(auto it=op_idx_tsepkop_map[o].begin(); it != op_idx_tsepkop_map[o].end(); ++it) tsep_k_op_max = std::max(tsep_k_op_max, it->second);
    }
    
    //Construct fit curves
    PlotCorrFuncType curve_ground;
    std::vector<std::vector<PlotCorrFuncType> > curve_exc(nop, std::vector<PlotCorrFuncType>(nstate-1));
    std::vector<PlotCorrFuncType> curve_sum(nop);

    int npoint = 60;
    double delta = double(tsep_k_op_max-tmin_op_snk)/(npoint - 1);
    for(int i=0;i<npoint;i++){
      double t_op_snk = tmin_op_snk + i*delta;
      DistributionTypeD gnd = M[0];
      curve_ground.push_back(t_op_snk, gnd);

      for(int o=0;o<nop;o++){
	DistributionTypeD sum = gnd;
	for(int s=1;s<nstate;s++){
	  DistributionTypeD term = A[o][s]*M[s]*exp(-(E[s]-E[0])*t_op_snk)/A[o][0];
	  curve_exc[o][s-1].push_back(t_op_snk, term);
	  sum = sum + term;
	}
	curve_sum[o].push_back(t_op_snk, sum);
      }
    }
  
    {
      //Plot
      MatPlotLibScriptGenerate plotter;
      typedef DataSeriesAccessor<PlotCorrFuncType, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<DistributionTypeD> > accessor;

      static const char colors[6] = { 'r','g','b','c','m','k' };
      assert(nop <= 6);
     
      std::vector<std::string> shortd(nop), latex(nop);
      for(int o=0;o<nop;o++){ shortd[o] = all_ops_descr_short[op_idx[o]];  latex[o] = all_ops_latex[op_idx[o]]; }

      //Weighted avgs
      for(int o=0;o<nop;o++){
	auto handle = plotter.plotData(accessor(qdata_wavg[o]),  {{"color",colors[o]}}, stringize("%s_wavg", shortd[o].c_str()) );
	plotter.setLegend(handle, stringize("%s weighted avg", latex[o].c_str()) );
      }

      //Individual tseps
      for(int o=0;o<nop;o++){
	for(int i=0;i<qdata_tsepkop[o].size();i++){
	  int tsep_k_op = op_idx_tsepkop_map[o][i];
	  auto handle = plotter.plotData(accessor(qdata_tsepkop[o][i]), {{"color",colors[o]}}, stringize("%s_tsep_k_pi%d", shortd[o].c_str(),tsep_k_op));
	  plotter.setLegend(handle, stringize("%s $t_{\\rm sep}^{K\\to{\\rm op}}=%d$", latex[o].c_str(),tsep_k_op));
	}
      }      

      //Ground constant
      plotter.setLegend(plotter.errorBand(accessor(curve_ground),{{"color",colors[0]}, {"alpha",0.5}},"fit_gnd"),
			"Gnd");
      
      //Fit curves
      for(int o=0;o<nop;o++){
	//Excited state contributions
	for(int es=0;es<nstate-1;es++){ //excited state index
	  auto handle = plotter.errorBand(accessor(curve_exc[o][es]),{{"color",colors[o]}, {"alpha",0.5}}, stringize("fit_exc%d_%s", es+1, shortd[o].c_str() ) );
	  plotter.setLegend(handle, stringize("%s Exc %d", latex[o].c_str(), es+1) );
	}
	//Total
	auto handle = plotter.errorBand(accessor(curve_sum[o]),{{"color",colors[o]}, {"alpha",0.5}},stringize("fit_sum_%s", shortd[o].c_str()) );
	plotter.setLegend(handle, stringize("%s Tot", latex[o].c_str()) );
      }

      plotter.setXlabel("$t$");
      std::ostringstream ylabel; ylabel << "$M^{1/2,\\ \\rm{lat}}_" << q+1 << "$";
      plotter.setYlabel(ylabel.str());
      
      plotter.createLegend();
      std::ostringstream filename_stub; filename_stub << "plot_errw_Q" << q+1 << "_flat";
      plotter.write( filename_stub.str()+".py", filename_stub.str()+".pdf");

      std::cout << "Finished plotting for Q" << q+1 << std::endl;
    }
  }
  std::cout << "Finished plotting" << std::endl;
}



CPSFIT_END_NAMESPACE

#endif
