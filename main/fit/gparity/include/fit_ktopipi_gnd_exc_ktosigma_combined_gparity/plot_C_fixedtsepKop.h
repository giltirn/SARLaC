#ifndef _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_PLOT_C_FIXEDTSEP_K_OP
#define _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_PLOT_C_FIXEDTSEP_K_OP

#include<config.h>
#include<utils/macros.h>

#include "resampled_data.h"
#include "simfit_common.h"

SARLAC_START_NAMESPACE

//Plot C_q for operator idx 1<=q<=10 separately for all tsep(K->op) for a fixed pipi operator
//Plotting as a function of t'=tsep-t
//Output file format will be  "<PATH>/plot_C_fixed_tsepKop_q<q>_<KTOOP>.py"
template<typename DistributionType>
void plotCfixedTsepKop(const ResampledData<DistributionType> &data,
		       PiPiOperator op, const int q,
		       const std::vector<int> &tsep_k_op,
		       const std::string &path = "."){
  int ntsep_k_op = tsep_k_op.size();
  
  //The data are ordered with tsep_k_op outer and t inner
  const correlationFunction<amplitudeDataCoord, DistributionType> &qdata = data(op)[q-1];
  
  std::vector<correlationFunction<int, DistributionType> > qdata_sep(ntsep_k_op);
  
  for(int i=0;i<qdata.size();i++){
    int tsep_k_op_val = qdata.coord(i).tsep_k_pi;
    //Check if in the list of vals to plot
    int tsep_k_op_idx = -1;
    for(int tt=0;tt<ntsep_k_op;tt++){
      if(tsep_k_op_val == tsep_k_op[tt]){
	tsep_k_op_idx = tt; break;
      }
    }
    if(tsep_k_op_idx == -1) continue; //skip if not

    int t_val = qdata.coord(i).t;
    int tprime_val = tsep_k_op_val - t_val;
    qdata_sep[tsep_k_op_idx].push_back(tprime_val, qdata.value(i));
  }
  
  typedef DataSeriesAccessor<correlationFunction<int, DistributionType>, ScalarCoordinateAccessor<int>, DistributionPlotAccessor<DistributionType> > Acc;

  MatPlotLibScriptGenerate plot;
  for(int i=0;i<ntsep_k_op;i++){
    //We need to reverse the data for each dataset to get t' in the right order
    qdata_sep[i].reverse();
    std::string tag = stringize("data_tsep_k_op_%d",tsep_k_op[i]);
    auto h = plot.plotData(Acc(qdata_sep[i]), tag);
    plot.setLegend(h, stringize("%d",tsep_k_op[i]));
  }
  std::string stub = stringize("%s/plot_C_fixed_tsepKop_q%d_%s", path.c_str(), q, opDescrFile(op).c_str());
  plot.createLegend();
  plot.setXlabel("$t^{\\prime}$");
  plot.setYlabel(stringize("$C_%d$", q));
  plot.write(stub + ".py", stub + ".pdf");
}




SARLAC_END_NAMESPACE

#endif
