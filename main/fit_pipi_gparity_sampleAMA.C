#include <fstream>
#include <algorithm>
#include <minimizer.h>
#include <cost_function.h>
#include <fitfunc.h>
#include <fitfunc_mapping.h>
#include <random.h>
#include <plot.h>
#include <distribution.h>
#include <data_series.h>
#include <parser.h>
#include <common_defs.h>
#include <fit_wrapper.h>
#include <sstream>
#include <boost/timer/timer.hpp>

#include <fit_pipi_gparity/args.h>
#include <fit_pipi_gparity/data_containers.h>
#include <fit_pipi_gparity/mom_data_containers.h>
#include <fit_pipi_gparity/read_data.h>
#include <fit_pipi_gparity/fitfunc.h>
#include <fit_pipi_gparity/cmdline.h>
#include <fit_pipi_gparity/fit.h>
#include <fit_pipi_gparity/plot.h>
#include <fit_pipi_gparity/main.h>

#include <fit_simple_sampleAMA/sampleAMA_resample.h>
#include <fit_pipi_gparity_sampleAMA/args.h>
#include <fit_pipi_gparity_sampleAMA/cmdline.h>

bubbleDataDoubleJackAllMomenta doubleJackknifeResampleBubbleSampleAMA(const bubbleDataAllMomenta &bubbles, const char ens, const int nS, const int nC){
  int Lt = bubbles.getLt();
  bubbleDataDoubleJackAllMomenta out(Lt,nS+nC);

  for(auto it = bubbles.begin(); it != bubbles.end(); it++){
    const threeMomentum & mom = it->first;
    const typename bubbleDataAllMomenta::ContainerType & raw = it->second;
    for(int t=0;t<Lt;t++)
      out(mom)(t) = ens == 'S' ? superDoubleJackknifeResampleS(raw(t),nS,nC) : superDoubleJackknifeResampleC(raw(t),nS,nC);
  }
  return out;
}
typedef bubbleDataBase<jackknifeDistributionD > bubbleDataJack;
typedef bubbleDataAllMomentaBase<bubbleDataJack> bubbleDataJackAllMomenta;

bubbleDataJackAllMomenta jackknifeResampleBubbleSampleAMA(const bubbleDataAllMomenta &bubbles, const char ens, const int nS, const int nC){
  int Lt = bubbles.getLt();
  bubbleDataJackAllMomenta out(Lt,nS+nC);

  for(auto it = bubbles.begin(); it != bubbles.end(); it++){
    const threeMomentum & mom = it->first;
    const typename bubbleDataAllMomenta::ContainerType & raw = it->second;
    for(int t=0;t<Lt;t++)
      out(mom)(t) = ens == 'S' ? superJackknifeResampleS(raw(t),nS,nC) : superJackknifeResampleC(raw(t),nS,nC);
  }
  return out;
}


void generateData(jackknifeCorrelationFunction &pipi_j, doubleJackCorrelationFunction &pipi_dj, const ArgsSampleAMA &args, const CMDlineSampleAMA &cmdline){
  bubbleDataAllMomenta raw_bubble_data_sloppy_S;
  rawCorrelationFunction pipi_raw_sloppy_S;
  getRawData(pipi_raw_sloppy_S, raw_bubble_data_sloppy_S, args.toArgs('S'), cmdline.toCMDline(Sloppy,'S') );

  bubbleDataAllMomenta raw_bubble_data_sloppy_C;
  rawCorrelationFunction pipi_raw_sloppy_C;
  getRawData(pipi_raw_sloppy_C, raw_bubble_data_sloppy_C, args.toArgs('C'), cmdline.toCMDline(Sloppy,'C') );

  bubbleDataAllMomenta raw_bubble_data_exact_C;
  rawCorrelationFunction pipi_raw_exact_C;
  getRawData(pipi_raw_exact_C, raw_bubble_data_exact_C, args.toArgs('C'), cmdline.toCMDline(Exact,'C') );
  
  const int nS = pipi_raw_sloppy_S.value(0).size();
  const int nC = pipi_raw_sloppy_C.value(0).size();
  
  typedef jackknifeCorrelationFunction::ElementType JackCorrElem;
  jackknifeCorrelationFunction pipi_sloppy_S_j(args.Lt, [&](const int t){ return JackCorrElem(pipi_raw_sloppy_S.coord(t), superJackknifeResampleS(pipi_raw_sloppy_S.value(t),nS,nC)); });
  jackknifeCorrelationFunction pipi_sloppy_C_j(args.Lt, [&](const int t){ return JackCorrElem(pipi_raw_sloppy_C.coord(t), superJackknifeResampleC(pipi_raw_sloppy_C.value(t),nS,nC)); });
  jackknifeCorrelationFunction pipi_exact_C_j(args.Lt, [&](const int t){ return JackCorrElem(pipi_raw_exact_C.coord(t), superJackknifeResampleC(pipi_raw_exact_C.value(t),nS,nC)); });


  typedef doubleJackCorrelationFunction::ElementType DJackCorrElem;
  doubleJackCorrelationFunction pipi_sloppy_S_dj(args.Lt, [&](const int t){ return DJackCorrElem(pipi_raw_sloppy_S.coord(t), superDoubleJackknifeResampleS(pipi_raw_sloppy_S.value(t),nS,nC)); });
  doubleJackCorrelationFunction pipi_sloppy_C_dj(args.Lt, [&](const int t){ return DJackCorrElem(pipi_raw_sloppy_C.coord(t), superDoubleJackknifeResampleC(pipi_raw_sloppy_C.value(t),nS,nC)); });
  doubleJackCorrelationFunction pipi_exact_C_dj(args.Lt, [&](const int t){ return DJackCorrElem(pipi_raw_exact_C.coord(t), superDoubleJackknifeResampleC(pipi_raw_exact_C.value(t),nS,nC)); });
  
  if(args.do_vacuum_subtraction){
    //Jackknife
    bubbleDataJackAllMomenta j_bubble_data_sloppy_S = jackknifeResampleBubbleSampleAMA(raw_bubble_data_sloppy_S, 'S', nS, nC);
    bubbleDataJackAllMomenta j_bubble_data_sloppy_C = jackknifeResampleBubbleSampleAMA(raw_bubble_data_sloppy_C, 'C', nS, nC);
    bubbleDataJackAllMomenta j_bubble_data_exact_C = jackknifeResampleBubbleSampleAMA(raw_bubble_data_exact_C, 'C', nS, nC);

    jackknifeCorrelationFunction A2_realavg_V_sloppy_S_j = computeVprojectA2sourceAvg(j_bubble_data_sloppy_S,args.tsep_pipi);
    jackknifeCorrelationFunction A2_realavg_V_sloppy_C_j = computeVprojectA2sourceAvg(j_bubble_data_sloppy_C,args.tsep_pipi);
    jackknifeCorrelationFunction A2_realavg_V_exact_C_j = computeVprojectA2sourceAvg(j_bubble_data_exact_C,args.tsep_pipi);

    pipi_sloppy_S_j = pipi_sloppy_S_j - 3*A2_realavg_V_sloppy_S_j;
    pipi_sloppy_C_j = pipi_sloppy_C_j - 3*A2_realavg_V_sloppy_C_j;
    pipi_exact_C_j = pipi_exact_C_j - 3*A2_realavg_V_exact_C_j;

    //Double-jackknife
    bubbleDataDoubleJackAllMomenta dj_bubble_data_sloppy_S = doubleJackknifeResampleBubbleSampleAMA(raw_bubble_data_sloppy_S, 'S', nS, nC);
    bubbleDataDoubleJackAllMomenta dj_bubble_data_sloppy_C = doubleJackknifeResampleBubbleSampleAMA(raw_bubble_data_sloppy_C, 'C', nS, nC);
    bubbleDataDoubleJackAllMomenta dj_bubble_data_exact_C = doubleJackknifeResampleBubbleSampleAMA(raw_bubble_data_exact_C, 'C', nS, nC);

    doubleJackCorrelationFunction A2_realavg_V_sloppy_S_dj = computeVprojectA2sourceAvg(dj_bubble_data_sloppy_S,args.tsep_pipi);
    doubleJackCorrelationFunction A2_realavg_V_sloppy_C_dj = computeVprojectA2sourceAvg(dj_bubble_data_sloppy_C,args.tsep_pipi);
    doubleJackCorrelationFunction A2_realavg_V_exact_C_dj = computeVprojectA2sourceAvg(dj_bubble_data_exact_C,args.tsep_pipi);

    pipi_sloppy_S_dj = pipi_sloppy_S_dj - 3*A2_realavg_V_sloppy_S_dj;
    pipi_sloppy_C_dj = pipi_sloppy_C_dj - 3*A2_realavg_V_sloppy_C_dj;
    pipi_exact_C_dj = pipi_exact_C_dj - 3*A2_realavg_V_exact_C_dj;
  }

  pipi_dj = pipi_sloppy_S_dj + pipi_exact_C_dj - pipi_sloppy_C_dj;
  pipi_j = pipi_sloppy_S_j + pipi_exact_C_j - pipi_sloppy_C_j;
  
  pipi_j = fold(pipi_j, args.tsep_pipi);
  pipi_dj = fold(pipi_dj, args.tsep_pipi);
}

void getData(jackknifeCorrelationFunction &pipi_j, doubleJackCorrelationFunction &pipi_dj, const ArgsSampleAMA &args, const CMDlineSampleAMA &cmdline){
  if(cmdline.load_combined_data){
#ifdef HAVE_HDF5
    HDF5reader reader(cmdline.load_combined_data_file);
    read(reader, pipi_j, "pipi_j");
    read(reader, pipi_dj, "pipi_dj");
#else
    error_exit("getData: Reading amplitude data requires HDF5\n");
#endif
  }else{
    generateData(pipi_j,pipi_dj,args,cmdline);
  }

  if(cmdline.save_combined_data){
#ifdef HAVE_HDF5
    HDF5writer writer(cmdline.save_combined_data_file);
    write(writer, pipi_j, "pipi_j");
    write(writer, pipi_dj, "pipi_dj");
#else
    error_exit("getData: Saving amplitude data requires HDF5\n");
#endif
  }
}

int main(const int argc, const char* argv[]){
  ArgsSampleAMA args;
  if(argc < 2){
    std::ofstream of("template.args");
    of << args;
    std::cout << "Wrote template argument file to template.args\n";
    return 0;
  }
  {
    const std::string arg_file = argv[1];
    parse(args,arg_file);
    std::cout << "Read arguments: \n" << args << std::endl;
  }
  
  CMDlineSampleAMA cmdline(argc,argv,2);

  //Read resampled
  jackknifeCorrelationFunction pipi_j;
  doubleJackCorrelationFunction pipi_dj;

  getData(pipi_j,pipi_dj,args,cmdline);
  
  
  //Filter out the data that is to be fitted
  const int nsample = pipi_j.value(0).size();
  
  filterXrange<double> trange(args.t_min,args.t_max);
  
  doubleJackCorrelationFunction pipi_dj_inrange;
  jackknifeCorrelationFunction pipi_j_inrange;
  for(int d=0;d<pipi_dj.size();d++)
    if(trange.accept(pipi_dj.coord(d),pipi_dj.value(d) )){
      pipi_dj_inrange.push_back(pipi_dj[d]);
      pipi_j_inrange.push_back(pipi_j[d]);
    }
  
  //Perform the fit
  Args args_gen = args.toArgs('S'); //traj info and data dirs no longer used
  CMDline cmdline_gen = cmdline.toCMDline(Sloppy,'S');

  std::pair<jackknifeDistributionD,jackknifeDistributionD> Epipi_and_const = fit(pipi_j_inrange,pipi_dj_inrange,args_gen,cmdline_gen);

  plot(pipi_j,Epipi_and_const.first,Epipi_and_const.second,args_gen,cmdline_gen);
  
  std::cout << "Done\n";
  return 0;
}
