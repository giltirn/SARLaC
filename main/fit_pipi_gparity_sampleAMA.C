#include <fstream>
#include <algorithm>
#include <sstream>
#include <boost/timer/timer.hpp>

#include <fit.h>
#include <random.h>
#include <plot.h>
#include <distribution.h>
#include <data_series.h>
#include <parser.h>
#include <common.h>

using namespace CPSfit;

#include <fit_pipi_gparity/args.h>
#include <fit_pipi_gparity/data_containers.h>
#include <fit_pipi_gparity/mom_data_containers.h>
#include <fit_pipi_gparity/read_data.h>
#include <fit_pipi_gparity/fitfunc.h>
#include <fit_pipi_gparity/cmdline.h>
#include <fit_pipi_gparity/fit.h>
#include <fit_pipi_gparity/plot.h>
#include <fit_pipi_gparity/mom_project.h>
#include <fit_pipi_gparity/main.h>

#include <fit_simple_sampleAMA/sampleAMA_resample.h>
#include <fit_pipi_gparity_sampleAMA/args.h>
#include <fit_pipi_gparity_sampleAMA/cmdline.h>

struct rawData{ //raw, unbinned data
  bubbleDataAllMomenta raw_bubble_data_sloppy_S;
  rawCorrelationFunction pipi_raw_sloppy_S;

  bubbleDataAllMomenta raw_bubble_data_sloppy_C;
  rawCorrelationFunction pipi_raw_sloppy_C;

  bubbleDataAllMomenta raw_bubble_data_exact_C;
  rawCorrelationFunction pipi_raw_exact_C;

  int nS; //unbinned
  int nC;

  rawData(const PiPiProject &proj_src, const PiPiProject &proj_snk, const PiPiMomAllow &allow, const int isospin, const ArgsSampleAMA &args, const CMDlineSampleAMA &cmdline){
    bubbleDataAllMomenta* raw_bubble_data[3] = { &raw_bubble_data_sloppy_S, &raw_bubble_data_sloppy_C, &raw_bubble_data_exact_C };
    rawCorrelationFunction* pipi_raw[3] = { &pipi_raw_sloppy_S, &pipi_raw_sloppy_C, &pipi_raw_exact_C };
    const char ens[3] = { 'S', 'C', 'C' };
    const SloppyExact se[3] = { Sloppy, Sloppy, Exact };

    for(int i=0;i<3;i++){
      Args argsi = args.toArgs(ens[i]);
      CMDline cmdlinei = cmdline.toCMDline(se[i], ens[i]);
      figureDataAllMomenta raw_data;
      readRawData(raw_data, *raw_bubble_data[i], argsi, cmdlinei);
      getRawPiPiCorrFunc(*pipi_raw[i], raw_data, *raw_bubble_data[i], proj_src, proj_snk, allow, isospin, argsi, cmdlinei, "", false);
    }

    nS = pipi_raw_sloppy_S.value(0).size();
    nC = pipi_raw_sloppy_C.value(0).size();
  }
};
    
template<typename DistributionType>
bubbleDataAllMomentaBase<bubbleDataBase<DistributionType> > resampleCorrectBubbleSampleAMA(const rawData &raw, const int nS, const int nC, const int bin_size){
  const int Lt = raw.raw_bubble_data_sloppy_S.getLt();

  const int nSb = nS/bin_size;
  const int nCb = nC/bin_size;

  bubbleDataAllMomentaBase<bubbleDataBase<DistributionType> > out(Lt,nSb+nCb);

  std::vector<threeMomentum> momenta;
  for(auto it = raw.raw_bubble_data_sloppy_S.begin(); it != raw.raw_bubble_data_sloppy_S.end(); it++) momenta.push_back(it->first);
  
  for(int i=0;i<momenta.size();i++)
    for(int t=0;t<Lt;t++){
      const threeMomentum &mom = momenta[i];
      DistributionType sloppy_S = sampleAMAresample<DistributionType>::resample(raw.raw_bubble_data_sloppy_S(mom)(t).bin(bin_size),'S',nSb,nCb);
      DistributionType sloppy_C = sampleAMAresample<DistributionType>::resample(raw.raw_bubble_data_sloppy_C(mom)(t).bin(bin_size),'C',nSb,nCb);
      DistributionType exact_C = sampleAMAresample<DistributionType>::resample(raw.raw_bubble_data_exact_C(mom)(t).bin(bin_size),'C',nSb,nCb);
      out(mom)(t) = sloppy_S + exact_C - sloppy_C;
    }
  return out;
}

template<typename DistributionType>
void resampleCombineData(correlationFunction<double,DistributionType> &pipi, 
			 const PiPiProject &proj_src, const PiPiProject &proj_snk, const PiPiMomAllow &allow, const int isospin,
			 const rawData &raw, const ArgsSampleAMA &args, const CMDlineSampleAMA &cmdline){
  typedef correlationFunction<double,DistributionType> CorrFunc;
  typedef typename CorrFunc::ElementType Elem;
  pipi.resize(args.Lt);
  const int nS = raw.nS;
  const int nC = raw.nC;

  const int nSb = nS/args.bin_size;
  const int nCb = nC/args.bin_size;

  DistributionType pipi_sloppy_S_r, pipi_sloppy_C_r, pipi_exact_C_r;
  for(int t=0;t<args.Lt;t++){
    pipi.coord(t) = t;

    pipi_sloppy_S_r = sampleAMAresample<DistributionType>::resample(raw.pipi_raw_sloppy_S.value(t).bin(args.bin_size),'S',nSb,nCb);
    pipi_sloppy_C_r = sampleAMAresample<DistributionType>::resample(raw.pipi_raw_sloppy_C.value(t).bin(args.bin_size),'C',nSb,nCb);
    pipi_exact_C_r = sampleAMAresample<DistributionType>::resample(raw.pipi_raw_exact_C.value(t).bin(args.bin_size),'C',nSb,nCb);

    pipi.value(t) = pipi_sloppy_S_r + pipi_exact_C_r - pipi_sloppy_C_r;
  }
  if(isospin == 0 && args.do_vacuum_subtraction){
    auto bubble_data_corrected_r = resampleCorrectBubbleSampleAMA<DistributionType>(raw,nS,nC,args.bin_size);
    CorrFunc A2_realavg_V_r = computeVprojectSourceAvg(bubble_data_corrected_r,args.tsep_pipi, proj_src, proj_snk, allow);
    pipi = pipi - 3*A2_realavg_V_r;
  }
}

void generateData(jackknifeCorrelationFunction &pipi_j, doubleJackCorrelationFunction &pipi_dj, 
		  const PiPiProject &proj_src, const PiPiProject &proj_snk, const PiPiMomAllow &allow, const int isospin,
		  const ArgsSampleAMA &args, const CMDlineSampleAMA &cmdline){
  rawData raw(proj_src, proj_snk, allow, isospin, args,cmdline);
  resampleCombineData<jackknifeDistributionD>(pipi_j,proj_src, proj_snk, allow, isospin, raw,args,cmdline);
  resampleCombineData<doubleJackknifeDistributionD>(pipi_dj,proj_src, proj_snk, allow, isospin, raw,args,cmdline);
  
  pipi_j = fold(pipi_j, args.tsep_pipi);
  pipi_dj = fold(pipi_dj, args.tsep_pipi);
}

void getData(const PiPiProject &proj_src, const PiPiProject &proj_snk, const PiPiMomAllow &allow, const int isospin, 
	     jackknifeCorrelationFunction &pipi_j, doubleJackCorrelationFunction &pipi_dj, const ArgsSampleAMA &args, const CMDlineSampleAMA &cmdline){
  if(cmdline.load_combined_data){
#ifdef HAVE_HDF5
    HDF5reader reader(cmdline.load_combined_data_file);
    read(reader, pipi_j, "pipi_j");
    read(reader, pipi_dj, "pipi_dj");
#else
    error_exit("getData: Reading amplitude data requires HDF5\n");
#endif
  }else{
    generateData(pipi_j,pipi_dj, proj_src, proj_snk, allow, isospin,args,cmdline);
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

  PiPiProject *proj_src = getProjector(args.proj_src);
  PiPiProject *proj_snk = getProjector(args.proj_snk);
  PiPiMomAllow* allow = getMomPairFilter(args.allowed_mom);

  //Read resampled
  jackknifeCorrelationFunction pipi_j;
  doubleJackCorrelationFunction pipi_dj;

  getData(*proj_src, *proj_snk, *allow, args.isospin, pipi_j,pipi_dj,args,cmdline);
  
  delete proj_src; delete proj_snk; delete allow;
  
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
