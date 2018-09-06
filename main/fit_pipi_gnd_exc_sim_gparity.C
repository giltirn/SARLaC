#include<pipi_common/simfit_generic.h>

#include <fit_pipi_sigma_sim_gparity/read_data.h>
#include <fit_pipi_sigma_sim_gparity/resampled_correlator.h>
using namespace CPSfit;

#include<fit_pipi_gnd_exc_sim_gparity/args.h>
#include<fit_pipi_gnd_exc_sim_gparity/cmdline.h>

struct readRawDataOptions{
  bool load_hdf5_data_checkpoint;
  std::string load_hdf5_data_checkpoint_stub;
  bool save_hdf5_data_checkpoint;
  std::string save_hdf5_data_checkpoint_stub;
  readRawDataOptions(): load_hdf5_data_checkpoint(false), save_hdf5_data_checkpoint(false){}
};

void readRawData(bubbleDataAllMomenta &raw_bubble_gnd_gnd, bubbleDataAllMomenta &raw_bubble_exc_exc, bubbleDataAllMomenta &raw_bubble_gnd_exc,
		 rawCorrelationFunction &raw_data_gnd_gnd, rawCorrelationFunction &raw_data_exc_exc, rawCorrelationFunction &raw_data_gnd_exc,
		 const std::string &data_dir, const std::string &figure_file_format, const std::string &bubble_file_format,
		 const int Lt, const int tsep_pipi, const int tstep_pipi, const std::vector<threeMomentum> &pion_mom,
		 const int traj_start, const int traj_inc, const int traj_lessthan,
		 const readRawDataOptions &opt = readRawDataOptions()){
  
  if(opt.load_hdf5_data_checkpoint){
    HDF5reader rd(opt.load_hdf5_data_checkpoint_stub);
    read(rd, raw_bubble_gnd_gnd, "raw_bubble_gnd_gnd");
    read(rd, raw_bubble_exc_exc, "raw_bubble_exc_exc");
    read(rd, raw_bubble_gnd_exc, "raw_bubble_gnd_exc");
    read(rd, raw_data_gnd_gnd, "raw_data_gnd_gnd");
    read(rd, raw_data_exc_exc, "raw_data_exc_exc");
    read(rd, raw_data_gnd_exc, "raw_data_gnd_exc");
  }else{
    figureData::useFileCache() = true;
    readPiPi2pt(raw_data_gnd_gnd, raw_bubble_gnd_gnd, data_dir, figure_file_format, bubble_file_format, tsep_pipi, tstep_pipi, Lt, traj_start, traj_inc, traj_lessthan, 
		pion_mom, PiPiProjector::A1momSet111, PiPiProjector::A1momSet111);
  
    figureData::getFileCache().clear();

    readPiPi2pt(raw_data_exc_exc, raw_bubble_exc_exc, data_dir, figure_file_format, bubble_file_format, tsep_pipi, tstep_pipi, Lt, traj_start, traj_inc, traj_lessthan, 
		pion_mom, PiPiProjector::A1momSet311, PiPiProjector::A1momSet311);

    figureData::getFileCache().clear();

    readPiPi2pt(raw_data_gnd_exc, raw_bubble_gnd_exc, data_dir, figure_file_format, bubble_file_format, tsep_pipi, tstep_pipi, Lt, traj_start, traj_inc, traj_lessthan, 
		pion_mom, PiPiProjector::A1momSet111, PiPiProjector::A1momSet311);

    figureData::getFileCache().clear();
  }

  if(opt.save_hdf5_data_checkpoint){
    HDF5writer wr(opt.save_hdf5_data_checkpoint_stub);
    write(wr, raw_bubble_gnd_gnd, "raw_bubble_gnd_gnd");
    write(wr, raw_bubble_exc_exc, "raw_bubble_exc_exc");
    write(wr, raw_bubble_gnd_exc, "raw_bubble_gnd_exc");
    write(wr, raw_data_gnd_gnd, "raw_data_gnd_gnd");
    write(wr, raw_data_exc_exc, "raw_data_exc_exc");
    write(wr, raw_data_gnd_exc, "raw_data_gnd_exc");
  }
}

struct generateResampledDataOptions{
  bool load_combined_data;
  std::string load_combined_data_file;

  bool save_combined_data;
  std::string save_combined_data_file;

  generateResampledDataOptions(): load_combined_data(false), save_combined_data(false){}
};

void generateResampledData(doubleJackCorrelationFunction &dj_data_gnd_gnd, doubleJackCorrelationFunction &dj_data_exc_exc, doubleJackCorrelationFunction &dj_data_gnd_exc,
			   const bubbleDataAllMomenta &raw_bubble_gnd_gnd, const bubbleDataAllMomenta &raw_bubble_exc_exc, const bubbleDataAllMomenta &raw_bubble_gnd_exc,
			   const rawCorrelationFunction &raw_data_gnd_gnd, const rawCorrelationFunction &raw_data_exc_exc, const rawCorrelationFunction &raw_data_gnd_exc,
			   const int tsep_pipi, const std::vector<threeMomentum> &pion_mom, const int bin_size, 
			   const bool do_vacuum_subtraction, const generateResampledDataOptions &opt = generateResampledDataOptions()){

  if(opt.load_combined_data){
    HDF5reader rd(opt.load_combined_data_file);
    read(rd, dj_data_gnd_gnd, "dj_data_gnd_gnd");
    read(rd, dj_data_exc_exc, "dj_data_exc_exc");
    read(rd, dj_data_gnd_exc, "dj_data_gnd_exc");
  }else{
    dj_data_gnd_gnd = binDoubleJackResample(raw_data_gnd_gnd, bin_size);
    dj_data_exc_exc = binDoubleJackResample(raw_data_exc_exc, bin_size);
    dj_data_gnd_exc = binDoubleJackResample(raw_data_gnd_exc, bin_size);
  
    //Compute vacuum subtractions
    if(do_vacuum_subtraction){
      std::cout << "Computing vacuum subtractions" << std::endl;
      doubleJackCorrelationFunction vac_sub_dj = computePiPi2ptVacSub(raw_bubble_gnd_gnd, bin_size, tsep_pipi, pion_mom, PiPiProjector::A1momSet111, PiPiProjector::A1momSet111);
      dj_data_gnd_gnd = dj_data_gnd_gnd - vac_sub_dj;

      vac_sub_dj = computePiPi2ptVacSub(raw_bubble_exc_exc, bin_size, tsep_pipi, pion_mom, PiPiProjector::A1momSet311, PiPiProjector::A1momSet311);
      dj_data_exc_exc = dj_data_exc_exc - vac_sub_dj;

      vac_sub_dj = computePiPi2ptVacSub(raw_bubble_gnd_exc, bin_size, tsep_pipi, pion_mom, PiPiProjector::A1momSet111, PiPiProjector::A1momSet311);
      dj_data_gnd_exc = dj_data_gnd_exc - vac_sub_dj;
    }

    std::cout << "Folding data" << std::endl;
    //Fold data
    dj_data_gnd_gnd = foldPiPi2pt(dj_data_gnd_gnd, tsep_pipi);
    dj_data_exc_exc = foldPiPi2pt(dj_data_exc_exc, tsep_pipi);
    dj_data_gnd_exc = foldPiPi2pt(dj_data_gnd_exc, tsep_pipi);
  }

  if(opt.save_combined_data){
    HDF5writer wr(opt.save_combined_data_file);
    write(wr, dj_data_gnd_gnd, "dj_data_gnd_gnd");
    write(wr, dj_data_exc_exc, "dj_data_exc_exc");
    write(wr, dj_data_gnd_exc, "dj_data_gnd_exc");
  }

}

template<typename FitFunc>
void fit_ff(jackknifeDistribution<typename FitFunc::Params> &params, jackknifeDistributionD &chisq, jackknifeDistributionD &chisq_per_dof,
	    const correlationFunction<SimFitCoordGen, jackknifeDistributionD> &corr_comb_j,
	    const correlationFunction<SimFitCoordGen, doubleJackknifeDistributionD> &corr_comb_dj,
	    const FitFunc &fitfunc){
    typedef typename composeFitPolicy<FitFunc, standardFitFuncPolicy, correlatedFitPolicy>::type FitPolicies;
  
    fitter<FitPolicies> fit;
    fit.importFitFunc(fitfunc);
    
    importCostFunctionParameters<correlatedFitPolicy, FitPolicies> import(fit, corr_comb_dj);

    fit.fit(params, chisq, chisq_per_dof, corr_comb_j);
}  

void fit(jackknifeDistribution<taggedValueContainer<double,std::string> > &params, jackknifeDistributionD &chisq, jackknifeDistributionD &chisq_per_dof,
	 const correlationFunction<SimFitCoordGen,  jackknifeDistributionD> &corr_comb_j,
	 const correlationFunction<SimFitCoordGen,  doubleJackknifeDistributionD> &corr_comb_dj,
	 FitFuncType ffunc, const std::unordered_map<std::string,size_t> &param_map,
	 const int Lt, const double Ascale, const double Cscale){

  if(ffunc == FitFuncType::FSimGenOneState){
    typedef FitSimGenOneState FitFunc;
    FitFunc fitfunc(Lt, param_map.size(), Ascale, Cscale);
    return fit_ff<FitFunc>(params, chisq, chisq_per_dof, corr_comb_j, corr_comb_dj, fitfunc);
  }else if(ffunc == FitFuncType::FSimGenTwoState){
    typedef FitSimGenTwoState FitFunc;
    FitFunc fitfunc(Lt, param_map.size(), Ascale, Cscale);
    return fit_ff<FitFunc>(params, chisq, chisq_per_dof, corr_comb_j, corr_comb_dj, fitfunc);
  }else{
    assert(0);
  }
}


int main(const int argc, const char* argv[]){
  Args args;
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

  CMDline cmdline(argc,argv,2);

  std::vector<threeMomentum> pion_mom = { {1,1,1}, {-1,-1,-1},
					  {-1,1,1}, {1,-1,-1},
					  {1,-1,1}, {-1,1,-1},
					  {1,1,-1}, {-1,-1,1},
  
					  {3,1,1}, {-3,-1,-1},
					  {-3,1,1}, {3,-1,-1},
					  {3,-1,1}, {-3,1,-1},
					  {3,1,-1}, {-3,-1,1},
					  
					  {1,3,1}, {-1,-3,-1},
					  {-1,3,1}, {1,-3,-1},
					  {1,-3,1}, {-1,3,-1},
					  {1,3,-1}, {-1,-3,1},
					  
					  {1,1,3}, {-1,-1,-3},
					  {-1,1,3}, {1,-1,-3},
					  {1,-1,3}, {-1,1,-3},
					  {1,1,-3}, {-1,-1,3}
  };

  bubbleDataAllMomenta raw_bubble_gnd_gnd, raw_bubble_exc_exc, raw_bubble_gnd_exc;
  rawCorrelationFunction raw_data_gnd_gnd, raw_data_exc_exc, raw_data_gnd_exc;
  
  readRawDataOptions ropt;
  ropt.load_hdf5_data_checkpoint = cmdline.load_hdf5_data_checkpoint;
  ropt.load_hdf5_data_checkpoint_stub = cmdline.load_hdf5_data_checkpoint_stub;
  ropt.save_hdf5_data_checkpoint = cmdline.save_hdf5_data_checkpoint;
  ropt.save_hdf5_data_checkpoint_stub = cmdline.save_hdf5_data_checkpoint_stub;
  
  if(!cmdline.load_combined_data) readRawData(raw_bubble_gnd_gnd, raw_bubble_exc_exc, raw_bubble_gnd_exc,
					      raw_data_gnd_gnd, raw_data_exc_exc, raw_data_gnd_exc,
					      args.data_dir, args.figure_file_format, args.bubble_file_format,
					      args.Lt, args.tsep_pipi, args.tstep_pipi, pion_mom,
					      args.traj_start, args.traj_inc, args.traj_lessthan, ropt);
  

  //Get double-jack data
  doubleJackCorrelationFunction dj_data_gnd_gnd, dj_data_exc_exc, dj_data_gnd_exc;

  generateResampledDataOptions gopt; 
  gopt.load_combined_data = cmdline.load_combined_data;
  gopt.load_combined_data_file = cmdline.load_combined_data_file;
  gopt.save_combined_data = cmdline.save_combined_data;
  gopt.save_combined_data_file = cmdline.save_combined_data_file;

  generateResampledData(dj_data_gnd_gnd, dj_data_exc_exc, dj_data_gnd_exc,
			raw_bubble_gnd_gnd, raw_bubble_exc_exc, raw_bubble_gnd_exc,
			raw_data_gnd_gnd, raw_data_exc_exc, raw_data_gnd_exc,
			args.tsep_pipi, pion_mom, args.bin_size, 
			args.do_vacuum_subtraction, gopt);

  int nsample = dj_data_gnd_gnd.value(0).size();


  typedef std::unordered_map<std::string, std::string> SubFitFuncParameterMap;
  typedef std::unordered_map<std::string,size_t> ParamTagIdxMap;

  SubFitFuncParameterMap pipi2pt_gnd_gnd_oneexp_pmap, pipi2pt_exc_exc_oneexp_pmap, pipi2pt_gnd_exc_oneexp_pmap;
  ParamTagIdxMap param_map;

  typedef taggedValueContainer<double,std::string> Params;

  if(args.fitfunc == FitFuncType::FSimGenOneState){
    pipi2pt_gnd_gnd_oneexp_pmap = {  {"Asrc","Apipi_gnd"}, {"Asnk", "Apipi_gnd"}, {"E", "Epipi"}, {"Csys", "Cpipi_gnd_gnd"} };
    pipi2pt_exc_exc_oneexp_pmap = {  {"Asrc","Apipi_exc"}, {"Asnk", "Apipi_exc"}, {"E", "Epipi"}, {"Csys", "Cpipi_exc_exc"} };
    pipi2pt_gnd_exc_oneexp_pmap = {  {"Asrc","Apipi_gnd"}, {"Asnk", "Apipi_exc"}, {"E", "Epipi"}, {"Csys", "Cpipi_gnd_exc"} };
    
    param_map = { {"Apipi_gnd", 0}, {"Apipi_exc", 1},  {"Epipi",2}, {"Cpipi_gnd_gnd", 3}, {"Cpipi_exc_exc", 4},  {"Cpipi_gnd_exc", 5} };
  }else if(args.fitfunc == FitFuncType::FSimGenTwoState){
    pipi2pt_gnd_gnd_oneexp_pmap = {  {"Asrc0","Apipi_gnd_0"}, {"Asnk0", "Apipi_gnd_0"}, {"E0", "Epipi"},
				     {"Asrc1","Apipi_gnd_1"}, {"Asnk1", "Apipi_gnd_1"}, {"E1", "Eexc"},
				     {"Csys", "Cpipi_gnd_gnd"} };

    pipi2pt_exc_exc_oneexp_pmap = {  {"Asrc0","Apipi_exc_0"}, {"Asnk0", "Apipi_exc_0"}, {"E0", "Epipi"},
				     {"Asrc1","Apipi_exc_1"}, {"Asnk1", "Apipi_exc_1"}, {"E1", "Eexc"},
				     {"Csys", "Cpipi_exc_exc"} };

    pipi2pt_gnd_exc_oneexp_pmap = {  {"Asrc0","Apipi_gnd_0"}, {"Asnk0", "Apipi_exc_0"}, {"E0", "Epipi"},
				     {"Asrc1","Apipi_gnd_1"}, {"Asnk1", "Apipi_exc_1"}, {"E1", "Eexc"},
				     {"Csys", "Cpipi_gnd_exc"} };
    
    param_map = { {"Apipi_gnd_0", 0}, {"Apipi_exc_0", 1},  {"Epipi",2}, 
		  {"Apipi_gnd_1", 3}, {"Apipi_exc_1", 4},  {"Eexc",5},
		  {"Cpipi_gnd_gnd", 6}, {"Cpipi_exc_exc", 7},  {"Cpipi_gnd_exc", 8} };
  }else{
    assert(0);
  }

  std::unordered_map< SubFitFuncParameterMap const*, std::string > type_names = { {&pipi2pt_gnd_gnd_oneexp_pmap, "PiPi 2pt gnd/gnd"},
										  {&pipi2pt_exc_exc_oneexp_pmap, "PiPi 2pt exc/exc"},
										  {&pipi2pt_gnd_exc_oneexp_pmap, "PiPi 2pt gnd/exc"}  }; //just used for printing input data

  std::cout << "Selecting data in fit range" << std::endl;
  correlationFunction<SimFitCoordGen,  doubleJackknifeDistributionD> corr_comb_dj;

  for(int t=args.t_min;t<=args.t_max;t++)
    corr_comb_dj.push_back( SimFitCoordGen(t, &pipi2pt_gnd_gnd_oneexp_pmap, 2*args.tsep_pipi), dj_data_gnd_gnd.value(t) );
  
  for(int t=args.t_min;t<=args.t_max;t++)
    corr_comb_dj.push_back( SimFitCoordGen(t, &pipi2pt_exc_exc_oneexp_pmap, 2*args.tsep_pipi), dj_data_exc_exc.value(t) );

  for(int t=args.t_min;t<=args.t_max;t++)
    corr_comb_dj.push_back( SimFitCoordGen(t, &pipi2pt_gnd_exc_oneexp_pmap, 2*args.tsep_pipi), dj_data_gnd_exc.value(t) );


  correlationFunction<SimFitCoordGen,  jackknifeDistributionD> corr_comb_j(corr_comb_dj.size(), 
    [&](const int i){
	return correlationFunction<SimFitCoordGen,  jackknifeDistributionD>::ElementType( corr_comb_dj.coord(i), corr_comb_dj.value(i).toJackknife() );
    });

  std::cout << "Data in fit:" << std::endl;
  for(int i=0;i<corr_comb_j.size();i++){
    std::cout << type_names[corr_comb_j.coord(i).param_map] << " " << corr_comb_j.coord(i).t << " " << corr_comb_j.value(i) << std::endl;
  }

  std::cout << "Setting up guess" << std::endl;

  Params guess(&param_map);
  if(args.fitfunc == FitFuncType::FSimGenOneState){
    guess("Apipi_gnd") = 0.2;
    guess("Apipi_exc") = 0.2;
    guess("Epipi") = 0.4;
    guess("Cpipi_gnd_gnd") = 0;
    guess("Cpipi_exc_exc") = 0;
  }else if(args.fitfunc == FitFuncType::FSimGenTwoState){
    guess("Apipi_gnd_0") = 0.2;
    guess("Apipi_exc_0") = 0.2;
    guess("Apipi_gnd_1") = 0.2;
    guess("Apipi_exc_1") = 0.2;
    guess("Epipi") = 0.4;
    guess("Eexc") = 0.7;
    guess("Cpipi_gnd_gnd") = 0;
    guess("Cpipi_exc_exc") = 0;
    guess("Cpipi_gnd_exc") = 0;
  }

  std::cout << "Performing fit" << std::endl;

  jackknifeDistribution<Params> params(nsample, guess);
  jackknifeDistributionD chisq(nsample), chisq_per_dof(nsample);

  fit(params, chisq, chisq_per_dof,
      corr_comb_j, corr_comb_dj, args.fitfunc, param_map,
      args.Lt, args.Ascale, args.Cscale);

  std::cout << "Params:\n" << params << std::endl;
  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
  
  writeParamsStandard(params, "params.hdf5");
  writeParamsStandard(chisq, "chisq.hdf5");
  writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");

  std::cout << "Done\n";
  return 0;
}

