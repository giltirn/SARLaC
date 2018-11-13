#include <pipi_common/pipi_common.h>
#include <pipi_common/analyze_chisq.h>

using namespace CPSfit;

#include <fit_pipi_gparity/args.h>
#include <fit_pipi_gparity/cmdline.h>

#include <fit_pipi_sigma_sim_gparity/fitfunc.h>
#include <fit_pipi_sigma_sim_gparity/fit.h>
#include <fit_pipi_sigma_sim_gparity/args.h>
#include <fit_pipi_sigma_sim_gparity/cmdline.h>
#include <fit_pipi_sigma_sim_gparity/read_data.h>

struct rawData{
  rawCorrelationFunction pipi_raw, pipi_to_sigma_raw, sigma2pt_raw;
  bubbleDataAllMomenta pipi_self_data;
  sigmaSelfContraction sigma_self_data;

  void readDataFromOrigFiles(const std::string &data_dir, 
		const std::string &pipi2pt_figure_file_fmt, 
		const std::string &sigma2pt_file_fmt, 
		const std::string &pipitosigma_file_fmt, 
		const std::string &pipi_bubble_file_fmt, 
		const std::string &sigma_bubble_file_fmt,
		const int tsep_pipi,
		const int tstep_pipi2pt, int tstep_pipitosigma,
		const int Lt, const int traj_start, const int traj_inc, const int traj_lessthan,
		bool compute_pipitosigma_disconn_ReRe){
        readData(pipi_raw,pipi_to_sigma_raw,sigma2pt_raw,pipi_self_data,sigma_self_data,
		 data_dir, 
		 pipi2pt_figure_file_fmt, sigma2pt_file_fmt, pipitosigma_file_fmt,
		 pipi_bubble_file_fmt, sigma_bubble_file_fmt,
		 tsep_pipi, tstep_pipi2pt, tstep_pipitosigma,
		 Lt, traj_start, traj_inc, traj_lessthan, compute_pipitosigma_disconn_ReRe);
  }
  void readDataFromCheckpoint(const std::string &file){
    readCheckpoint(pipi_raw,pipi_to_sigma_raw,sigma2pt_raw,pipi_self_data,sigma_self_data,file);
  }
  void saveCheckpoint(const std::string &file) const{
    writeCheckpoint(file,pipi_raw,pipi_to_sigma_raw,sigma2pt_raw,pipi_self_data,sigma_self_data);
  }
};

void getResampledData(simFitCorrFuncJ &corr_comb_j, simFitCorrFuncDJ &corr_comb_dj, const rawData &raw_data,
		      const int Lt, const int tsep_pipi, const int bin_size, bool do_vacuum_subtraction,
		      bool include_pipi_2pt, bool include_pipi_to_sigma, bool include_sigma_2pt){

  //Get double-jack data
  auto pipi_to_sigma_j = binResample<jackknifeCorrelationFunction>(raw_data.pipi_to_sigma_raw, bin_size);
  auto pipi_to_sigma_dj = binResample<doubleJackCorrelationFunction>(raw_data.pipi_to_sigma_raw, bin_size);

  auto sigma2pt_j = binResample<jackknifeCorrelationFunction>(raw_data.sigma2pt_raw, bin_size);
  auto sigma2pt_dj = binResample<doubleJackCorrelationFunction>(raw_data.sigma2pt_raw, bin_size);

  auto pipi_j = binResample<jackknifeCorrelationFunction>(raw_data.pipi_raw, bin_size);
  auto pipi_dj = binResample<doubleJackCorrelationFunction>(raw_data.pipi_raw, bin_size);
  
  //Compute vacuum subtractions
  if(do_vacuum_subtraction){
    { //Pipi->sigma
      PiPiProjectorA1Basis111 proj_pipi;
      bubbleData pipi_self_proj = projectSourcePiPiBubble(raw_data.pipi_self_data, proj_pipi);
      pipi_to_sigma_j = pipi_to_sigma_j - computePiPiToSigmaVacSub<jackknifeCorrelationFunction>(raw_data.sigma_self_data, pipi_self_proj, bin_size);
      pipi_to_sigma_dj = pipi_to_sigma_dj - computePiPiToSigmaVacSub<doubleJackCorrelationFunction>(raw_data.sigma_self_data, pipi_self_proj, bin_size);
    }
    { //sigma->sigma
      sigma2pt_j = sigma2pt_j - computeSigmaVacSub<jackknifeCorrelationFunction>(raw_data.sigma_self_data, bin_size);
      sigma2pt_dj = sigma2pt_dj - computeSigmaVacSub<doubleJackCorrelationFunction>(raw_data.sigma_self_data, bin_size);
    }
    { //Pipi->pipi
      pipi_j = pipi_j - computePiPi2ptVacSub<jackknifeCorrelationFunction>(raw_data.pipi_self_data, bin_size, tsep_pipi);
      pipi_dj = pipi_dj - computePiPi2ptVacSub<doubleJackCorrelationFunction>(raw_data.pipi_self_data, bin_size, tsep_pipi);
    }
  }
  
  //Fold data
  pipi_to_sigma_j = fold(pipi_to_sigma_j, tsep_pipi);
  pipi_to_sigma_dj = fold(pipi_to_sigma_dj, tsep_pipi);

  sigma2pt_j = fold(sigma2pt_j, 0);
  sigma2pt_dj = fold(sigma2pt_dj, 0);

  pipi_j = fold(pipi_j, 2*tsep_pipi);
  pipi_dj = fold(pipi_dj, 2*tsep_pipi);
  
  //Build the combined data set
  jackknifeCorrelationFunction const* dsets_j[3] = { &pipi_j, &pipi_to_sigma_j, &sigma2pt_j };
  doubleJackCorrelationFunction const* dsets_dj[3] = { &pipi_dj, &pipi_to_sigma_dj, &sigma2pt_dj };

  SimFitType dtype[3] = { SimFitType::PiPi2pt, SimFitType::PiPiToSigma, SimFitType::Sigma2pt };
  bool dincl[3] = { include_pipi_2pt, include_pipi_to_sigma, include_sigma_2pt };

  for(int d=0;d<3;d++)
    if(dincl[d])
      for(int t=0; t<Lt; t++){ 
	corr_comb_j.push_back(SimFitCoord(dtype[d],t), dsets_j[d]->value(t));
	corr_comb_dj.push_back(SimFitCoord(dtype[d],t), dsets_dj[d]->value(t));
      }
}		      
		      



int main(const int argc, const char* argv[]){
  PiPiSigmaSimArgs args;
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

  PiPiSigmaSimCMDline cmdline(argc,argv,2);

  simFitCorrFuncJ corr_comb_j_all;
  simFitCorrFuncDJ corr_comb_dj_all;

  if(cmdline.load_resampled_data){
    HDF5reader rd(cmdline.load_resampled_data_file);
    read(rd,corr_comb_j_all, "corr_comb_j_all");
    read(rd,corr_comb_dj_all, "corr_comb_dj_all");
  }else{
    rawData raw_data;
    if(cmdline.load_checkpoint){
      raw_data.readDataFromCheckpoint(cmdline.load_checkpoint_file);
    }else{
      raw_data.readDataFromOrigFiles(args.data_dir, 
				     args.pipi2pt_figure_file_fmt, args.sigma2pt_file_fmt, args.pipitosigma_file_fmt,
				     args.pipi_bubble_file_fmt, args.sigma_bubble_file_fmt,
				     args.tsep_pipi, args.tstep_pipi2pt, args.tstep_pipitosigma,
				     args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan, !cmdline.use_pipitosigma_disconn_complex_prod);
    }
    if(cmdline.save_checkpoint){
      raw_data.saveCheckpoint(cmdline.save_checkpoint_file);
    }
    getResampledData(corr_comb_j_all, corr_comb_dj_all, raw_data, args.Lt, args.tsep_pipi, args.bin_size, args.do_vacuum_subtraction,
		     cmdline.include_pipi_2pt, cmdline.include_pipi_to_sigma, cmdline.include_sigma_2pt);
  }
  if(cmdline.save_resampled_data){
    HDF5writer wr(cmdline.save_resampled_data_file);
    write(wr,corr_comb_j_all, "corr_comb_j_all");
    write(wr,corr_comb_dj_all, "corr_comb_dj_all");
  }

  //Get data in fit range
  bool include_type[3] = {false,false,false};
  include_type[(int)SimFitType::PiPi2pt] = cmdline.include_pipi_2pt;
  include_type[(int)SimFitType::PiPiToSigma] = cmdline.include_pipi_to_sigma;
  include_type[(int)SimFitType::Sigma2pt] = cmdline.include_sigma_2pt;
  
  simFitCorrFuncJ corr_comb_j;
  simFitCorrFuncDJ corr_comb_dj;
  for(int tt=0;tt<corr_comb_j_all.size();tt++){
    auto coord = corr_comb_j_all.coord(tt);
    if(include_type[(int)coord.type] && 
       coord.t >= args.t_min && 
       coord.t <= args.t_max){
      corr_comb_j.push_back(corr_comb_j_all[tt]);
      corr_comb_dj.push_back(corr_comb_dj_all[tt]);
    }
  }

  std::cout << "Data in fit range:\n";
  for(int i=0;i<corr_comb_j.size();i++){
    std::cout << corr_comb_j.coord(i).type << " " << corr_comb_j.coord(i).t << " " << corr_comb_j.value(i) << std::endl;
  }

  //Run the fit
  SimFitArgs fargs;
  args.transfer(fargs);
  cmdline.transfer(fargs);

  fit(corr_comb_j, corr_comb_dj, fargs);

  std::cout << "Done\n";
  return 0;
}
