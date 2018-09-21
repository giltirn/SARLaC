#include <fstream>
#include <regex>
#include <vector>
#include <map>
#include <set>

#include <utils.h>
#include <common.h>

#include <pipi_common/pipi_common.h>

#include <fit_pipi_sigma_sim_gparity/fitfunc.h>
#include <fit_pipi_sigma_sim_gparity/fit.h>

#include<fit_pipi_sigma_sim_gparity_sampleAMA/args.h>
#include<fit_pipi_sigma_sim_gparity_sampleAMA/cmdline.h>

#include<fit_pipi_sigma_sim_gparity_sampleAMA/data_map.h>
#include<fit_pipi_sigma_sim_gparity_sampleAMA/set_operations.h>
#include<fit_pipi_sigma_sim_gparity_sampleAMA/data_map_parse.h>
#include<fit_pipi_sigma_sim_gparity_sampleAMA/data_map_utils.h>
#include<fit_pipi_sigma_sim_gparity_sampleAMA/enums.h>
#include<fit_pipi_sigma_sim_gparity_sampleAMA/utils.h>

#include<fit_pipi_sigma_sim_gparity_sampleAMA/pipi2pt_datasets.h>
#include<fit_pipi_sigma_sim_gparity_sampleAMA/sigma2pt_datasets.h>
#include<fit_pipi_sigma_sim_gparity_sampleAMA/pipitosigma_datasets.h>

#include<fit_pipi_sigma_sim_gparity_sampleAMA/read_pipibubble.h>
#include<fit_pipi_sigma_sim_gparity_sampleAMA/read_sigmabubble.h>
#include<fit_pipi_sigma_sim_gparity_sampleAMA/read_pipitosigma.h>
#include<fit_pipi_sigma_sim_gparity_sampleAMA/read_sigma2pt.h>
#include<fit_pipi_sigma_sim_gparity_sampleAMA/read_pipi2pt.h>

//TODO: Rename outer config to outer sample in code; former is used only when checking the config file
//      typedef int to OuterSample to make code easier to read
using namespace CPSfit;

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


  assert(args.bin_size == 1); //doesn't yet support non-unit bin size

  DataMap dmap(args.config_map);
  const int nsample_full = dmap.size();

  /*-------------------------------------------
    Ground operator pipi 2pt and bubble
    -------------------------------------------*/

  std::vector<threeMomentum> gnd_pion_mom = { {1,1,1}, {-1,-1,-1},
					  {-1,1,1}, {1,-1,-1},
					  {1,-1,1}, {-1,1,-1},
					  {1,1,-1}, {-1,-1,1} };

  std::map<SubensTag, std::set<int> > ground_pipi_subsets = getGroundPiPiSubsets(dmap);
  std::map<DataTag, std::map<int, DataLocationInfo const*> > ground_pipi_data_map = getGroundPiPiDataSubsets(ground_pipi_subsets, dmap);

  bubbleDataDoubleJackAllMomenta ground_pipi_bub;
  getPiPiBubble(ground_pipi_bub, args.Lt, args.tsep_pipi, ground_pipi_data_map, nsample_full, PiPiProjector::A1momSet111, PiPiProjector::A1momSet111);

  doubleJackCorrelationFunction ground_pipi2pt;
  getPiPi2pt(ground_pipi2pt, args.Lt, args.tsep_pipi, args.tstep_pipi2pt, ground_pipi_data_map, nsample_full, PiPiProjector::A1momSet111, PiPiProjector::A1momSet111);

  /*-------------------------------------------
    Ground operator pipi->sigma
    -------------------------------------------*/

  std::map<SubensTag, std::set<int> > pipi_to_sigma_subsets = getPiPiToSigmaSubsets(dmap);
  std::map<DataTag, std::map<int, DataLocationInfo const*> > pipi_to_sigma_data_map = getPiPiToSigmaDataSubsets(pipi_to_sigma_subsets, dmap);

  doubleJackCorrelationFunction pipi_to_sigma;
  getPiPiToSigma(pipi_to_sigma, args.Lt, args.tsep_pipi, args.tstep_pipitosigma, pipi_to_sigma_data_map, nsample_full);

  /*-------------------------------------------
    Sigma -> Sigma  and sigma bubble
    -------------------------------------------*/

  std::map<SubensTag, std::set<int> > sigma2pt_subsets = getSigma2ptSubsets(dmap);
  std::map<DataTag, std::map<int, DataLocationInfo const*> > sigma2pt_data_map = getSigma2ptDataSubsets(sigma2pt_subsets, dmap);

  doubleJackCorrelationFunction sigma2pt;
  getSigma2pt(sigma2pt, args.Lt, sigma2pt_data_map, nsample_full);

  sigmaSelfContractionDoubleJack sigma_bub;
  getSigmaBubble(sigma_bub, args.Lt, sigma2pt_data_map, nsample_full);

  
  //Do vacuum subtractions
  if(args.do_vacuum_subtraction){
    performPiPi2ptVacuumSubtraction(ground_pipi2pt, ground_pipi2pt, ground_pipi_bub, args.Lt, args.tsep_pipi, PiPiProjector::A1momSet111, PiPiProjector::A1momSet111);
    performPiPiToSigmaVacuumSubtraction(pipi_to_sigma, pipi_to_sigma, sigma_bub, ground_pipi_bub, args.Lt, gnd_pion_mom);
    performSigma2ptVacuumSubtraction(sigma2pt, sigma2pt, sigma_bub);
  }

  //Fold data

  ground_pipi2pt = fold(ground_pipi2pt, 2*args.tsep_pipi);
  pipi_to_sigma = fold(pipi_to_sigma, args.tsep_pipi);
  sigma2pt = fold(sigma2pt, 0);


  //Build the combined data set
  simFitCorrFuncDJ corr_comb_dj;
  doubleJackCorrelationFunction* dsets[3] = { &ground_pipi2pt, &pipi_to_sigma, &sigma2pt };
  SimFitType dtype[3] = { SimFitType::PiPi2pt, SimFitType::PiPiToSigma, SimFitType::Sigma2pt };
  bool dincl[3] = { true, true, true };

  for(int d=0;d<3;d++)
    if(dincl[d])
      for(int t=args.t_min; t<=args.t_max; t++) 
	corr_comb_dj.push_back( simFitCorrFuncDJ::ElementType(  SimFitCoord(dtype[d],t), dsets[d]->value(t) ) );
  
  simFitCorrFuncJ corr_comb_j(corr_comb_dj.size(), [&](const int i){ return simFitCorrFuncJ::ElementType(corr_comb_dj.coord(i), corr_comb_dj.value(i).toJackknife()); });

  std::cout << "Data in fit range:\n";
  for(int i=0;i<corr_comb_j.size();i++){
    std::cout << corr_comb_j.coord(i).type << " " << corr_comb_j.coord(i).t << " " << corr_comb_j.value(i) << std::endl;
  }

  //Run the fit
  SimFitArgs fargs;
  args.transfer(fargs);
  cmdline.transfer(fargs);

  fit(corr_comb_j, corr_comb_dj, fargs);

  std::cout << "Finished" << std::endl;
  return 0;
}


  // std::vector<SubensTag> data_ens = { SubensTag::AsymmOnly, SubensTag::Correction, SubensTag::SymmOnly };
  
  // for(auto etag = data_ens.begin(); etag != data_ens.end(); etag++){
  //   std::cout << "Pruning ground pipi subens " << *etag << " for bin size " << bin_size << std::endl;
  //   ground_pipi_subsets[*etag] = binningPrune(ground_pipi_subsets[*etag], bin_size);
  // }


  // //TEST
  // {
  //   for(auto ensset = sigma2pt_subsets.begin(); ensset != sigma2pt_subsets.end(); ensset++){
  //     std::set<int> &theset = ensset->second;
  //     int i=0;
  //     auto it = theset.begin();
  //     while(it != theset.end()){
  // 	if(i >= 100)
  // 	  it = theset.erase(it);	
  // 	else
  // 	  ++it;	
  // 	++i;
  //     }
  //   }

  //   for(auto it = sigma2pt_subsets.begin(); it != sigma2pt_subsets.end(); ++it)
  //     std::cout << "sigma 2pt restricted " << it->first << " " << it->second.size() << std::endl;
  
  // }
  // //TEST
  // figureData pipi_to_sigma_sloppy_C_raw = readData(pipi_to_sigma_correction, "sigma", dmap, pipi_to_sigma_asymm_parser);
  
  
  // figureData pipi_to_sigma_exact_C_raw = readData(pipi_to_sigma_correction, "extended", dmap, pipi_to_sigma_symm_parser);

  // rawCorrelationFunction pipi_to_sigma_sloppy_C_srcavg_raw = sourceAverage(pipi_to_sigma_sloppy_C_raw);
  // rawCorrelationFunction pipi_to_sigma_exact_C_srcavg_raw = sourceAverage(pipi_to_sigma_sloppy_C_raw);



  // StandardParseReal parse("traj_<CONF>_pipitosigma_sigmawdagmom111_pionmom222_v2", Lt);
  // NumericTensor<rawDataDistributionD,2> pipi_to_sigma_sloppy_C = readData(pipi_to_sigma_correction, "sigma", dmap, parse);



  


  // //Get the outer configurations for the symmpi sample-AMA correction
  // std::vector<int> sloppy_configs;
  // std::vector<int> exact_configs;
  
  // std::string sloppy_tag = "orig";
  // std::string exact_tag = "correction";

  // getSampleAMAconfigurations(sloppy_configs, exact_configs, sloppy_tag, exact_tag, dmap);
  
  // std::cout << "Exact configs:\n";
  // for(int i=0;i<exact_configs.size();i++){
  //   std::cout << dmap(exact_configs[i], exact_tag) << std::endl;
  // }
  // std::cout << "Found " << exact_configs.size() << " exact configs\n";

  // std::cout << "Sloppy configs:\n";
  // for(int i=0;i<sloppy_configs.size();i++){
  //   std::cout << dmap(sloppy_configs[i], sloppy_tag) << std::endl;
  // }
  // std::cout << "Found " << sloppy_configs.size() << " sloppy configs\n";



  // for(int i=0;i<dmap.size();i++){
  //   std::cout << i << " ";
  //   for(int j=0;j<dmap[i].size();j++)
  //     std::cout << "(" << dmap[i][j] << ") ";
  //   std::cout << std::endl;

  // }

  // const std::vector< std::pair<int,int> > &extended = dmap.getLocationsWithTag("extended");
  
  // for(int i=0;i<extended.size();i++)
  //   std::cout << dmap(extended[i]) << std::endl;


//   std::cout << "Finished" << std::endl;
//   return 0;
// }

//using namespace CPSfit;
