#include <fstream>
#include <algorithm>
#include <sstream>
#include <boost/timer/timer.hpp>

#include<fit.h>
#include<plot.h>
#include<data_series.h>
#include<common.h>
#include<parser.h>
#include <random.h>

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
#include <compare_asymm_symm_pipi/cmdline.h>
#include <compare_asymm_symm_pipi/args.h>
#include <compare_simple_correlators/compare.h>

int main(const int argc, const char* argv[]){
  ComparisonArgs args;
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
  
  ComparisonCMDline cmdline(argc,argv,2);

  PiPiProject *proj_src = getProjector(args.proj_src);
  PiPiProject *proj_snk = getProjector(args.proj_snk);
  PiPiMomAllow* allow = getMomPairFilter(args.allowed_mom);

  //Get the double-jackknife resampled data
  doubleJackCorrelationFunction pipi_dj_asymm, pipi_dj_symm;
  {
    CMDline c = cmdline.toCMDline(Asymmetric);
    Args a = args.toArgs(Asymmetric);
    pipi_dj_asymm = getData(*proj_src,*proj_snk,*allow,args.isospin,a,c);
  }
  {
    CMDline c = cmdline.toCMDline(Symmetric);
    Args a = args.toArgs(Symmetric);
    pipi_dj_symm = getData(*proj_src,*proj_snk,*allow,args.isospin,a,c);
  }
  
  delete proj_src; delete proj_snk; delete allow;

  //Convert to single-jackknife
  jackknifeCorrelationFunction pipi_j_asymm(pipi_dj_asymm.size(),
					    [&pipi_dj_asymm](const int i){
					      return typename jackknifeCorrelationFunction::ElementType(double(pipi_dj_asymm.coord(i)), pipi_dj_asymm.value(i).toJackknife());					      
					    });
  jackknifeCorrelationFunction pipi_j_symm(pipi_dj_symm.size(),
					   [&pipi_dj_symm](const int i){
					     return typename jackknifeCorrelationFunction::ElementType(double(pipi_dj_symm.coord(i)), pipi_dj_symm.value(i).toJackknife());					      
					   }); 
  
  compareRelativeDifferences(pipi_j_asymm,pipi_j_symm);
  
  std::cout << "Done\n";
  return 0;
}
