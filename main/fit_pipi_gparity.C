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
#include <fit_pipi_gparity/main.h>

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

  //Get the double-jackknife resampled data
  doubleJackCorrelationFunction pipi_dj = getData(args,cmdline);

  //Convert to single-jackknife
  jackknifeCorrelationFunction pipi_j(pipi_dj.size(),
				      [&pipi_dj](const int i){
					return typename jackknifeCorrelationFunction::ElementType(double(pipi_dj.coord(i)), pipi_dj.value(i).toJackknife());
				      });
  
  
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
  std::pair<jackknifeDistributionD,jackknifeDistributionD> Epipi_and_const = fit(pipi_j_inrange,pipi_dj_inrange,args,cmdline);

  plot(pipi_j,Epipi_and_const.first,Epipi_and_const.second,args,cmdline);
  
  std::cout << "Done\n";
  return 0;
}
