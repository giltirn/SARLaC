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
  doubleJackCorrelationFunction pipi_dj_vacsubbed = getData(args,cmdline);
  
  //Filter out the data that is to be fitted
  const int nsample = pipi_dj_vacsubbed.value(0).size();
  
  filterXrange<double> trange(args.t_min,args.t_max);
  
  doubleJackCorrelationFunction pipi_dj_vacsubbed_inrange;
  for(int d=0;d<pipi_dj_vacsubbed.size();d++)
    if(trange.accept(pipi_dj_vacsubbed.coord(d),pipi_dj_vacsubbed.value(d) ))
      pipi_dj_vacsubbed_inrange.push_back(pipi_dj_vacsubbed[d]);
    
  const int ndata_fit = pipi_dj_vacsubbed_inrange.size();
  
  //Get the single-jackknife data
  jackknifeCorrelationFunction pipi_j_vacsubbed_inrange(ndata_fit,
							[&pipi_dj_vacsubbed_inrange](const int i)
							{
							  return typename jackknifeCorrelationFunction::ElementType(double(pipi_dj_vacsubbed_inrange.coord(i)), pipi_dj_vacsubbed_inrange.value(i).toJackknife());
							}
							);

  //Perform the fit
  fit(pipi_j_vacsubbed_inrange,pipi_dj_vacsubbed_inrange,args,cmdline);

  std::cout << "Done\n";
  return 0;
}
