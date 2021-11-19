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

  //Resampled data
  doubleJackknifeCorrelationFunctionD pipi_dj;
  jackknifeCorrelationFunctionD pipi_j;

  getData(pipi_j, pipi_dj, args,cmdline);
  
  //Filter out the data that is to be fitted
  const int nsample = pipi_j.value(0).size();
  
  filterXrange<double> trange(args.t_min,args.t_max);
  
  doubleJackknifeCorrelationFunctionD pipi_dj_inrange;
  jackknifeCorrelationFunctionD pipi_j_inrange;

  std::cout << "Data in fit range:\n";
  for(int d=0;d<pipi_dj.size();d++)
    if(trange.accept(pipi_dj.coord(d),pipi_dj.value(d) )){
      pipi_dj_inrange.push_back(pipi_dj[d]);
      pipi_j_inrange.push_back(pipi_j[d]);
      std::cout << pipi_j_inrange.coord(pipi_j_inrange.size()-1) << " " << pipi_j_inrange.value(pipi_j_inrange.size()-1) << std::endl;
    }
  
  //Perform the fit
  pipiFitOptions opt; cmdline.Export(opt);
  std::pair<jackknifeDistributionD,jackknifeDistributionD> Epipi_and_const = fit(pipi_j_inrange, pipi_dj_inrange, args.fitfunc, args.correlated, args.Lt, args.tsep_pipi, args.Ascale, args.Cscale, opt);

  plot(pipi_j, Epipi_and_const.first, Epipi_and_const.second, 
       args.t_min, args.t_max, args.effective_energy, args.Lt, args.tsep_pipi, args.Ascale, args.Cscale);
  
  std::cout << "Done\n";
  return 0;
}
