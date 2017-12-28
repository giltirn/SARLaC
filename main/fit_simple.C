#include<common_defs.h>
#include<fit_wrapper.h>
#include<parser.h>
#include<expression_parse.h>
#include<fitfunc.h>
#include<effective_mass.h>
#include<plot.h>

#include<fit_simple/cmdline.h>
#include<fit_simple/args.h>
#include<fit_simple/read_data.h>
#include<fit_simple/fit.h>
#include<fit_simple/main.h>

//Basic fitting
int main(const int argc, const char** argv){
  CMDline cmdline(argc,argv,2);

  Args args;
  if(argc < 2){
    std::ofstream of("template.args");
    (std::cout << "No parameter file provided: writing template to 'template.args' and exiting\n").flush();
    of << args;
    return 1;
  }    
  
  parse(args, argv[1]);

  std::vector<rawDataCorrelationFunctionD> channels(args.data.size());
  for(int i=0;i<channels.size();i++) readData(channels[i], args.data[i], args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan);

  rawDataCorrelationFunctionD data;
  applyCombination(data,channels,args.combination);
  applyTimeDep(data, args.outer_time_dep, args.Lt);

  bin(data, args.bin_size);

  fit(data, args, cmdline);
}

