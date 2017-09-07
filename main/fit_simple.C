#include<common_defs.h>
#include<fit_wrapper.h>
#include<parser.h>
#include<expression_parse.h>
#include<fitfunc.h>
#include<effective_mass.h>
#include<plot.h>

#include<fit_simple/args.h>
#include<fit_simple/read_data.h>
#include<fit_simple/fit.h>
#include<fit_simple/main.h>

//Basic fitting
int main(const int argc, const char** argv){
  assert(argc == 2);

  Args args;
  parse(args, argv[1]);

  std::vector<rawDataCorrelationFunctionD> channels(args.data.size());
  for(int i=0;i<channels.size();i++) readData(channels[i], args.data[i], args);

  rawDataCorrelationFunctionD data;
  applyCombination(data,channels,args.combination);
  applyTimeDep(data, args.outer_time_dep, args);

  fit(data, args);
}

