#include<common_defs.h>
#include<fit_wrapper.h>
#include<parser.h>
#include<expression_parse.h>
#include<fitfunc.h>
#include<effective_mass.h>
#include<plot.h>

#include<fit_simple/read_data.h>
#include<compare_simple_correlators/args.h>
#include<compare_simple_correlators/compare.h>

int main(const int argc, const char** argv){
  Args args;
  if(argc < 2){
    std::ofstream of("template.args");
    (std::cout << "No parameter file provided: writing template to 'template.args' and exiting\n").flush();
    of << args;
    return 1;
  }    
  
  parse(args, argv[1]);

  rawDataCorrelationFunctionD data[2];
  
  for(int aa=0;aa<2;aa++){
    const std::vector<DataInfo> &data_args = aa == 0 ? args.data_A : args.data_B;
    
    std::vector<rawDataCorrelationFunctionD> channels(data_args.size());
    for(int i=0;i<channels.size();i++) readData(channels[i], data_args[i], args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan);

    applyCombination(data[aa],channels, aa == 0 ? args.combination_A : args.combination_B);
    applyTimeDep(data[aa], aa == 0 ? args.outer_time_dep_A : args.outer_time_dep_B, args.Lt);
  }

  if(data[0].size() != data[1].size()) error_exit(std::cout << "Error: correlation functions must be of the same size");
  int sz = data[0].size();

  jackknifeCorrelationFunctionD data_j[2];
  for(int aa=0;aa<2;aa++)
    data_j[aa] = jackknifeCorrelationFunctionD(sz, [&](const int i){ return typename jackknifeCorrelationFunctionD::ElementType(data[aa].coord(i), jackknifeDistributionD(data[aa].value(i))); });

  compareRelativeDifferences(data_j[0],data_j[1]);

  std::cout << "Done\n";
  return 0;
}

