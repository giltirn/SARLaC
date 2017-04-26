#include <fstream>
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
#include <sstream>

#include <generic_ET.h>

//Fit the pion mass from a simultaneous fit to multiple pseudoscalar two-point functions
#define SLOPPY_EXACT_MEMBERS \
  ELEM( std::string, sloppy_fmt )     \
  ELEM( std::string, exact_fmt )      \
  ELEM( bool, include_data )

struct SloppyExact{
#define STRUCT_ARGS SLOPPY_EXACT_MEMBERS
#include<struct_gen.incl>

  SloppyExact(): include_data(false){}
};
#define STRUCT_TYPE SloppyExact
#define STRUCT_ARGS SLOPPY_EXACT_MEMBERS
#include<parser_gen.incl>


#define TWOPOINTFUNCTION_MEMBERS \
  ELEM( SloppyExact, FF_data )     \
  ELEM( SloppyExact, BB_data )

struct TwoPointFunction{
#define STRUCT_ARGS TWOPOINTFUNCTION_MEMBERS
#include<struct_gen.incl>

};
#define STRUCT_TYPE TwoPointFunction
#define STRUCT_ARGS TWOPOINTFUNCTION_MEMBERS
#include<parser_gen.incl>


#define ARGS_MEMBERS \
  ELEM( TwoPointFunction, PP_LW )   \
  ELEM( TwoPointFunction, AP_LW )   \
  ELEM( int, traj_start ) \
  ELEM( int, traj_inc ) \
  ELEM( int, traj_lessthan )

struct Args{
#define STRUCT_ARGS ARGS_MEMBERS
#include<struct_gen.incl>

  Args(): traj_start(0), traj_inc(1), traj_lessthan(2){}
};
#define STRUCT_TYPE Args
#define STRUCT_ARGS ARGS_MEMBERS
#include<parser_gen.incl>




//Coordinate, parameters and param derivatives for aggregate fit func
struct Coord{
  double t; //time
  int idx; //data type index
};
std::ostream & operator<<(std::ostream &os, const Coord &c){
  os << "(" << c.t << "," << c.idx << ")";
  return os;
}

int main(const int argc, const char** argv){
  Args args;
  if(argc != 2){
    std::ofstream of("template.args");
    of << args;
    std::cout << "Wrote template argument file to template.args\n";
    return 0;
  }
  
  std::string arg_file = argv[1];

  std::ifstream arg_f(arg_file.c_str());
  
  arg_f >> args;

  std::cout << "Read arguments: \n" << args << std::endl;

  return 0;
}
