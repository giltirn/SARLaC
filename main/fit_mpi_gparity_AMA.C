//#define BOOST_SPIRIT_X3_DEBUG
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
#include <sstream>

#include <generic_ET.h>

typedef NumericMatrix<distribution<double> > distributionMatrix;
typedef NumericVector<distribution<double> > distributionVector;

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
  ELEM( int, Lt) \
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

struct configData{
  distributionMatrix *all_data;
  int conf;

  configData(distributionMatrix &ad, const int c): all_data(&ad), conf(c){}
  inline int getLt() const{ return all_data->size(); }
};
namespace configData_grammar{
  namespace ascii = boost::spirit::x3::ascii;
  namespace x3 = boost::spirit::x3;
  using boost::fusion::at_c;

  struct set_entry{
    template <typename Context> void operator()(Context const& ctx) const{
      const auto &attr = x3::_attr(ctx);
      configData &con = x3::_val(ctx);
      distribution<double> &d = con.all_data->operator()(at_c<0>(attr),at_c<1>(attr));
      d.sample(con.conf) = at_c<2>(attr); //just the real part
    }
  };
  
  auto const line = x3::lexeme[x3::int_ > ' ' > x3::int_ > ' ' > x3::double_ > ' ' > x3::double_];
  x3::rule<struct line_rule_, configData> const line_rule = "line_rule";
  auto const line_rule_def = line[set_entry()];
  BOOST_SPIRIT_DEFINE(line_rule);
};

std::istream & operator>>(std::istream &is, configData &s){
  namespace ascii = boost::spirit::x3::ascii;
  namespace x3 = boost::spirit::x3;
  
  is >> std::noskipws;

  boost::spirit::istream_iterator f(is), l;
  const int nelems = s.getLt()*s.getLt();
  bool r = x3::phrase_parse(f, l, x3::repeat(nelems)[configData_grammar::line_rule], ascii::space, s);
  
  if(!r){
    throw std::runtime_error("Parsing ConfigData failed\n");
  }
  return is;
}

void read(configData &into, const std::string &fmt, const int traj){
  const std::string file = subsIdx(fmt,traj);
  std::ifstream is(file.c_str());
  if(!is.good()){ std::cout << "Could not open file \"" << file << "\"\n"; std::cout.flush(); exit(-1); }
  is >> into;
  if(is.fail() || is.bad()){ std::cout << "Error reading file \"" << file << "\"\n"; std::cout.flush(); exit(-1); }
  is.close();
}
void read(distributionMatrix &exact, distributionMatrix &sloppy, const SloppyExact &args, const int traj, const int sample_idx){
  if(args.include_data){
    configData e(exact, sample_idx);
    configData s(sloppy, sample_idx);    
    read(e, args.exact_fmt, traj);
    read(s, args.sloppy_fmt, traj);
  }
}


//Coordinate, parameters and param derivatives for aggregate fit func
struct Coord{
  double t; //time
  int idx; //data type index
};
std::ostream & operator<<(std::ostream &os, const Coord &c){
  os << "(" << c.t << "," << c.idx << ")";
  return os;
}

bool isZero(const distribution<double> &d){
  for(int i=0;i<d.size();i++) if(d.sample(i) != 0.) return false;
  return true;
}

std::vector<int> nonZeroSourceTimeSlices(const distributionMatrix &M, const int conf){
  const int Lt = M.size();
  std::vector<int> nonzero_slices;
  for(int tsrc=0;tsrc<Lt;tsrc++){
    bool allzero = true;
    for(int tsep=0;tsep<Lt;tsep++)
      if(M(tsrc,tsep).sample(conf) != 0. ){ allzero = false; break; }
    if(!allzero) nonzero_slices.push_back(tsrc);
  }
  return nonzero_slices;
}

distributionVector computeAMAcorrection(const distributionMatrix &sloppy, const distributionMatrix &exact){
  const int Lt = exact.size();
  const int nsample = sloppy(0,0).size();
  distributionVector out(Lt,distributionD(nsample,0.));

  for(int conf=0;conf<nsample;conf++){
    std::vector<int> nonzero_slices = nonZeroSourceTimeSlices(exact, conf);
    std::cout << "On conf " << conf << " exact is nonzero on tsrc={";
    for(int i=0;i<nonzero_slices.size();i++){
      const int tsrc = nonzero_slices[i];
      for(int tsep=0;tsep<Lt;tsep++)
	out[tsep].sample(conf) = out[tsep].sample(conf) + exact(tsrc,tsep).sample(conf) - sloppy(tsrc,tsep).sample(conf);      
      std::cout << tsrc << " ";
    }
    for(int tsep=0;tsep<Lt;tsep++)
      out[tsep].sample(conf) = out[tsep].sample(conf) / double(nonzero_slices.size());    
    std::cout << "}\n";
  }
  return out;
}

distributionMatrix timeReflect(const distributionMatrix &m){
  const int Lt = m.size();
  distributionMatrix out(Lt);
  for(int tsrc=0;tsrc<Lt;tsrc++)
    for(int tsep=0;tsep<Lt;tsep++)
      out(tsrc,Lt-tsep-1) = m(tsrc,tsep);
  return out;
}
distributionVector sourceTimeSliceAverage(const distributionMatrix &m){
  const int Lt = m.size();
  const int nsample = m(0,0).size();
  distributionVector out(Lt,distributionD(nsample));
  for(int tsrc=0;tsrc<Lt;tsrc++)
    for(int tsep=0;tsep<Lt;tsep++)
      out(tsep) = out(tsep) + m(tsrc,tsep);
    
  return out/double(Lt);
}



int main(const int argc, const char** argv){
  Args args;
  if(argc != 2){
    std::ofstream of("template.args");
    of << args;
    std::cout << "Wrote template argument file to template.args\n";
    return 0;
  }
  
  const std::string arg_file = argv[1];

  std::ifstream arg_f(arg_file.c_str());
  
  arg_f >> args;

  std::cout << "Read arguments: \n" << args << std::endl;

  const int ntraj = (args.traj_lessthan - args.traj_start)/args.traj_inc;
  assert(ntraj > 0);
  
  distributionMatrix PP_LW_FF_exact_data(args.Lt, distributionD(ntraj));
  distributionMatrix PP_LW_FF_sloppy_data(args.Lt, distributionD(ntraj));
  distributionMatrix PP_LW_BB_exact_data(args.Lt, distributionD(ntraj));
  distributionMatrix PP_LW_BB_sloppy_data(args.Lt, distributionD(ntraj));
  
#pragma omp parallel for
  for(int i=0;i<ntraj;i++){
    const int c = args.traj_start + i*args.traj_inc;
    read(PP_LW_FF_exact_data, PP_LW_FF_sloppy_data, args.PP_LW.FF_data, c, i);
    read(PP_LW_BB_exact_data, PP_LW_BB_sloppy_data, args.PP_LW.BB_data, c, i);  
  }
  PP_LW_BB_exact_data = timeReflect(PP_LW_BB_exact_data);
  PP_LW_BB_sloppy_data = timeReflect(PP_LW_BB_sloppy_data);

  distributionVector PP_LW_FF_sloppy_avg = sourceTimeSliceAverage(PP_LW_FF_sloppy_data);
  distributionVector PP_LW_BB_sloppy_avg = sourceTimeSliceAverage(PP_LW_BB_sloppy_data);

  std::cout << "Computing PP_LW_FF_correction\n";
  distributionVector PP_LW_FF_correction = computeAMAcorrection(PP_LW_FF_sloppy_data, PP_LW_FF_exact_data);
  std::cout << "Computing PP_LW_BB_correction\n";
  distributionVector PP_LW_BB_correction = computeAMAcorrection(PP_LW_BB_sloppy_data, PP_LW_BB_exact_data);
  publicationPrint printer;
  for(int t=0;t<args.Lt;t++){
    printer << "PP_LW_FF_sloppy_avg[" << t << "] = " << PP_LW_FF_sloppy_avg[t] << std::endl;
    printer << "PP_LW_BB_sloppy_avg[" << t << "] = " << PP_LW_BB_sloppy_avg[t] << std::endl;
    printer << "PP_LW_FF_correction[" << t << "] = " << PP_LW_FF_correction[t] << std::endl;
    printer << "PP_LW_BB_correction[" << t << "] = " << PP_LW_BB_correction[t] << std::endl;
  }

  return 0;
}
