#ifndef _FIT_SIMPLE_SAMPLE_AMA_ARGS_H_
#define _FIT_SIMPLE_SAMPLE_AMA_ARGS_H_

#include<fit_simple/data_info.h>
#include<fit_simple/args.h>

enum SloppyExact {Sloppy, Exact};

//Ensembles are labelled S for the sloppy ensemble and C for the ensemble upon which the corrections are computed

#define DATA_INFO_SAMPLE_AMA_MEMBERS \
  ( ParserType, parser )	       \
  ( std::string, operation )   \
  ( TimeDependence, time_dep ) \
  ( std::string, sloppy_file_fmt_S ) \
  ( std::string, sloppy_file_fmt_C ) \
  ( std::string, exact_file_fmt_C )
  
//*_file_fmt should contain a '%d' which is replaced by the trajectory index
//operation is a math expression. Use x to represent the data. Use an empty string to leave as-is

struct DataInfoSampleAMA{
  GENERATE_MEMBERS(DATA_INFO_SAMPLE_AMA_MEMBERS);
  DataInfoSampleAMA(): parser(ParserStandard), operation(""), time_dep(TimeDepNormal), sloppy_file_fmt_S("sloppy_data_S.%d"), sloppy_file_fmt_C("sloppy_data_C.%d"), exact_file_fmt_C("exact_data_C.%d"){}

  DataInfo toDataInfo(const SloppyExact &se, const char ens) const{
    DataInfo out;
    out.parser = parser;
    out.operation = operation;
    out.time_dep = time_dep;
    if(se == Exact){
      assert(ens == 'C'); out.file_fmt = exact_file_fmt_C;
    }else{
      out.file_fmt = ens == 'S' ? sloppy_file_fmt_S : sloppy_file_fmt_C;
    }
    return out;
  }
    
};
GENERATE_PARSER(DataInfoSampleAMA, DATA_INFO_SAMPLE_AMA_MEMBERS);


#define ARGS_SAMPLE_AMA_MEMBERS \
  ( std::vector<DataInfoSampleAMA>, data )		\
  ( Combination, combination )			\
  ( TimeDependence, outer_time_dep )		\
  ( bool, correlated )				\
  ( FitFuncType, fitfunc)			\
  ( int, Lt)					\
  ( int, t_min)					\
  ( int, t_max)					\
  ( int, bin_size)				\
  ( int, traj_inc )				\
  ( int, sloppy_traj_start )				\
  ( int, sloppy_traj_lessthan )				\
  ( int, correction_traj_start )				\
  ( int, correction_traj_lessthan )

//Ranges of sloppy and exact data should not overlap

struct ArgsSampleAMA{
  GENERATE_MEMBERS(ARGS_SAMPLE_AMA_MEMBERS);

  ArgsSampleAMA(): Lt(64), combination(CombinationAverage), outer_time_dep(TimeDepNormal), correlated(false), t_min(0), t_max(32), data(1), bin_size(1), traj_inc(1),
		   sloppy_traj_start(0), sloppy_traj_lessthan(10), correction_traj_start(11), correction_traj_lessthan(15){}

  Args toArgs() const{
    Args out;
    out.data.resize(0);
    out.combination = combination;
    out.outer_time_dep = outer_time_dep;
    out.correlated = correlated;
    out.fitfunc = fitfunc;
    out.Lt = Lt;
    out.t_min = t_min;
    out.t_max = t_max;
    out.bin_size = bin_size;
    out.traj_inc = traj_inc;
    out.traj_start = 0; //these should not be used
    out.traj_lessthan = 0;
    return out;
  }
  
};
GENERATE_PARSER(ArgsSampleAMA, ARGS_SAMPLE_AMA_MEMBERS);

#endif
