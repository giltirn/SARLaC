#ifndef _FIT_KTOPIPI_GPARITY_SAMPLEAMA_ARGS_H_
#define _FIT_KTOPIPI_GPARITY_SAMPLEAMA_ARGS_H_

#define SAMPLEAMA_ARGS_MEMBERS						\
  ( std::string, data_dir_S )					\
  ( std::string, data_dir_C )					\
  ( int, Lt)							\
  ( int, tsep_pipi)						\
  ( std::vector<int>, tsep_k_pi)				\
  ( KtoPiPiFitFunc, fitfunc)					\
  ( bool, correlated )						\
  ( int, tmin_k_op)						\
  ( int, tmin_op_pi)						\
  ( int, bin_size )						\
  ( int, traj_inc )						\
  ( int, traj_start_S )						\
  ( int, traj_lessthan_S )					\
  ( int, traj_start_C )						\
  ( int, traj_lessthan_C )


struct SampleAMAargs{
  GENERATE_MEMBERS(SAMPLEAMA_ARGS_MEMBERS);

  SampleAMAargs(): data_dir_S("data_S"), data_dir_C("data_C"), fitfunc(KtoPiPiFitFunc::FitSeparate), correlated(false), Lt(64), tsep_pipi(4), tsep_k_pi(1,10), tmin_k_op(6), tmin_op_pi(4), traj_inc(1), bin_size(1), 
    traj_start_S(0), traj_lessthan_S(2), traj_start_C(0), traj_lessthan_C(2){}
} ;
GENERATE_PARSER(SampleAMAargs, SAMPLEAMA_ARGS_MEMBERS);


#endif
