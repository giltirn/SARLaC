#ifndef _FIT_KTOPIPI_GPARITY_ARGS_H_
#define _FIT_KTOPIPI_GPARITY_ARGS_H_

#define ARGS_MEMBERS						\
  ( std::string, data_dir )					\
  ( int, Lt)							\
  ( int, tsep_pipi)						\
  ( std::vector<int>, tsep_k_pi)				\
  ( int, tmin_k_op)						\
  ( int, tmin_op_pi)						\
  ( int, traj_start )						\
  ( int, traj_inc )						\
  ( int, traj_lessthan )					\
  ( double, AKscale)						\
  ( double, Apipiscale)

struct Args{
  GENERATE_MEMBERS(ARGS_MEMBERS)

  Args(): data_dir("data"), Lt(64), tsep_pipi(4), tsep_k_pi(1,10), tmin_k_op(6), tmin_op_pi(4), traj_start(0), traj_inc(1), traj_lessthan(2), AKscale(1), Apipiscale(1e13) {}
};
GENERATE_PARSER(Args, ARGS_MEMBERS)

#endif
