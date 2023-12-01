#if 1

//Code needs to be updated and fixed

int main(void){
  return 0;
}


#else

#include <fstream>
#include <algorithm>
#include <sstream>
#include <boost/timer/timer.hpp>

#include<fit.h>
#include<plot.h>
#include<data_series.h>
#include<common.h>
#include<parser.h>
#include <random.h>

#include <pipi_common/pipi_common.h>

using namespace CPSfit;

#include <compare_asymm_symm_pipi/cmdline.h>
#include <compare_asymm_symm_pipi/args.h>
#include <compare_simple_correlators/compare.h>

int main(const int argc, const char* argv[]){
  ComparisonArgs args;
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
  
  ComparisonCMDline cmdline(argc,argv,2);

  std::unique_ptr<PiPiProjectorBase> proj_src( getProjector(args.proj_src, {0,0,0}) );
  std::unique_ptr<PiPiProjectorBase> proj_snk( getProjector(args.proj_snk, {0,0,0}) );

  getResampledPiPi2ptDataOpts opts;
#define CP(A) opts.A = cmdline.A
  CP(load_hdf5_data_checkpoint);
  CP(save_hdf5_data_checkpoint);    
#undef CP
  //Get the double-jackknife resampled data
  jackknifeCorrelationFunction pipi_j_asymm, pipi_j_symm;
  doubleJackCorrelationFunction pipi_dj_asymm, pipi_dj_symm;
  {
    readFigureStationaryPolicy ffn(false);
    readBubbleStationaryPolicy bfn_src(false,Source), bfn_snk(false,Sink);

    opts.load_hdf5_data_checkpoint_stub = cmdline.load_hdf5_data_checkpoint_asymm_stub;
    opts.save_hdf5_data_checkpoint_stub = cmdline.save_hdf5_data_checkpoint_asymm_stub;

    getResampledPiPi2ptData(pipi_j_asymm, pipi_dj_asymm,
			    *proj_src,*proj_snk, args.isospin, args.Lt, args.tsep_pipi, args.tstep_pipi, args.do_vacuum_subtraction,
			    args.data_dir_asymm, args.traj_start, args.traj_inc, args.traj_lessthan, args.bin_size,
			    ffn, bfn_src, bfn_snk, "", opts); 
  }
  {
    readFigureStationaryPolicy ffn(true);
    readBubbleStationaryPolicy bfn_src(true,Source), bfn_snk(true,Sink);

    opts.load_hdf5_data_checkpoint_stub = cmdline.load_hdf5_data_checkpoint_symm_stub;
    opts.save_hdf5_data_checkpoint_stub = cmdline.save_hdf5_data_checkpoint_symm_stub;

    getResampledPiPi2ptData(pipi_j_symm, pipi_dj_symm,
			    *proj_src,*proj_snk, args.isospin, args.Lt, args.tsep_pipi, args.tstep_pipi, args.do_vacuum_subtraction,
			    args.data_dir_symm, args.traj_start, args.traj_inc, args.traj_lessthan, args.bin_size,
			    ffn, bfn_src, bfn_snk, "", opts); 
  }
 
  compareRelativeDifferences(pipi_j_asymm,pipi_j_symm);
  
  std::cout << "Done\n";
  return 0;
}

#endif
