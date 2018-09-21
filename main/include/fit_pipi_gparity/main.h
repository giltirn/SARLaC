#ifndef _PIPI_MAIN_H__
#define _PIPI_MAIN_H__

#include <pipi_common/pipi_common.h>

using namespace CPSfit;

#include "args.h"
#include "cmdline.h"

//User provides a set of total source pipi momenta. The data are read, rotational-state projected, resampled then averaged over the provided set of total momenta. The resulting correlation function is returned
doubleJackCorrelationFunction generateData(const Args &args, const CMDline &cmdline){
  getResampledPiPi2ptDataOpts opts;
#define CP(A) opts.A = cmdline.A
  CP(load_hdf5_data_checkpoint);
  CP(load_hdf5_data_checkpoint_stub);
  CP(save_hdf5_data_checkpoint);
  CP(save_hdf5_data_checkpoint_stub);
#undef CP

  doubleJackCorrelationFunction out;
  for(int p=0;p<args.total_mom.size();p++){
    std::unique_ptr<PiPiProjectorBase> proj_src( getProjector(args.proj_src, args.total_mom[p]) );
    std::unique_ptr<PiPiProjectorBase> proj_snk( getProjector(args.proj_snk, -args.total_mom[p]) );

    PiPiSymmetrySubsetFigureFileMapping ffn(args.data_dir, args.figure_file_format, args.traj_start, args.tsep_pipi, getSrcSnkMomentumSet(*proj_src, *proj_snk), args.total_mom[p]);

    bubbleFilenamePolicyGeneric bfn_src(args.bubble_file_format, args.total_mom[p], Source);
    bubbleFilenamePolicyGeneric bfn_snk(args.bubble_file_format, args.total_mom[p], Sink);

    doubleJackCorrelationFunction outp = getResampledPiPi2ptData(*proj_src,*proj_snk, args.isospin, args.Lt, args.tsep_pipi, args.tstep_pipi, args.do_vacuum_subtraction,
								 args.data_dir, args.traj_start, args.traj_inc, args.traj_lessthan, args.bin_size,
								 ffn, bfn_src, bfn_snk, std::string("ptot") + momStr(args.total_mom[p]), opts); 
    if(p==0) out = std::move(outp);
    else out = out + outp;
  }
  if(args.total_mom.size() > 1) out = out * (1./args.total_mom.size());
  return out;
}

//Load a checkpoint of precomputed resampled correlation function or create it from raw data
doubleJackCorrelationFunction getData(const Args &args, const CMDline &cmdline){
  doubleJackCorrelationFunction data;
  if(cmdline.load_combined_data){
#ifdef HAVE_HDF5
    HDF5reader reader(cmdline.load_combined_data_file);
    read(reader, data, "data");
#else
    error_exit("getData: Reading amplitude data requires HDF5\n");
#endif
  }else{
    data = generateData(args,cmdline);
  }

  if(cmdline.save_combined_data){
#ifdef HAVE_HDF5
    HDF5writer writer(cmdline.save_combined_data_file);
    write(writer, data, "data");
#else
    error_exit("getData: Saving amplitude data requires HDF5\n");
#endif
  }
  return std::move(data);
}



#endif
