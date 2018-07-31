#ifndef _FIT_KTOPIPI_READ_DATA_H
#define _FIT_KTOPIPI_READ_DATA_H

#include<config.h>
#include<utils/macros.h>

#include<algorithm>
#include<parser.h>

#include<pipi_common/read_data.h>

#include "data_containers.h"

CPSFIT_START_NAMESPACE

//Will sum over all pion mom provided with the coefficient provided to project onto the desired rotational state (only applies to type1 ; other types the average is performed at runtime)
//For original calculation use  {  { {1,1,1}, 1.0/8.0 }, { {-1,-1,-1}, 1.0/8.0 },  { {-1,1,1}, 3.0/8.0 }, { {1,-1,-1}, 3.0/8.0 }  }
type1234Data readType(const int type, const int traj_start, const int traj_inc, const int traj_lessthan,
		      const int tsep_k_pi, const int tsep_pipi, const int Lt, const std::vector<std::pair<threeMomentum, double> > &type1_pimom_proj,
		      const std::string &data_dir, const std::string &file_fmt){
  int nsample = (traj_lessthan - traj_start)/traj_inc;
  if(type == 1){
    //File format expects <TRAJ> <TSEP_K_PI> <TSEP_PIPI> <MOM>   where mom is the momentum of pi_1^snk
    static std::vector<subStringSpecify> keys { subStringSpecify("<TRAJ>"), subStringSpecify("<TSEP_K_PI>"), subStringSpecify("<TSEP_PIPI>"), subStringSpecify("<MOM>") };
    subStringReplace repl(file_fmt, keys);
    
    std::vector<type1234Data> type1_momdata(type1_pimom_proj.size(), type1234Data(1,Lt,nsample) );
#pragma omp parallel for
    for(int i=0;i<nsample;i++){
      const int traj = traj_start + i*traj_inc;
      for(int pp=0;pp<type1_pimom_proj.size();pp++){
	threeMomentum p = type1_pimom_proj[pp].first;
#ifdef DAIQIAN_PION_PHASE_CONVENTION
	p = -p;
#endif
	std::ostringstream fname;
	fname << data_dir << '/';
	repl.replace(fname, { anyToStr(traj), anyToStr(tsep_k_pi), anyToStr(tsep_pipi), momStr(p) });
	type1_momdata[pp].parse(fname.str(), i);
      }
    }
    type1234Data out = type1_momdata[0] * type1_pimom_proj[0].second;
    for(int pp=1;pp<type1_pimom_proj.size();pp++)
      out = out +  type1_momdata[pp] * type1_pimom_proj[pp].second;

    return out;
  }else if(type == 2 || type == 3){
    //File format expects <TRAJ> <TSEP_K_PI> <TSEP_PIPI>
    static std::vector<subStringSpecify> keys { subStringSpecify("<TRAJ>"), subStringSpecify("<TSEP_K_PI>"), subStringSpecify("<TSEP_PIPI>") };
    subStringReplace repl(file_fmt, keys);

    type1234Data typedata(type,Lt,nsample);
#pragma omp parallel for
    for(int i=0;i<nsample;i++){
      const int traj = traj_start + i*traj_inc;

      std::ostringstream fname;
      fname << data_dir << '/';
      repl.replace(fname, { anyToStr(traj), anyToStr(tsep_k_pi), anyToStr(tsep_pipi) });
      typedata.parse(fname.str(),i);
    }
#ifdef DAIQIAN_COMPATIBILITY_MODE
    if(type == 2 || type == 3) typedata = typedata * 0.5; //correct for missing coefficient
#endif
    return typedata;
  }else if(type == 4){
    //File format expects <TRAJ>
    static std::vector<subStringSpecify> keys { subStringSpecify("<TRAJ>") };
    subStringReplace repl(file_fmt, keys);

    type1234Data typedata(type,Lt,nsample);
#pragma omp parallel for
    for(int i=0;i<nsample;i++){
      const int traj = traj_start + i*traj_inc;

      std::ostringstream fname;
      fname << data_dir << '/';
      repl.replace(fname, { anyToStr(traj) });
      typedata.parse(fname.str(),i);
    }
    return typedata;
  }else{
    error_exit(std::cout << "readType invalid type " << type << std::endl);
  }
}

//For bubble_pimom_proj in original calculation use  {  { {1,1,1}, 1.0/8.0 }, { {-1,-1,-1}, 1.0/8.0 },  
//                                                      { {-1,1,1}, 1.0/8.0 }, { {1,-1,-1}, 1.0/8.0 }, 
//                                                      { {1,-1,1}, 1.0/8.0 }, { {-1,1,-1}, 1.0/8.0 }, 
//                                                      { {1,1,-1}, 1.0/8.0 }, { {-1,-1,1}, 1.0/8.0 } }
NumericTensor<rawDataDistributionD,1> readProjectedBubble(const std::string &data_dir, const std::string &file_fmt,
							  const int traj_start, const int traj_inc, const int traj_lessthan, 
							  const int Lt, const int tsep_pipi, 
							  const std::vector<std::pair<threeMomentum, double> > &bubble_pimom_proj){
  bubbleDataAllMomenta raw_bubble_data;

  int nsample = (traj_lessthan - traj_start)/traj_inc;
  raw_bubble_data.setup(Lt,tsep_pipi,nsample);
  
  std::vector<threeMomentum> pion_momenta(bubble_pimom_proj.size());
  for(int i=0;i<bubble_pimom_proj.size();i++) 
    pion_momenta[i] = bubble_pimom_proj[i].first;

  bubbleFilenamePolicyGeneric fn(file_fmt, {0,0,0}, Sink);
  readBubble(raw_bubble_data, data_dir, tsep_pipi, Lt, traj_start, traj_inc, traj_lessthan, fn, pion_momenta, Sink);

  NumericTensor<rawDataDistributionD,1> out({Lt}, rawDataDistributionD(nsample,0.));
  for(int t=0;t<Lt;t++)
    for(int p=0;p<bubble_pimom_proj.size();p++)
      out({t}) = out({t}) + raw_bubble_data(Sink,bubble_pimom_proj[p].first)(t) * bubble_pimom_proj[p].second;

  return out;
}


#ifdef HAVE_HDF5
void writeProjectedBubble(HDF5writer &writer, const NumericTensor<rawDataDistributionD,1> &value, const std::string &tag){
  boost::timer::auto_cpu_timer total_time("writeBubble(HDF5writer &, const NumericTensor<rawDataDistributionD,1> &, const std::string &)  %w s\n");
  writer.enter(tag);
  int size = value.size(0);
  int nsample = value({0}).size();
  
  write(writer, size, "size");
  write(writer, nsample, "nsample");

  std::vector<double> data(size * nsample);
  for(int i=0;i<size;i++)
    for(int s=0;s<nsample;s++)
      data[s + nsample*i] = value(&i).sample(s);
  write(writer, data, "data");
  writer.leave();
}
void readProjectedBubble(HDF5reader &reader, NumericTensor<rawDataDistributionD,1> &value, const std::string &tag){
  boost::timer::auto_cpu_timer total_time("readBubble(HDF5writer &, NumericTensor<rawDataDistributionD,1> &, const std::string &)  %w s\n");
  reader.enter(tag);
  int size;
  int nsample;
  
  read(reader, size, "size");
  read(reader, nsample, "nsample");

  std::vector<double> data(size * nsample);
  read(reader, data, "data");

  value.resize(&size, rawDataDistributionD(nsample));
  
  for(int i=0;i<size;i++)
    for(int s=0;s<nsample;s++)
      value(&i).sample(s) = data[s + nsample*i];

  reader.leave();
}
#endif


CPSFIT_END_NAMESPACE

#endif
