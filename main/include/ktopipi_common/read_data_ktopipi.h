#ifndef _FIT_KTOPIPI_READ_DATA_H
#define _FIT_KTOPIPI_READ_DATA_H

#include<config.h>
#include<utils/macros.h>

#include<algorithm>
#include<parser.h>

#include<pipi_common/read_data.h>

#include "data_containers.h"

//Parse K->pipi type 1-4 data files and A1 projected pipi bubble into raw data structures

CPSFIT_START_NAMESPACE

struct KtoPiPiFilenamePolicyGen{
  subStringReplace repl_type1; //File format expects <TRAJ> <TSEP_K_PI> <TSEP_PIPI> <MOM>   where mom is the momentum of pi_1^snk
  subStringReplace repl_type2; //File format expects <TRAJ> <TSEP_K_PI> <TSEP_PIPI>
  subStringReplace repl_type3; //Same as type 2
  subStringReplace repl_type4; //File format expects <TRAJ>
  
  KtoPiPiFilenamePolicyGen(const std::string &type1_fmt, 
			   const std::string &type2_fmt, 
			   const std::string &type3_fmt,
			   const std::string &type4_fmt):
    repl_type1(type1_fmt, { subStringSpecify("<TRAJ>"), subStringSpecify("<TSEP_K_PI>"), subStringSpecify("<TSEP_PIPI>"), subStringSpecify("<MOM>") }),
    repl_type2(type2_fmt, { subStringSpecify("<TRAJ>"), subStringSpecify("<TSEP_K_PI>"), subStringSpecify("<TSEP_PIPI>") }),
    repl_type3(type3_fmt, { subStringSpecify("<TRAJ>"), subStringSpecify("<TSEP_K_PI>"), subStringSpecify("<TSEP_PIPI>") }),
    repl_type4(type4_fmt, { subStringSpecify("<TRAJ>") })
  {}
  
  inline std::string type1filename(const std::string &data_dir, const int traj, const int tsep_k_pi, const int tsep_pipi, threeMomentum p) const{
#ifdef DAIQIAN_PION_PHASE_CONVENTION
    p = -p;
#endif
    std::ostringstream fname;
    fname << data_dir << '/';
    repl_type1.replace(fname, { anyToStr(traj), anyToStr(tsep_k_pi), anyToStr(tsep_pipi), momStr(p) });
    return fname.str();
  }
  inline std::string type2filename(const std::string &data_dir, const int traj, const int tsep_k_pi, const int tsep_pipi) const{
    std::ostringstream fname;
    fname << data_dir << '/';
    repl_type2.replace(fname, { anyToStr(traj), anyToStr(tsep_k_pi), anyToStr(tsep_pipi) });
    return fname.str();
  }
  inline std::string type3filename(const std::string &data_dir, const int traj, const int tsep_k_pi, const int tsep_pipi) const{
    std::ostringstream fname;
    fname << data_dir << '/';
    repl_type3.replace(fname, { anyToStr(traj), anyToStr(tsep_k_pi), anyToStr(tsep_pipi) });
    return fname.str();
  }
  inline std::string type4filename(const std::string &data_dir, const int traj) const{
    std::ostringstream fname;
    fname << data_dir << '/';
    repl_type4.replace(fname, { anyToStr(traj) });
    return fname.str();
  }
};


template<typename FilenamePolicy>
class BasicKtoPiPiReadPolicy{
  const std::string data_dir;
  const int traj_start;
  const int traj_inc;
  const int traj_lessthan;
  const FilenamePolicy &fp;

  int traj(const int sample) const{ return traj_start + sample*traj_inc; }

public:
  int nSample() const{ return (traj_lessthan - traj_start)/traj_inc; }

  inline std::string type1filename(const int sample, const int tsep_k_pi, const int tsep_pipi, const threeMomentum &p) const{
    return fp.type1filename(data_dir, traj(sample), tsep_k_pi, tsep_pipi, p);
  }
  inline std::string type2filename(const int sample, const int tsep_k_pi, const int tsep_pipi) const{
    return fp.type2filename(data_dir, traj(sample), tsep_k_pi, tsep_pipi);
  }
  inline std::string type3filename(const int sample, const int tsep_k_pi, const int tsep_pipi) const{
    return fp.type3filename(data_dir, traj(sample), tsep_k_pi, tsep_pipi);
  }
  inline std::string type4filename(const int sample) const{
    return fp.type4filename(data_dir, traj(sample));
  }

  BasicKtoPiPiReadPolicy(const std::string &data_dir, int traj_start, int traj_inc, int traj_lessthan, const FilenamePolicy &fp): 
    data_dir(data_dir), traj_start(traj_start), traj_inc(traj_inc), traj_lessthan(traj_lessthan), fp(fp){}
};

//Will sum over all pion mom provided with the coefficient provided to project onto the desired rotational state (only applies to type1 ; other types the average is performed at runtime)
//For original calculation use  {  { {1,1,1}, 1.0/8.0 }, { {-1,-1,-1}, 1.0/8.0 },  { {-1,1,1}, 3.0/8.0 }, { {1,-1,-1}, 3.0/8.0 }  }
template<typename ReadPolicy>
type1234Data readKtoPiPiType(const int type, const int tsep_k_pi, const int tsep_pipi, const int Lt, const std::vector<std::pair<threeMomentum, double> > &type1_pimom_proj,
			     const ReadPolicy &rp){
  int nsample = rp.nSample();
  if(type == 1){
    std::vector<type1234Data> type1_momdata(type1_pimom_proj.size(), type1234Data(1,Lt,nsample) );
#pragma omp parallel for
    for(int sample=0;sample<nsample;sample++){
      for(int pp=0;pp<type1_pimom_proj.size();pp++){
	threeMomentum p = type1_pimom_proj[pp].first;
	std::string filename = rp.type1filename(sample, tsep_k_pi, tsep_pipi, p);
	type1_momdata[pp].parse(filename, sample);
      }
    }
    type1234Data out = type1_momdata[0] * type1_pimom_proj[0].second;
    for(int pp=1;pp<type1_pimom_proj.size();pp++)
      out = out +  type1_momdata[pp] * type1_pimom_proj[pp].second;

    return out;
  }else if(type == 2 || type == 3){
    type1234Data typedata(type,Lt,nsample);
#pragma omp parallel for
    for(int sample=0;sample<nsample;sample++){
      std::string filename = type ==2 ? rp.type2filename(sample, tsep_k_pi, tsep_pipi) : rp.type3filename(sample, tsep_k_pi, tsep_pipi);
      typedata.parse(filename,sample);
    }
#ifdef DAIQIAN_COMPATIBILITY_MODE
    if(type == 2 || type == 3) typedata = typedata * 0.5; //correct for missing coefficient
#endif
    return typedata;
  }else if(type == 4){
    type1234Data typedata(type,Lt,nsample);
#pragma omp parallel for
    for(int sample=0;sample<nsample;sample++){
      std::string filename = rp.type4filename(sample);
      typedata.parse(filename,sample);
    }
    return typedata;
  }else{
    error_exit(std::cout << "readKtoPiPiType invalid type " << type << std::endl);
  }
}

template<typename FilenamePolicy>
inline type1234Data readKtoPiPiType(const int type, const int traj_start, const int traj_inc, const int traj_lessthan,
				    const int tsep_k_pi, const int tsep_pipi, const int Lt, const std::vector<std::pair<threeMomentum, double> > &type1_pimom_proj,
				    const std::string &data_dir, const FilenamePolicy &fp){
  BasicKtoPiPiReadPolicy<FilenamePolicy> rp(data_dir, traj_start, traj_inc, traj_lessthan, fp);
  return readKtoPiPiType(type, tsep_k_pi, tsep_pipi, Lt, type1_pimom_proj, rp);
}


//For bubble_pimom_proj in original calculation use  {  { {1,1,1}, 1.0/8.0 }, { {-1,-1,-1}, 1.0/8.0 },  
//                                                      { {-1,1,1}, 1.0/8.0 }, { {1,-1,-1}, 1.0/8.0 }, 
//                                                      { {1,-1,1}, 1.0/8.0 }, { {-1,1,-1}, 1.0/8.0 }, 
//                                                      { {1,1,-1}, 1.0/8.0 }, { {-1,-1,1}, 1.0/8.0 } }
NumericTensor<rawDataDistributionD,1> readProjectedPiPiBubble(const std::string &data_dir, const std::string &file_fmt,
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
