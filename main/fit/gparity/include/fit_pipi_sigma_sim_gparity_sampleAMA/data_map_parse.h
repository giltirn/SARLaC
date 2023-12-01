#ifndef _GLOBAL_DATA_MAP_PARSE_H
#define _GLOBAL_DATA_MAP_PARSE_H

#include<config.h>
#include<utils/macros.h>
#include<tensors/numeric_tensor.h>
#include<common.h>

#include <pipi_common/pipi_common.h>
#include "data_map.h"

SARLAC_START_NAMESPACE

//Example parser
//File format  <t> <re> <im>
//Expect file name format including substring <CONF> which is replaced by the config index
struct StandardParseReal{
  typedef NumericTensor<rawDataDistributionD,2> OutputType;
  int Lt;
  subStringReplace fmt;

  StandardParseReal(const std::string &format, int Lt): Lt(Lt){
    static std::vector<subStringSpecify> spec = { subStringSpecify("<CONF>") };
    fmt.chunkString(format, spec);
  }
  
  OutputType initialize(const int nsample) const{
    return OutputType({Lt,Lt}, [&](int const* tt){ return rawDataDistributionD(nsample); });
  }
  
  void parse(OutputType &out, const int sample, const std::string &directory, const int config) const{
    std::ostringstream fname;
    fname << directory << '/';
    fmt.replace(fname, {anyToStr(config)});
    
    if(!fileExists(fname.str())) error_exit(std::cout << "File \"" << fname.str() << "\" does not exist!\n");

    std::cout << "Parsing " << fname.str() << std::endl;

    std::ifstream f(fname.str());
    for(int tsrc=0;tsrc<Lt;tsrc++){
      for(int tsep=0;tsep<Lt;tsep++){	
	assert(f.good());
	int tsrc_read, tsep_read;
	double re, im;
	f >> tsrc_read >> tsep_read;
	assert(tsrc_read == tsrc && tsep_read == tsep);
	f >> re >> im;
	out({tsrc,tsep}).sample(sample) = re;
      }
    }
    assert(!f.fail() && !f.bad());
  }
};
	
template<typename FileParser>
typename FileParser::OutputType readData(const std::set<int> &outer_configs, const std::string &tag, const DataMap &dmap, const FileParser &parser){
  std::vector<std::pair<std::string, int> > toread = getDirectoryConfigList(outer_configs, tag, dmap);
  typename FileParser::OutputType out = parser.initialize(outer_configs.size());
#pragma omp parallel for
  for(int i=0;i<toread.size();i++){
    parser.parse(out, i, toread[i].first, toread[i].second);
  }
  return out;
}

//Parse pipi->sigma data and project onto s-wave
struct PiPiToSigmaFileParseProject{
  typedef figureData OutputType;
  subStringReplace fmt;
  int Lt;

  PiPiToSigmaFileParseProject(const std::string &format, int Lt): Lt(Lt){
    static std::vector<subStringSpecify> spec = { subStringSpecify("<CONF>"), subStringSpecify("<MOM_QUARK_SIGMA>"), subStringSpecify("<MOM_PI>") };
    fmt.chunkString(format, spec);
  }
  
  OutputType initialize(const int nsample) const{
    figureData out; 
    out.setup(Lt,rawDataDistributionD(nsample));
    return out;
  }
  
  void parse(OutputType &out, const int sample, const std::string &directory, const int config) const{
    static std::vector<threeMomentum> quark_mom = { {1,1,1}, {-1,-1,-1},
						    {-3,1,1}, {3,-1,-1},
						    {1,-3,1}, {-1,3,-1},
						    {1,1,-3}, {-1,-1,3} };
    
  
    static std::vector<threeMomentum> pion_mom = { {2,2,2}, {-2,-2,-2},
						   {-2,2,2}, {2,-2,-2},
						   {2,-2,2}, {-2,2,-2},
						   {2,2,-2}, {-2,-2,2} };

    figureData accum_raw_data(Lt,rawDataDistributionD(1)), tmp_raw_data(Lt,rawDataDistributionD(1)); 
    accum_raw_data.zero();

    for(int ppiidx = 0 ; ppiidx < 8 ; ppiidx++){
      for(int psigqidx = 0 ; psigqidx < 8 ; psigqidx++){
	std::ostringstream filename;
	filename << directory << '/';
	fmt.replace(filename, {anyToStr(config), momStr(quark_mom[psigqidx]), momStr(pion_mom[ppiidx])} );
	std::cout << "Parsing " << filename.str() << std::endl;
	tmp_raw_data.parseCDR(filename.str(), 0);
	accum_raw_data = accum_raw_data + tmp_raw_data;
      }
    }
    
    for(int t1=0;t1<Lt;t1++)
      for(int t2=0;t2<Lt;t2++)
	out(t1,t2).sample(sample) = accum_raw_data(t1,t2).sample(0)/64.;
  }
};


SARLAC_END_NAMESPACE

#endif
