#ifndef _BOOTSTRAP_RESAMPLE_TABLE_H_
#define _BOOTSTRAP_RESAMPLE_TABLE_H_

#include<distribution/bootstrap/class.h>
#include<distribution/raw_data_distribution/class.h>
#include<random/utils.h>
#include<parser/parser.h>

SARLAC_START_NAMESPACE

template<typename RNGstoreType>
struct RNGwrap{};
template<>
struct RNGwrap<RNGstore>{
  RNGstore &rng;
  RNGwrap(RNGstore &rng): rng(rng){
    assert(rng.isInitialized());
  }

  RNGstore &getRNG(const int thr) const{ return rng; }

  template<typename BootOp>
  inline void bootLoop(int nboots, const BootOp &op) const{
    for(int b=0;b<nboots;b++)
      op(b, rng);
  }
};
template<>
struct RNGwrap<threadRNGstore>{
  threadRNGstore &rng;
  RNGwrap(threadRNGstore &rng): rng(rng){
    for(int t=0;t<rng.size();t++) assert(rng(t).isInitialized());
  }

  RNGstore &getRNG(const int thr) const{ return rng(thr); }

  template<typename BootOp>
  inline void bootLoop(int nboots, const BootOp &op) const{
#pragma omp parallel for schedule(static)
    for(int b=0;b<nboots;b++)
      op(b, rng());
  }
};

///////////////////////////////////////////////////////////////////////////////////
//Non-overlapping block bootstrap
///////////////////////////////////////////////////////////////////////////////////

template<typename RNGwrapType>
std::vector<std::vector<int> > nonoverlappingBlockResampleTable(const RNGwrapType &wrng, const int nsample, const int block_size,
								const int nboots = bootstrapDistributionOptions::defaultBoots()){
  int nblock = nsample / block_size;
  int nsample_b = nblock * block_size; //just in case block size doesn't divide the number of samples equally
    
  std::uniform_int_distribution<> dis(0,nblock-1);
  std::vector<std::vector<int> > out(nboots, std::vector<int>(nsample_b));

  wrng.bootLoop(nboots, [&](const int b, RNGstore &brng){
    for(int i=0;i<nblock;i++){
      int block_idx = dis(brng());

      int to_off = block_size*i;
      int from_off = block_size*block_idx;

      for(int ss=0;ss<block_size;ss++)
	out[b][to_off + ss] = from_off + ss;
    }
    });

  return out;
}

inline std::vector<std::vector<int> > nonoverlappingBlockResampleTable(RNGstore &brng, const int nsample, const int block_size,
								       const int nboots = bootstrapDistributionOptions::defaultBoots()){
  return nonoverlappingBlockResampleTable(RNGwrap<RNGstore>(brng), nsample, block_size, nboots);
}
inline std::vector<std::vector<int> > nonoverlappingBlockResampleTable(threadRNGstore &brng, const int nsample, const int block_size,
								       const int nboots = bootstrapDistributionOptions::defaultBoots()){
  return nonoverlappingBlockResampleTable(RNGwrap<threadRNGstore>(brng), nsample, block_size, nboots);
}

inline static std::vector<std::vector<int> > nonoverlappingBlockResampleTable(const int nsample, const int block_size,
									      const int nboots = bootstrapDistributionOptions::defaultBoots(),
									      const int seed = bootstrapDistributionOptions::defaultSeed()){
  RNGstore brng(seed);
  return nonoverlappingBlockResampleTable(brng, nsample, block_size, nboots);
}

///////////////////////////////////////////////////////////////////////////////////
//Overlapping block bootstrap
///////////////////////////////////////////////////////////////////////////////////


template<typename RNGwrapType>
std::vector<std::vector<int> > overlappingBlockResampleTable(const RNGwrapType &wrng, const int nsample, const int block_size,
								    const int nboots = bootstrapDistributionOptions::defaultBoots()){
  //Here we use overlapping blocks of size block_size, of which there are nsample-block_size+1. These are the blocks that are resampled
  //We draw floor(n/block_size) of these with replacement and set them down in order, which gives n samples back
  int nblock_ov = nsample-block_size+1;
  int nresample = nsample/block_size;
  int nsample_b = nresample * block_size; //just in case nsample is not evenly divisible by block size

  std::uniform_int_distribution<> dis(0,nblock_ov-1);

  std::vector<std::vector<int> > out(nboots, std::vector<int>(nsample_b));
  wrng.bootLoop(nboots, [&](const int b, RNGstore &brng){
    for(int i=0;i<nresample;i++){
      int block_idx = dis(brng());

      int to_off = block_size*i;
      int from_off = block_idx;

      for(int ss=0;ss<block_size;ss++)
	out[b][to_off + ss] = from_off + ss;
    }
    });
  return out;
}

inline std::vector<std::vector<int> > overlappingBlockResampleTable(RNGstore &brng, const int nsample, const int block_size,
								    const int nboots = bootstrapDistributionOptions::defaultBoots()){
  return overlappingBlockResampleTable(RNGwrap<RNGstore>(brng), nsample, block_size, nboots);
}
inline std::vector<std::vector<int> > overlappingBlockResampleTable(threadRNGstore &brng, const int nsample, const int block_size,
								    const int nboots = bootstrapDistributionOptions::defaultBoots()){
  return overlappingBlockResampleTable(RNGwrap<threadRNGstore>(brng), nsample, block_size, nboots);
}

inline static std::vector<std::vector<int> > overlappingBlockResampleTable(const int nsample, const int block_size,
									   const int nboots = bootstrapDistributionOptions::defaultBoots(), 
									   const int seed = bootstrapDistributionOptions::defaultSeed()){
  RNGstore brng(seed);
  return overlappingBlockResampleTable(brng, nsample, block_size, nboots);
}


///////////////////////////////////////////////////////////////////////////////////
//Circular overlapping block bootstrap
///////////////////////////////////////////////////////////////////////////////////

//One problem with the overlapping block bootstrap is that the underlying samples do not appear in the blocks with equal frequency. For example, sample 0 only appears in block 0, whereas sample 1 appears in 0 and 1, sample 2 in 0,1,2 and so on up to the block size. This can apparently skew the mean away from that of the underlying ensemble. In the circular version the configurations are laid out on a circle
template<typename RNGwrapType>
std::vector<std::vector<int> > circularOverlappingBlockResampleTable(const RNGwrapType &wrng, const int nsample, const int block_size,
								     const int nboots = bootstrapDistributionOptions::defaultBoots()){
  //Here we use circularly overlapping blocks of size block_size, of which there are nsample. These are the blocks that are resampled
  //We draw floor(n/block_size) of these with replacement and set them down in order, which gives n samples back
  int nblock_ov = nsample;
  int nresample = nsample/block_size;
  int nsample_b = nresample * block_size; //just in case nsample is not evenly divisible by block size

  std::uniform_int_distribution<> dis(0,nblock_ov-1);

  std::vector<std::vector<int> > out(nboots, std::vector<int>(nsample_b));
  wrng.bootLoop(nboots, [&](const int b, RNGstore &brng){
    for(int i=0;i<nresample;i++){
      int block_idx = dis(brng());

      int to_off = block_size*i;
      int from_off = block_idx;

      for(int ss=0;ss<block_size;ss++)
	out[b][to_off + ss] = (from_off + ss) % nsample;
    }
    });
  return out;
}

inline std::vector<std::vector<int> > circularOverlappingBlockResampleTable(RNGstore &brng, const int nsample, const int block_size,
									    const int nboots = bootstrapDistributionOptions::defaultBoots()){
   return circularOverlappingBlockResampleTable(RNGwrap<RNGstore>(brng), nsample, block_size, nboots);
}
inline std::vector<std::vector<int> > circularOverlappingBlockResampleTable(threadRNGstore &brng, const int nsample, const int block_size,
								    const int nboots = bootstrapDistributionOptions::defaultBoots()){
  return circularOverlappingBlockResampleTable(RNGwrap<threadRNGstore>(brng), nsample, block_size, nboots);
}



inline static std::vector<std::vector<int> > circularOverlappingBlockResampleTable(const int nsample, const int block_size,
										   const int nboots = bootstrapDistributionOptions::defaultBoots(), 
										   const int seed = bootstrapDistributionOptions::defaultSeed()){
  RNGstore brng(seed);
  return circularOverlappingBlockResampleTable(brng, nsample, block_size, nboots);
}

///////////////////////////////////////////////////////////////////////////////////
//Balanced non-overlapping block bootstrap
///////////////////////////////////////////////////////////////////////////////////

template<typename RNGwrapType>
std::vector<std::vector<int> > balancedNonoverlappingBlockResampleTable(const RNGwrapType &wrng, const int nsample, const int block_size,
									       const int nboots = bootstrapDistributionOptions::defaultBoots()){
  int nblock = nsample / block_size;
  int nsample_b = nblock * block_size; //just in case block size doesn't divide the number of samples equally
    
  //Generate B copies of the indices 0...nblock-1
  std::vector<int> idxB( nboots * nblock );
  int q=0;
  for(int b=0;b<nboots;b++)
    for(int i=0;i<nblock;i++)
      idxB[q++] = i;

  //Randomly permute
  idxB = randomPermutation(idxB, wrng.getRNG(0));

  //Generate resample table from the result
  std::vector<std::vector<int> > out(nboots, std::vector<int>(nsample_b));
  
  wrng.bootLoop(nboots, [&](const int b, RNGstore &brng){
    for(int i=0;i<nblock;i++){
      int block_idx = idxB[i+nblock*b];

      int to_off = block_size*i;
      int from_off = block_size*block_idx;

      for(int ss=0;ss<block_size;ss++)
	out[b][to_off + ss] = from_off + ss;
    }
    });
  return out;
}

inline std::vector<std::vector<int> > balancedNonoverlappingBlockResampleTable(RNGstore &brng, const int nsample, const int block_size,
									    const int nboots = bootstrapDistributionOptions::defaultBoots()){
   return balancedNonoverlappingBlockResampleTable(RNGwrap<RNGstore>(brng), nsample, block_size, nboots);
}
inline std::vector<std::vector<int> > balancedNonoverlappingBlockResampleTable(threadRNGstore &brng, const int nsample, const int block_size,
								    const int nboots = bootstrapDistributionOptions::defaultBoots()){
  return balancedNonoverlappingBlockResampleTable(RNGwrap<threadRNGstore>(brng), nsample, block_size, nboots);
}

inline static std::vector<std::vector<int> > balancedNonoverlappingBlockResampleTable(const int nsample, const int block_size,
										      const int nboots = bootstrapDistributionOptions::defaultBoots(),
										      const int seed = bootstrapDistributionOptions::defaultSeed()){
  RNGstore brng(seed);
  return balancedNonoverlappingBlockResampleTable(brng, nsample, block_size, nboots);
}

///////////////////////////////////////////////////////////////////////////////////
//Original
///////////////////////////////////////////////////////////////////////////////////
  
//Generate the mapping between the resampled ensembles and the original. Output indices are [boot][sample]
template<typename RNGwrapType>
std::vector<std::vector<int> > resampleTable(const RNGwrapType &wrng, const int nsample, 
					     const int nboots = bootstrapDistributionOptions::defaultBoots()){
  std::uniform_int_distribution<> dis(0,nsample-1);
  std::vector<std::vector<int> > out(nboots, std::vector<int>(nsample));
  wrng.bootLoop(nboots, [&](const int b, RNGstore &brng){
    for(int i=0;i<nsample;i++)
      out[b][i] = dis(brng());
    });
  return out;
}

inline std::vector<std::vector<int> > resampleTable(RNGstore &brng, const int nsample, 
						    const int nboots = bootstrapDistributionOptions::defaultBoots()){
  return resampleTable(RNGwrap<RNGstore>(brng), nsample, nboots);
}
inline std::vector<std::vector<int> > resampleTable(threadRNGstore &brng, const int nsample, 
						    const int nboots = bootstrapDistributionOptions::defaultBoots()){
  return resampleTable(RNGwrap<threadRNGstore>(brng), nsample, nboots);
}

inline static std::vector<std::vector<int> > resampleTable(const int nsample, 
							   const int nboots = bootstrapDistributionOptions::defaultBoots(), 
							   const int seed = bootstrapDistributionOptions::defaultSeed()){
  RNGstore brng(seed);
  return resampleTable(brng, nsample, nboots);
}


  
///////////////////////////////////////////////////////////////////////////////////
//Wrapper function with enum
///////////////////////////////////////////////////////////////////////////////////


GENERATE_ENUM_AND_PARSER(BootResampleTableType, (Basic)(NonOverlappingBlock)(OverlappingBlock)(CircularOverlappingBlock)(BalancedNonOverlappingBlock) );
GENERATE_HDF5_ENUM_SERIALIZE(BootResampleTableType);

struct resampleTableOptions{
  bool read_from_file;
  std::string read_file;

  bool write_to_file;
  std::string write_file;

  resampleTableOptions(): read_from_file(false), write_to_file(false){}
};

template<typename RNGwrapType>
std::vector<std::vector<int> > generateResampleTable(const size_t nsample, const size_t nboot, 
						     const BootResampleTableType table_type,
						     const size_t block_size, const RNGwrapType &wrng, const resampleTableOptions &opt = resampleTableOptions()){
  std::vector<std::vector<int> > otable;  //[b][s]

  if(opt.read_from_file){ //overrides table_type
    std::cout << "Reading resample table from " << opt.read_file << std::endl;

    HDF5reader rd(opt.read_file);
    BootResampleTableType tt;
    read(rd, tt, "table_type");
    
    if(tt != table_type) std::cout << "Warning: Resample table type of input file differs from requested type" << std::endl;
    
    read(rd, otable, "resample_table");
  }else{
    std::cout << "Generating resample table" << std::endl;

    switch(table_type){
    case BootResampleTableType::Basic:
      otable = std::move(resampleTable(wrng,nsample,nboot)); break;
    case BootResampleTableType::NonOverlappingBlock:
      otable = std::move(nonoverlappingBlockResampleTable(wrng,nsample,block_size, nboot)); break;
    case BootResampleTableType::OverlappingBlock:
      otable = std::move(overlappingBlockResampleTable(wrng,nsample,block_size, nboot)); break;
    case BootResampleTableType::CircularOverlappingBlock:
      otable = std::move(circularOverlappingBlockResampleTable(wrng,nsample,block_size, nboot)); break;
    case BootResampleTableType::BalancedNonOverlappingBlock:
      otable = std::move(balancedNonoverlappingBlockResampleTable(wrng,nsample,block_size, nboot)); break;
    default:
      assert(0);
    }
  }
  
  if(otable[0].size() != nsample)
    std::cout << "Samples " << nsample << " truncated to " << otable[0].size() << " due to blocking" << std::endl;

  if(opt.write_to_file){
    std::cout << "Writing resample table to " << opt.write_file << std::endl;
    HDF5writer wr(opt.write_file);
    write(wr, table_type, "table_type");
    write(wr, otable, "resample_table");
  }
  return otable;
}

inline std::vector<std::vector<int> > generateResampleTable(const size_t nsample, const size_t nboot, 
						     const BootResampleTableType table_type,
						     const size_t block_size, RNGstore &rng=RNG, const resampleTableOptions &opt = resampleTableOptions()){
  return generateResampleTable(nsample, nboot, table_type, block_size, RNGwrap<RNGstore>(rng), opt);
}
inline std::vector<std::vector<int> > generateResampleTable(const size_t nsample, const size_t nboot, 
						     const BootResampleTableType table_type,
						     const size_t block_size, threadRNGstore &rng=threadRNG, const resampleTableOptions &opt = resampleTableOptions()){
  return generateResampleTable(nsample, nboot, table_type, block_size, RNGwrap<threadRNGstore>(rng), opt);
}

//Generate resampled ensemble 'ens_idx' from the input ensemble 'in' and the resample table 'rtable'
template<typename T>
inline rawDataDistribution<T> resampledEnsemble(const rawDataDistribution<T> &in, const int ens_idx, const std::vector<std::vector<int> > &rtable){
  int nsample_reduced = rtable[0].size();
  assert(nsample_reduced <= in.size());
  rawDataDistribution<T> out(nsample_reduced);
  for(int s=0;s<nsample_reduced;s++)
    out.sample(s) = in.sample(rtable[ens_idx][s]);
  return out;
}


SARLAC_END_NAMESPACE

#endif
