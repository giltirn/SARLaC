#ifndef _BOOTSTRAP_RESAMPLE_TABLE_H_
#define _BOOTSTRAP_RESAMPLE_TABLE_H_

#include<distribution/bootstrap/class.h>
#include<random/utils.h>
#include<parser/parser.h>

SARLAC_START_NAMESPACE

static std::vector<std::vector<int> > nonoverlappingBlockResampleTable(RNGstore &brng, const int nsample, const int block_size,
								       const int nboots = bootstrapDistributionOptions::defaultBoots()){
  int nblock = nsample / block_size;
  int nsample_b = nblock * block_size; //just in case block size doesn't divide the number of samples equally
    
  std::uniform_int_distribution<> dis(0,nblock-1);
  std::vector<std::vector<int> > out(nboots, std::vector<int>(nsample_b));
  for(int b=0;b<nboots;b++){
    for(int i=0;i<nblock;i++){
      int block_idx = dis(brng());

      int to_off = block_size*i;
      int from_off = block_size*block_idx;

      for(int ss=0;ss<block_size;ss++)
	out[b][to_off + ss] = from_off + ss;
    }
  }
  return out;
}

inline static std::vector<std::vector<int> > nonoverlappingBlockResampleTable(const int nsample, const int block_size,
									      const int nboots = bootstrapDistributionOptions::defaultBoots(),
									      const int seed = bootstrapDistributionOptions::defaultSeed()){
  RNGstore brng(seed);
  return nonoverlappingBlockResampleTable(brng, nsample, block_size, nboots);
}

static std::vector<std::vector<int> > overlappingBlockResampleTable(RNGstore &brng, const int nsample, const int block_size,
								    const int nboots = bootstrapDistributionOptions::defaultBoots()){
  //Here we use overlapping blocks of size block_size, of which there are nsample-block_size+1. These are the blocks that are resampled
  //We draw floor(n/block_size) of these with replacement and set them down in order, which gives n samples back
  int nblock_ov = nsample-block_size+1;
  int nresample = nsample/block_size;
  int nsample_b = nresample * block_size; //just in case nsample is not evenly divisible by block size

  std::uniform_int_distribution<> dis(0,nblock_ov-1);

  std::vector<std::vector<int> > out(nboots, std::vector<int>(nsample_b));    
  for(int b=0;b<nboots;b++){
    for(int i=0;i<nresample;i++){
      int block_idx = dis(brng());

      int to_off = block_size*i;
      int from_off = block_idx;

      for(int ss=0;ss<block_size;ss++)
	out[b][to_off + ss] = from_off + ss;
    }
  }
  return out;
}

inline static std::vector<std::vector<int> > overlappingBlockResampleTable(const int nsample, const int block_size,
									   const int nboots = bootstrapDistributionOptions::defaultBoots(), 
									   const int seed = bootstrapDistributionOptions::defaultSeed()){
  RNGstore brng(seed);
  return overlappingBlockResampleTable(brng, nsample, block_size, nboots);
}


//One problem with the overlapping block bootstrap is that the underlying samples do not appear in the blocks with equal frequency. For example, sample 0 only appears in block 0, whereas sample 1 appears in 0 and 1, sample 2 in 0,1,2 and so on up to the block size. This can apparently skew the mean away from that of the underlying ensemble. In the circular version the configurations are laid out on a circle
static std::vector<std::vector<int> > circularOverlappingBlockResampleTable(RNGstore &brng, const int nsample, const int block_size,
									    const int nboots = bootstrapDistributionOptions::defaultBoots()){
  //Here we use circularly overlapping blocks of size block_size, of which there are nsample. These are the blocks that are resampled
  //We draw floor(n/block_size) of these with replacement and set them down in order, which gives n samples back
  int nblock_ov = nsample;
  int nresample = nsample/block_size;
  int nsample_b = nresample * block_size; //just in case nsample is not evenly divisible by block size

  std::uniform_int_distribution<> dis(0,nblock_ov-1);

  std::vector<std::vector<int> > out(nboots, std::vector<int>(nsample_b));    
  for(int b=0;b<nboots;b++){
    for(int i=0;i<nresample;i++){
      int block_idx = dis(brng());

      int to_off = block_size*i;
      int from_off = block_idx;

      for(int ss=0;ss<block_size;ss++)
	out[b][to_off + ss] = (from_off + ss) % nsample;
    }
  }
  return out;
}

inline static std::vector<std::vector<int> > circularOverlappingBlockResampleTable(const int nsample, const int block_size,
										   const int nboots = bootstrapDistributionOptions::defaultBoots(), 
										   const int seed = bootstrapDistributionOptions::defaultSeed()){
  RNGstore brng(seed);
  return circularOverlappingBlockResampleTable(brng, nsample, block_size, nboots);
}


static std::vector<std::vector<int> > balancedNonoverlappingBlockResampleTable(RNGstore &brng, const int nsample, const int block_size,
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
  idxB = randomPermutation(idxB, brng);

  //Generate resample table from the result
  std::vector<std::vector<int> > out(nboots, std::vector<int>(nsample_b));

  q = 0;
  for(int b=0;b<nboots;b++){
    for(int i=0;i<nblock;i++){
      int block_idx = idxB[q++];

      int to_off = block_size*i;
      int from_off = block_size*block_idx;

      for(int ss=0;ss<block_size;ss++)
	out[b][to_off + ss] = from_off + ss;
    }
  }
  return out;
}

inline static std::vector<std::vector<int> > balancedNonoverlappingBlockResampleTable(const int nsample, const int block_size,
										      const int nboots = bootstrapDistributionOptions::defaultBoots(),
										      const int seed = bootstrapDistributionOptions::defaultSeed()){
  RNGstore brng(seed);
  return balancedNonoverlappingBlockResampleTable(brng, nsample, block_size, nboots);
}


  
//Generate the mapping between the resampled ensembles and the original. Output indices are [boot][sample]
static std::vector<std::vector<int> > resampleTable(RNGstore &brng, const int nsample, 
						    const int nboots = bootstrapDistributionOptions::defaultBoots()){
  std::uniform_int_distribution<> dis(0,nsample-1);
  std::vector<std::vector<int> > out(nboots, std::vector<int>(nsample));
  for(int b=0;b<nboots;b++)
    for(int i=0;i<nsample;i++)
      out[b][i] = dis(brng());
  return out;
}

inline static std::vector<std::vector<int> > resampleTable(const int nsample, 
							   const int nboots = bootstrapDistributionOptions::defaultBoots(), 
							   const int seed = bootstrapDistributionOptions::defaultSeed()){
  RNGstore brng(seed);
  return resampleTable(brng, nsample, nboots);
}

//Threaded version
static std::vector<std::vector<int> > resampleTable(threadRNGstore &tbrng, const int nsample, 
						    const int nboots = bootstrapDistributionOptions::defaultBoots()){
  std::uniform_int_distribution<> dis(0,nsample-1);
  std::vector<std::vector<int> > out(nboots, std::vector<int>(nsample));
#pragma omp parallel for
  for(int b=0;b<nboots;b++)
    for(int i=0;i<nsample;i++)
      out[b][i] = dis(tbrng()());
  return out;
}
  



GENERATE_ENUM_AND_PARSER(BootResampleTableType, (Basic)(NonOverlappingBlock)(OverlappingBlock)(CircularOverlappingBlock)(BalancedNonOverlappingBlock) );
GENERATE_HDF5_ENUM_SERIALIZE(BootResampleTableType);

struct resampleTableOptions{
  bool read_from_file;
  std::string read_file;

  bool write_to_file;
  std::string write_file;

  resampleTableOptions(): read_from_file(false), write_to_file(false){}
};

std::vector<std::vector<int> > generateResampleTable(const size_t nsample, const size_t nboot, 
						     const BootResampleTableType table_type,
						     const size_t block_size, RNGstore &rng=RNG, const resampleTableOptions &opt = resampleTableOptions()){
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
    assert(rng.isInitialized());

    switch(table_type){
    case BootResampleTableType::Basic:
      otable = resampleTable(rng,nsample,nboot); break;
    case BootResampleTableType::NonOverlappingBlock:
      otable = nonoverlappingBlockResampleTable(rng,nsample,block_size, nboot); break;
    case BootResampleTableType::OverlappingBlock:
      otable = overlappingBlockResampleTable(rng,nsample,block_size, nboot); break;
    case BootResampleTableType::CircularOverlappingBlock:
      otable = circularOverlappingBlockResampleTable(rng,nsample,block_size, nboot); break;
    case BootResampleTableType::BalancedNonOverlappingBlock:
      otable = balancedNonoverlappingBlockResampleTable(rng,nsample,block_size, nboot); break;
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



SARLAC_END_NAMESPACE

#endif
