#ifndef _SAMPLEAMA_RESAMPLE_H
#define _SAMPLEAMA_RESAMPLE_H

#include<config.h>
#include<utils/macros.h>

#include<distribution/raw_data_distribution.h>
#include<distribution/jackknife.h>
#include<distribution/double_jackknife.h>

CPSFIT_START_NAMESPACE


//Generalized version of the superjackknife resampling procedure supporting an arbitrary number of ensembles
template<typename T> 
jackknifeDistribution<T> superJackknifeResampleGen(const rawDataDistribution<T> &data, const int ens, const std::vector<int> sizes){
  assert(data.size() == sizes[ens]);
  int ntot = 0; for(int i=0;i<sizes.size();i++) ntot += sizes[i];
  int npre = 0; for(int i=0;i<ens;i++) npre += sizes[i];

  const int nens = sizes[ens];

  const T Smean = data.mean();
  const T Ssum = nens*Smean;
  const double nrm = 1./(nens - 1);

  jackknifeDistribution<double> out(ntot, Smean);
  for(int s=0;s<nens;s++) out.sample(npre + s) = nrm * ( Ssum - data.sample(s) );

  return out;
}

template<typename T> 
doubleJackknifeDistribution<T> superDoubleJackknifeResampleGen(const rawDataDistribution<T> &data, const int ens, const std::vector<int> sizes){
  assert(data.size() == sizes[ens]);
  int ntot = 0; for(int i=0;i<sizes.size();i++) ntot += sizes[i];
  int npre = 0; for(int i=0;i<ens;i++) npre += sizes[i];
  const int nens = sizes[ens];

  const T Smean = data.mean();
  const T Ssum = nens*Smean;
  const double invNm2 = 1./(nens - 2);
  const double invNm1 = 1./(nens - 1);

  doubleJackknifeDistribution<double> out(ntot, Smean);

  for(int i=0;i<ntot;i++){
    bool i_in_range = i>=npre && i<npre + nens;
    int jj = 0;
    for(int j=0;j<ntot;j++){
      if(j==i) continue;
     
      bool j_in_range = j>=npre && j<npre + nens;

      if(i_in_range){ //i is in range 	
	if(j_in_range){ //If i and j in range for ensemble require the double-jack value
	  out.sample(i).sample(jj) = invNm2 * ( Ssum - data.sample(i-npre) - data.sample(j-npre) );
	}else{ //if j not in range then use single-jack value
	  out.sample(i).sample(jj) = invNm1 * ( Ssum - data.sample(i-npre) );
	}
      }else{ //i is not in range
	if(j_in_range){ //If j is in range use single-jack value
	  out.sample(i).sample(jj) = invNm1 * ( Ssum - data.sample(j-npre) );
	}
	//if i and j not in range then use mean (default)
      }      
      ++jj;
    }
  }
  return out;
}


//We consider a quantity Y that resides on nS samples, and a quantity Z that resides on a disjoint set of nC samples
//We define a super-ensemble of size nS+nC for which Y resides only on the first nS and Z only on the last nC
//The algorithm for jackknife and double-jackknife then proceeds as normal by defining reduced and doubly-reduced ensembles and averaging over the remaining data
template<typename T>
inline jackknifeDistribution<T> superJackknifeResampleS(const rawDataDistribution<T> &data, const int nS, const int nC){
  return superJackknifeResampleGen(data,0,{nS,nC});
}

template<typename T>
inline jackknifeDistribution<T> superJackknifeResampleC(const rawDataDistribution<T> &data, const int nS, const int nC){
  return superJackknifeResampleGen(data,1,{nS,nC});
}

template<typename T>
inline doubleJackknifeDistribution<T> superDoubleJackknifeResampleS(const rawDataDistribution<T> &data, const int nS, const int nC){
  return superDoubleJackknifeResampleGen(data, 0, {nS,nC});
}

template<typename T>
inline doubleJackknifeDistribution<T> superDoubleJackknifeResampleC(const rawDataDistribution<T> &data, const int nS, const int nC){
  return superDoubleJackknifeResampleGen(data, 1, {nS,nC});
}

template<typename DistributionType>
struct sampleAMAresample{};

template<>
struct sampleAMAresample<jackknifeDistribution<double> >{
  static inline jackknifeDistribution<double> resample(const rawDataDistribution<double> &in, const char ens, const int nS, const int nC){
    return ens == 'S' ? superJackknifeResampleS(in,nS,nC) : superJackknifeResampleC(in,nS,nC);
  }
};
template<>
struct sampleAMAresample<doubleJackknifeDistribution<double> >{
  static inline doubleJackknifeDistribution<double> resample(const rawDataDistribution<double> &in, const char ens, const int nS, const int nC){
    return ens == 'S' ? superDoubleJackknifeResampleS(in,nS,nC) : superDoubleJackknifeResampleC(in,nS,nC);
  }
};

//Class version of the above that can be used to resample many data
class sampleAMA_resampler{
  char ens;
  int nS;
  int nC;
public:
  sampleAMA_resampler(){}
  sampleAMA_resampler(char _ens, int _nS, int _nC): ens(_ens), nS(_nS), nC(_nC){}

  template<typename DistributionType>
  inline void resample(DistributionType &out, const rawDataDistribution<double> &in) const{ 
    out = sampleAMAresample<DistributionType>::resample(in,ens,nS,nC);
  }
  template<typename DistributionType>
  inline DistributionType resample(const rawDataDistribution<double> &in) const{ 
    return sampleAMAresample<DistributionType>::resample(in,ens,nS,nC);
  }

  inline int getnS() const{ return nS; }
  inline int getnC() const{ return nC; }
  inline char getEns() const{ return ens; }
};

//Function to resample raw data and perform the sampleAMA correction together
template<typename resampledDistributionType>
resampledDistributionType sampleAMAresampleCorrect(const rawDataDistribution<double> &sloppy_S, const rawDataDistribution<double> &sloppy_C, const rawDataDistribution<double> &exact_C, 
						   const sampleAMA_resampler &resampler_S, const sampleAMA_resampler &resampler_C, const std::string &descr = ""){
  resampledDistributionType out_r, sloppy_S_r, sloppy_C_r, exact_C_r;
  resampler_S.resample(sloppy_S_r, sloppy_S);
  resampler_C.resample(sloppy_C_r, sloppy_C);
  resampler_C.resample(exact_C_r, exact_C);
  
  out_r = sloppy_S_r + exact_C_r - sloppy_C_r;

  if(descr != ""){
    resampledDistributionType diff = out_r - sloppy_S_r;
    resampledDistributionType reldiff = (out_r - sloppy_S_r)/sloppy_S_r;
    std::cout << descr << " corrected:" << out_r << " sloppy:" << sloppy_S_r << " diff:" << diff << " reldiff:" << reldiff << std::endl;
  }

  return out_r;
}


struct sampleAMA_resamplers{
  int nS;
  int nC;
  sampleAMA_resampler resampler_S;
  sampleAMA_resampler resampler_C;

  void setup(const int traj_start_S, const int traj_lessthan_S,
	     const int traj_start_C, const int traj_lessthan_C,
	     const int traj_inc, const int bin_size){
    nS = (traj_lessthan_S - traj_start_S)/traj_inc/bin_size;
    nC = (traj_lessthan_C - traj_start_C)/traj_inc/bin_size;

    resampler_S = sampleAMA_resampler('S',nS,nC);
    resampler_C = sampleAMA_resampler('C',nS,nC);
  }
  sampleAMA_resamplers(){}
  
  sampleAMA_resamplers(const int traj_start_S, const int traj_lessthan_S,
		      const int traj_start_C, const int traj_lessthan_C,
		       const int traj_inc, const int bin_size){
    setup(traj_start_S, traj_lessthan_S,
	  traj_start_C, traj_lessthan_C,
	  traj_inc, bin_size);
  }
};


//Perform the sampleAMA correction given the sloppy and exact data
template<typename resampledDistributionType>
inline resampledDistributionType sampleAMAresampleCorrect(const rawDataDistribution<double> &sloppy_S, const rawDataDistribution<double> &sloppy_C, const rawDataDistribution<double> &exact_C, 
							  const sampleAMA_resamplers &resamplers, const std::string &descr = ""){
  return sampleAMAresampleCorrect<resampledDistributionType>(sloppy_S,sloppy_C,exact_C,resamplers.resampler_S,resamplers.resampler_C,descr);
}




//Boost a resampled distribution from the S or C ensemble onto a sample-AMA distribution. 
//If the input distribution is a linear function of the distribution obtained from the raw data this will give exactly the same result as if the raw data had been sample-AMA resampled from the start
//otherwise it will differ slightly by 1/N^2 effects due to having to obtain the average from the resampled distribution 
template<typename DistributionType>
struct sampleAMA_boost{};

template<typename T>
struct sampleAMA_boost<jackknifeDistribution<T> >{
  static jackknifeDistribution<T> boost(const jackknifeDistribution<T> &v, const char ens, const int nS, const int nC){
    assert(ens == 'S' || ens == 'C');
    int cp_start = ens == 'S' ? 0 : nS;
    int ncp = ens == 'S' ? nS : nC;
    assert(v.size() == ncp);

    T mean = v.mean();
    jackknifeDistribution<T> out(nS+nC, mean);
    for(int i=0;i<ncp;i++) out.sample(cp_start+i) = v.sample(i);    
    return out;
  }  
};


template<typename T>
struct sampleAMA_boost<doubleJackknifeDistribution<T> >{
  static doubleJackknifeDistribution<T> boost(const doubleJackknifeDistribution<T> &v, const char ens, const int nS, const int nC){
    assert(ens == 'S' || ens == 'C');
    jackknifeDistribution<T> vjack = v.toJackknife();
    T vmean = vjack.mean();

    doubleJackknifeDistribution<T> out(nS+nC);

    if(ens == 'S'){
      assert(v.size() == nS);
      
      for(int i=0;i<nS;i++){
	int jj = 0;
	for(int j=0;j<nS;j++){ //i<nS j<nS
	  if(j!=i){
	    out.sample(i).sample(jj) = v.sample(i).sample(jj);
	    ++jj;
	  }
	}
	const T val = vjack.sample(i);
	for(int jj=nS-1;jj<nS+nC-1;jj++){ //i<nS j>=nS  (jj always j-1)
	  out.sample(i).sample(jj) = val;
	}
      }
      for(int i=nS;i<nS+nC;i++){
	for(int jj=0;jj<nS;jj++){ //i>=nS j<nS  (jj==j)
	  out.sample(i).sample(jj) = vjack.sample(jj);
	}
	for(int jj=nS;jj<nS+nC-1;jj++){ //i>=nS j>=nS (nC-1 values of j!=i)
	  out.sample(i).sample(jj) = vmean;
	}
      }
      return out;
    }else{ //ens == 'C'
      
      assert(v.size() == nC);

      for(int i=0;i<nS;i++){
	for(int jj=0;jj<nS-1;jj++){ //i<nS j<nS   (nS-1 samples where j!=i)
	  out.sample(i).sample(jj) = vmean;
	}

	int jj = nS-1;
	for(int j=nS;j<nS+nC;j++){ //i<nS j>=nS
	  out.sample(i).sample(jj++) = vjack.sample(j-nS);
	}
      }
      for(int i=nS;i<nS+nC;i++){
	const T val = vjack.sample(i-nS);
	for(int j=0;j<nS;j++){ //i>=nS j<nS    (jj=j because i>=nS and j<nS)
	  out.sample(i).sample(j) = val;
	}

	int jj=nS;
	for(int j=nS;j<nS+nC;j++){ //i>=nS j>=nS
	  if(j!=i){
	    out.sample(i).sample(jj) = v.sample(i-nS).sample(jj-nS);
	    ++jj;
	  }
	}
      }
      return out;
    }
  }
};



  





CPSFIT_END_NAMESPACE

#endif
