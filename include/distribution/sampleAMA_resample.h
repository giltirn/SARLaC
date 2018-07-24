#ifndef _SAMPLEAMA_RESAMPLE_H
#define _SAMPLEAMA_RESAMPLE_H

#include<config.h>
#include<utils/macros.h>

#include<distribution/raw_data_distribution.h>
#include<distribution/jackknife.h>
#include<distribution/double_jackknife.h>

CPSFIT_START_NAMESPACE

//We consider a quantity Y that resides on nS samples, and a quantity Z that resides on a disjoint set of nC samples
//We define a super-ensemble of size nS+nC for which Y resides only on the first nS and Z only on the last nC
//The algorithm for jackknife and double-jackknife then proceeds as normal by defining reduced and doubly-reduced ensembles and averaging over the remaining data

jackknifeDistribution<double> superJackknifeResampleS(const rawDataDistribution<double> &data, const int nS, const int nC){
  assert(data.size() == nS);
  double Smean = data.mean();
  double Ssum = nS*Smean;
  
  return jackknifeDistribution<double>(nS+nC,
				[&](const int s){
				  if(s<nS) return 1./(nS-1)*(Ssum - data.sample(s));
				  else return Smean;
				}
				);
}
jackknifeDistribution<double> superJackknifeResampleC(const rawDataDistribution<double> &data, const int nS, const int nC){
  assert(data.size() == nC);
  double Cmean = data.mean();
  double Csum = nC*Cmean;
  
  return jackknifeDistribution<double>(nS+nC,
				[&](const int s){
				  if(s<nS) return Cmean;
				  else return 1./(nC-1)*(Csum - data.sample(s-nS));
				}
				);
}

doubleJackknifeDistribution<double> superDoubleJackknifeResampleS_orig(const rawDataDistribution<double> &data, const int nS, const int nC){
  assert(data.size() == nS);
  double Smean = data.mean();
  double Ssum = nS*Smean;
  
  doubleJackknifeDistribution<double> out(nS+nC);
  for(int i=0;i<nS+nC;i++){
    int jj=0;
    for(int j=0;j<nS+nC;j++){
      if(j==i) continue;

      if(i<nS && j<nS) out.sample(i).sample(jj) = 1./(nS-2)*(Ssum - data.sample(i) - data.sample(j)); 
      else if(i<nS && j>= nS) out.sample(i).sample(jj) = 1./(nS-1)*(Ssum - data.sample(i));
      else if(i>=nS && j< nS) out.sample(i).sample(jj) = 1./(nS-1)*(Ssum - data.sample(j));
      else out.sample(i).sample(jj) = Smean;
      ++jj;
    }
  }
  return out;
}

//A faster implementation of the above
doubleJackknifeDistribution<double> superDoubleJackknifeResampleS(const rawDataDistribution<double> &data, const int nS, const int nC){
  assert(data.size() == nS);
  const double Smean = data.mean();
  const double Ssum = nS*Smean;
  const double invnSm2 = 1./(nS-2);
  const double invnSm1 = 1./(nS-1);

  doubleJackknifeDistribution<double> out(nS+nC);
  for(int i=0;i<nS;i++){
    int jj = 0;
    for(int j=0;j<nS;j++){ //i<nS j<nS
      if(j!=i) out.sample(i).sample(jj++) = invnSm2*(Ssum - data.sample(i) - data.sample(j)); 
    }
    const double val = invnSm1*(Ssum - data.sample(i));
    for(int jj=nS-1;jj<nS+nC-1;jj++){ //i<nS j>=nS  (jj always j-1)
      out.sample(i).sample(jj) = val;
    }
  }
  for(int i=nS;i<nS+nC;i++){
    for(int jj=0;jj<nS;jj++){ //i>=nS j<nS  (jj==j)
      out.sample(i).sample(jj) = invnSm1*(Ssum - data.sample(jj));
    }
    for(int jj=nS;jj<nS+nC-1;jj++){ //i>=nS j>=nS (nC-1 values of j!=i)
      out.sample(i).sample(jj) = Smean;
    }
  }
  return out;
}





doubleJackknifeDistribution<double> superDoubleJackknifeResampleC_orig(const rawDataDistribution<double> &data, const int nS, const int nC){
  assert(data.size() == nC);
  double Cmean = data.mean();
  double Csum = nC*Cmean;
  
  doubleJackknifeDistribution<double> out(nS+nC);
  for(int i=0;i<nS+nC;i++){
    int jj=0;
    for(int j=0;j<nS+nC;j++){
      if(j==i) continue;

      if(i<nS && j<nS) out.sample(i).sample(jj) = Cmean;
      else if(i<nS && j>= nS) out.sample(i).sample(jj) = 1./(nC-1)*(Csum - data.sample(j-nS));
      else if(i>=nS && j< nS) out.sample(i).sample(jj) = 1./(nC-1)*(Csum - data.sample(i-nS));
      else out.sample(i).sample(jj) = 1./(nC-2)*(Csum - data.sample(i-nS) - data.sample(j-nS)); 
      ++jj;
    }
  }
  return out;
}

//A faster implementation of the above
doubleJackknifeDistribution<double> superDoubleJackknifeResampleC(const rawDataDistribution<double> &data, const int nS, const int nC){
  assert(data.size() == nC);
  const double Cmean = data.mean();
  const double Csum = nC*Cmean;
  const double invnCm2 = 1./(nC-2);
  const double invnCm1 = 1./(nC-1);

  doubleJackknifeDistribution<double> out(nS+nC);
  for(int i=0;i<nS;i++){
    for(int jj=0;jj<nS-1;jj++){ //i<nS j<nS   (nS-1 samples where j!=i)
      out.sample(i).sample(jj) = Cmean;
    }

    int jj = nS-1;
    for(int j=nS;j<nS+nC;j++){ //i<nS j>=nS
      out.sample(i).sample(jj++) = invnCm1*(Csum - data.sample(j-nS));
    }
  }
  for(int i=nS;i<nS+nC;i++){
    const double val = invnCm1*(Csum - data.sample(i-nS));
    for(int j=0;j<nS;j++){ //i>=nS j<nS    (jj=j because i>=nS and j<nS)
      out.sample(i).sample(j) = val;
    }

    int jj=nS;
    for(int j=nS;j<nS+nC;j++){ //i>=nS j>=nS
      if(j!=i) out.sample(i).sample(jj++) = invnCm2*(Csum - data.sample(i-nS) - data.sample(j-nS));
    }
  }
  return out;
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


template<typename resampledDistributionType>
inline resampledDistributionType sampleAMAresampleCorrect(const rawDataDistribution<double> &sloppy_S, const rawDataDistribution<double> &sloppy_C, const rawDataDistribution<double> &exact_C, 
							  const sampleAMA_resamplers &resamplers, const std::string &descr = ""){
  return sampleAMAresampleCorrect<resampledDistributionType>(sloppy_S,sloppy_C,exact_C,resamplers.resampler_S,resamplers.resampler_C,descr);
}

CPSFIT_END_NAMESPACE

#endif
