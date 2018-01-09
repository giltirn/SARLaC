#ifndef _SAMPLEAMA_RESAMPLE_H
#define _SAMPLEAMA_RESAMPLE_H

//We consider a quantity Y that resides on nS samples, and a quantity Z that resides on a disjoint set of nC samples
//We define a super-ensemble of size nS+nC for which Y resides only on the first nS and Z only on the last nC
//The algorithm for jackknife and double-jackknife then proceeds as normal by defining reduced and doubly-reduced ensembles and averaging over the remaining data

jackknifeDistributionD superJackknifeResampleS(const rawDataDistributionD &data, const int nS, const int nC){
  assert(data.size() == nS);
  double Smean = data.mean();
  double Ssum = nS*Smean;
  
  return jackknifeDistributionD(nS+nC,
				[&](const int s){
				  if(s<nS) return 1./(nS-1)*(Ssum - data.sample(s));
				  else return Smean;
				}
				);
}
jackknifeDistributionD superJackknifeResampleC(const rawDataDistributionD &data, const int nS, const int nC){
  assert(data.size() == nC);
  double Cmean = data.mean();
  double Csum = nC*Cmean;
  
  return jackknifeDistributionD(nS+nC,
				[&](const int s){
				  if(s<nS) return Cmean;
				  else return 1./(nC-1)*(Csum - data.sample(s-nS));
				}
				);
}

doubleJackknifeDistributionD superDoubleJackknifeResampleS(const rawDataDistributionD &data, const int nS, const int nC){
  assert(data.size() == nS);
  double Smean = data.mean();
  double Ssum = nS*Smean;
  
  doubleJackknifeDistributionD out(nS+nC);
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

doubleJackknifeDistributionD superDoubleJackknifeResampleC(const rawDataDistributionD &data, const int nS, const int nC){
  assert(data.size() == nC);
  double Cmean = data.mean();
  double Csum = nC*Cmean;
  
  doubleJackknifeDistributionD out(nS+nC);
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


template<typename DistributionType>
struct sampleAMAresample{};

template<>
struct sampleAMAresample<jackknifeDistributionD>{
  static inline jackknifeDistributionD resample(const rawDataDistributionD &in, const char ens, const int nS, const int nC){
    return ens == 'S' ? superJackknifeResampleS(in,nS,nC) : superJackknifeResampleC(in,nS,nC);
  }
};
template<>
struct sampleAMAresample<doubleJackknifeDistributionD>{
  static inline doubleJackknifeDistributionD resample(const rawDataDistributionD &in, const char ens, const int nS, const int nC){
    return ens == 'S' ? superDoubleJackknifeResampleS(in,nS,nC) : superDoubleJackknifeResampleC(in,nS,nC);
  }
};


#endif
