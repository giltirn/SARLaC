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

doubleJackknifeDistributionD superDoubleJackknifeResampleS_orig(const rawDataDistributionD &data, const int nS, const int nC){
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

doubleJackknifeDistributionD superDoubleJackknifeResampleS(const rawDataDistributionD &data, const int nS, const int nC){
  assert(data.size() == nS);
  const double Smean = data.mean();
  const double Ssum = nS*Smean;
  const double invnSm2 = 1./(nS-2);
  const double invnSm1 = 1./(nS-1);

  doubleJackknifeDistributionD out(nS+nC);
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





doubleJackknifeDistributionD superDoubleJackknifeResampleC_orig(const rawDataDistributionD &data, const int nS, const int nC){
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


doubleJackknifeDistributionD superDoubleJackknifeResampleC(const rawDataDistributionD &data, const int nS, const int nC){
  assert(data.size() == nC);
  const double Cmean = data.mean();
  const double Csum = nC*Cmean;
  const double invnCm2 = 1./(nC-2);
  const double invnCm1 = 1./(nC-1);

  doubleJackknifeDistributionD out(nS+nC);
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
