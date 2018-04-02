#ifndef _CPSFIT_SUPERJACKKNIFE_QDP_BINARY_READ_H_
#define _CPSFIT_SUPERJACKKNIFE_QDP_BINARY_READ_H_

#include<config.h>
#include<utils/macros.h>
#include<serialize/qdp_binary_read.h>
#include<distribution/superjackknife/class.h>

CPSFIT_START_NAMESPACE

template<typename T>
void read(QDPbinaryReader &rd, superJackknifeDistribution<T> &into){
  //Get basic info and the samples
  int sample_type;
  read(rd, sample_type);
  if(sample_type != 4) error_exit(std::cout << "read(QDPbinaryReader &, superJackknifeDistribution<T> &) failed: Expected sample_type=4, got " << sample_type << std::endl);
  int nmeas;
  read(rd, nmeas);

  std::vector<T> samples(nmeas);
  for(int i=0;i<nmeas;i++) read(rd,samples[i]);
  T avg;
  read(rd,avg);

  //Read the layout information
  static std::vector<std::unique_ptr<superJackknifeLayout> > layouts; //all will be deleted at the end
  layouts.push_back(std::unique_ptr<superJackknifeLayout>(new superJackknifeLayout));  
  superJackknifeLayout* layout = layouts.back().get();

  int nens;
  read(rd, nens);

  for(int e=0;e<nens;e++){ //check all sub-ensembles are jackknife
    int sub_type;
    read(rd, sub_type);
    if(sub_type != 2) error_exit(std::cout << "read(QDPbinaryReader &, superJackknifeDistribution<T> &) failed: Expected sub-ensemble sample_type=2, got " << sub_type << std::endl);
  }

  int infosz;
  read(rd, infosz);
  assert(infosz == nens);

  for(int e=0;e<nens;e++){ //setup the layout
    int nraw, nbinned, nresampled, binningfactor,sub_type;
    read(rd, nraw);
    read(rd, nbinned);
    read(rd, nresampled);
    read(rd, binningfactor);
    read(rd, sub_type);
    if(sub_type != 2) error_exit(std::cout << "read(QDPbinaryReader &, superJackknifeDistribution<T> &) failed: Expected sub-ensemble sample_type=2, got " << sub_type << std::endl);
    std::string tag;
    read(rd,tag);

    layout->addEnsemble(tag, nresampled);
  }

  //Setup the superjack
  into.setLayout(*layout);
  assert(nmeas == into.size());
  for(int s=0;s<nmeas;s++) into.sample(s) = samples[s];
  into.best() = avg;
}

CPSFIT_END_NAMESPACE

#endif
