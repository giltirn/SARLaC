#ifndef SUPERMULTI_XML_READ_H_
#define SUPERMULTI_XML_READ_H_

//Support read from UKfit format xml files of superjackknife data
#include<distribution/superjackknife/xml_read.h>
#include<distribution/supermulti/class.h>

CPSFIT_START_NAMESPACE

//If try_reuse_layout, the reader will attempt to reuse layouts in the global layout manager if suitable
void read(XMLreader &reader, superMultiDistribution<double> &v, const std::string &tag, bool try_reuse_layout = false){
  reader.enter(tag);
#define GETIT(type,name) type name; read(reader,name,#name)

  GETIT(std::string, SampleType);
  assert(SampleType == "SuperJackBoot");
  
  GETIT(int, Nmeas);
  GETIT(int, Nensembles);
#undef GETIT

  std::vector<EnsembleData> ens;
  read(reader,ens,"Ensembles");
  
  reader.leave();

  //Setup the layout and superjack
  assert(ens.size() == Nensembles);
  
  /* static std::vector<std::unique_ptr<superMultiLayout> > layouts; //all will be deleted at the end */
  /* layouts.push_back(std::unique_ptr<superMultiLayout>(new superMultiLayout));   */
  /* superMultiLayout* layout = layouts.back().get(); */

  superMultiLayout* layout = new superMultiLayout;

  int meas=0;
  double avg;
  for(int e=0;e<ens.size();e++){
    if(e==0) avg = ens[e].avg;
    else if(ens[e].avg != avg) error_exit(std::cout << "read(XMLreader &reader, superMultiDistribution<double> &v, const std::string &tag) expect average on ensemble " << ens[e].tag << ", " << ens[e].avg << " to be equal to average of other ensembles, " << avg << std::endl);
    
    meas += ens[e].EnsembleSize;

    MultiType type;
    if(ens[e].SampleType == "Jackknife") type = MultiType::Jackknife;
    else if(ens[e].SampleType == "Bootstrap") type = MultiType::Bootstrap;
    else assert(0);

    layout->addEnsemble(ens[e].tag, type, ens[e].EnsembleSize);    
  }
  assert(meas == Nmeas);

  v = superMultiDistribution<double>(*layout, avg);
  for(int e=0;e<ens.size();e++){
    if(layout->ensType(e) == MultiType::Jackknife){
      jackknifeDistribution<double> tmp(ens[e].EnsembleSize);
      for(int s=0;s<ens[e].EnsembleSize;s++) tmp.sample(s) = ens[e].values[s];
      v.setEnsembleDistribution(e,tmp);
    }else if(layout->ensType(e) == MultiType::Bootstrap){
      bootstrapDistribution<double> tmp(bootstrapInitType(ens[e].EnsembleSize));
      tmp.best() = avg;
      for(int s=0;s<ens[e].EnsembleSize;s++) tmp.sample(s) = ens[e].values[s];
      v.setEnsembleDistribution(e,tmp);
    }else assert(0);
  }

  bool manager_append_layout = true;
  if(try_reuse_layout){
    //Look for an equivalent layout we could use to avoid always spawning new layouts for every read
    superMultiLayout const* equiv_layout = superMultiLayoutManager().findEquivalent(layout);
    if(equiv_layout != nullptr){
      v.setLayout(*equiv_layout);
      delete layout;
      manager_append_layout = false;
    }
  }
  if(manager_append_layout){
    static int tag_idx = 0;
    std::ostringstream ss; 
    ss << "xml_read_layout_" << (tag_idx++);
    superMultiLayoutManager().addLayout(layout, ss.str());
  }

}


CPSFIT_END_NAMESPACE
#endif
