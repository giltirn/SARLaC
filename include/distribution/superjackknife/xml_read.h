#ifndef SUPERJACKKNIFE_XML_READ_H_
#define SUPERJACKKNIFE_XML_READ_H_

//Support read from UKfit format xml files of superjackknife data

#include<distribution/superjackknife/class.h>
#include<serialize/xml_serialize.h>

CPSFIT_START_NAMESPACE


struct EnsembleData{
  std::string tag;
  std::string SampleType;
  int EnsembleSize;
  double avg;
  std::vector<double> values;
};
void read(XMLreader &reader, EnsembleData &v, const std::string &tag){
  //std::cout << "Reading EnsembleData with tag '" << tag << "'. Context contains:\n" << reader.printGroupEntries() << std::endl;
  reader.enter(tag);
  //std::cout << "Entered '" << tag << "'. Context now contains\n" << reader.printGroupEntries() << std::endl;
  read(reader,v.tag,"tag");
  read(reader,v.SampleType,"SampleType");
  read(reader,v.EnsembleSize,"EnsembleSize");
  read(reader,v.avg,"avg");
  read(reader,v.values,"values");
  reader.leave();
}

void read(XMLreader &reader, superJackknifeDistribution<double> &v, const std::string &tag){
  reader.enter(tag);
#define GETIT(type,name) type name; read(reader,name,#name)

  GETIT(std::string, SampleType);
  assert(SampleType == "SuperJackBoot");
  
  GETIT(int, Nmeas);
  GETIT(int, Nensembles);

  std::vector<EnsembleData> ens;
  read(reader,ens,"Ensembles");
  
  reader.leave();

  //Setup the layout and superjack
  assert(ens.size() == Nensembles);
  
  static std::vector<std::unique_ptr<superJackknifeLayout> > layouts; //all will be deleted at the end
  layouts.push_back(std::unique_ptr<superJackknifeLayout>(new superJackknifeLayout));  
  superJackknifeLayout* layout = layouts.back().get();

  int meas=0;
  double avg;
  for(int e=0;e<ens.size();e++){
    if(e==0) avg = ens[e].avg;
    else if(ens[e].avg != avg) error_exit(std::cout << "parseDoubleJackknifeXML expect average on ensemble " << ens[e].tag << ", " << ens[e].avg << " to be equal to average of other ensembles, " << avg << std::endl);
    
    meas += ens[e].EnsembleSize;
    if(ens[e].SampleType != "Jackknife") error_exit(std::cout << "read(XMLreader &reader, superJackknifeDistribution<double> &v, const std::string &tag) does not presently support sub-distributions that aren't jackknife\n");

    layout->addEnsemble(ens[e].tag,ens[e].EnsembleSize);    
  }
  assert(meas == Nmeas);

  v = superJackknifeDistribution<double>(*layout, avg);
  for(int e=0;e<ens.size();e++){
    jackknifeDistribution<double> tmp(ens[e].EnsembleSize);
    for(int s=0;s<ens[e].EnsembleSize;s++) tmp.sample(s) = ens[e].values[s];
    v.setEnsembleJackknife(e,tmp);
  }
#undef GETIT
}


CPSFIT_END_NAMESPACE
#endif
