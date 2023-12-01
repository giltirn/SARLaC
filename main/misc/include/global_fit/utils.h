#ifndef _GLOBALFIT_UTILS_H
#define _GLOBALFIT_UTILS_H_

//Convert double to string without trailing zeros
inline std::string dstr(const double v){
  std::ostringstream ss; ss << v;
  std::string vs = ss.str();
  int nzro = 0;
  int i=vs.size()-1;
  while(vs[i] == '0'){
    ++nzro;
    --i;
  }
  return vs.erase(vs.size()-nzro,std::string::npos);
}

//Convert the layout of a superjackknife with a provided mapping between the ensemble tags
superJackknifeDistribution<double> remapLayout(const superJackknifeDistribution<double> &sj, const superJackknifeLayout &target, const std::map<std::string, std::string> &map_orig_new){
  superJackknifeDistribution<double> out(target, sj.best());

  const superJackknifeLayout &layout_in = sj.getLayout();
  for(int orig_idx=0;orig_idx<layout_in.nEnsembles();orig_idx++){
    const std::string &orig_tag = layout_in.ensTag(orig_idx);
    int orig_ens_size = layout_in.nSamplesEns(orig_idx);

    //Find the new tag
    std::map<std::string, std::string>::const_iterator mpit = map_orig_new.find(orig_tag);
    if(mpit == map_orig_new.end()) error_exit(std::cout << "remapLayout(..) Could not find mapping for ensemble " << orig_tag << std::endl);

    //Check remapped ensemble exists in target and has the same size
    const std::string &rmp_tag = mpit->second;
    int targ_idx = target.ensIdx(rmp_tag);
    if(targ_idx == -1) error_exit(std::cout << "remapLayout(..) Mapped ensemble tag " << rmp_tag << " does not exist in target layout " << target << std::endl);
    
    int targ_ens_size = target.nSamplesEns(targ_idx);
    if(targ_ens_size != orig_ens_size) error_exit(std::cout << "remapLayout(..) Mapped ensemble " << rmp_tag << " size " << targ_ens_size << " does not match original " << orig_tag << " with size " << orig_ens_size << std::endl);

    //Poke sub-ensemble jackknife into output
    out.setEnsembleJackknife(targ_idx, sj.getEnsembleJackknife(orig_idx));
  }
  return out;
}

#define _ENUMERATED_STRUCT_BOOST_MEMBER(R,DUMMY,IDX,ELEM) T ELEM = b(p. ELEM, IDX);
#define ENUMERATED_STRUCT_BOOST_MEMBERS(DEF) BOOST_PP_SEQ_FOR_EACH_I(_ENUMERATED_STRUCT_BOOST_MEMBER, , _ENUMERATED_STRUCT_DEF_SEQ(DEF))



#endif
