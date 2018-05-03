#ifndef _GLOBALFIT_READ_REWEIGHTED_H_
#define _GLOBALFIT_READ_REWEIGHTED_H_

struct ReweightInfo{
  double mh_start;
  double mh_end;
  double mh_incr;
  double mh_data;
  ReweightInfo(double mh_start, double mh_end, double mh_incr, double mh_data): mh_start(mh_start), mh_end(mh_end), mh_incr(mh_incr), mh_data(mh_data){}
};

void addReweightedData(DataSeriesType &data, 
		       const std::vector< superJackknifeDistribution<double> > &rw_data, const std::vector< superJackknifeDistribution<double> > &rw_mres,
		       const DataParams &base_params, const ReweightInfo &rw_info, const std::vector<int> &rw_idx_keep){
  assert(rw_mres.size() >= rw_idx_keep.size());
  for(int i=0;i<rw_idx_keep.size();i++){
    DataParams iparams(base_params);
    iparams.mh = rw_info.mh_start + rw_idx_keep[i]*rw_info.mh_incr;
    
    std::cout << "mh="<<iparams.mh << " value=" << rw_data[ rw_idx_keep[i] ] << " (idx=" << rw_idx_keep[i] << ")" << std::endl;
    data.push_back(iparams + rw_mres[ rw_idx_keep[i] ], rw_data[ rw_idx_keep[i] ]);
  }
}

struct readAddDoNothing{
  inline void operator()(jackknifeDistribution<double> &val) const{}
};


template<typename Op = readAddDoNothing>
void readAddReweighedDataQDP(DataSeriesType &data, const std::string &file, const int parameter_idx, const DataParams &base_params, const std::string &ens, const ReweightInfo &rw_info, const std::vector<int> &rw_idx_keep, const std::vector< superJackknifeDistribution<double> > &rw_mres, const Op &op = readAddDoNothing()){
  assert(rw_mres.size() == rw_idx_keep.size() && rw_mres.size() > 0);
  const superJackknifeLayout &layout = rw_mres[0].getLayout();

  std::cout << "Reading data from file " << file << std::endl;
  QDPbinaryReader rd(file);
  std::vector<std::vector<jackknifeDistribution<double> > > data_j;
  read(rd,data_j);

  for(int i=0;i<rw_idx_keep.size();i++){
    int idx = rw_idx_keep[i];

    op(data_j[parameter_idx][idx]);
    
    DataParams iparams(base_params);
    iparams.mh = rw_info.mh_start + idx*rw_info.mh_incr;
    
    superJackknifeDistribution<double> val(layout,0.);
    val.setEnsembleJackknife(ens, data_j[parameter_idx][idx]);
    val.best() = data_j[parameter_idx][idx].best();

    std::cout << "mh="<<iparams.mh << " value=" << val << " (idx=" << idx << ")" << std::endl;

    //data.push_back(iparams + rw_mres[ rw_idx_keep[i] ], val);
    data.push_back(iparams + rw_mres[i], val);
  }
}

// std::vector<superJackknifeDistribution<double> > readReweightedMresQDP(const std::string &file, const int outer_param = 1){
//   std::vector<std::vector<superJackknifeDistribution<double> > > tmp;
//   QDPbinaryReader rd(file);
//   read(rd,tmp);
//   return std::move(tmp[outer_param]);
// }
std::vector<superJackknifeDistribution<double> > readReweightedMresQDP(const std::string &file, const std::vector<int> &rw_idx_keep, const int outer_param = 1){
  std::cout << "Reading mres from " << file << std::endl;
  std::vector<std::vector<superJackknifeDistribution<double> > > tmp;
  QDPbinaryReader rd(file);
  read(rd,tmp);

  std::vector<superJackknifeDistribution<double> > out(rw_idx_keep.size());
  for(int i=0;i<rw_idx_keep.size();i++)
    out[i] = std::move(tmp[outer_param][ rw_idx_keep[i] ]);

  return out;
}
#endif
