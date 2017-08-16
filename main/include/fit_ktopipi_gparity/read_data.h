#ifndef _FIT_KTOPIPI_READ_DATA_H
#define _FIT_KTOPIPI_READ_DATA_H

inline std::string typeFile(const int traj, const int type, const int tsep_k_pi, const int tsep_pipi, const std::string &data_dir, const threeMomentum &mom = {0,0,0}){
  std::ostringstream os;
  os << data_dir << "/traj_" << traj << "_type" << type;
  if(type != 4) os << "_deltat_" << tsep_k_pi << "_sep_" << tsep_pipi;
  if(type == 1) os << "_mom" << momStr(mom);
  return os.str();
}

type1234Data readType(const int type, const int traj_start, const int traj_inc, const int traj_lessthan, const int tsep_k_pi, const int tsep_pipi, const int Lt, const std::string &data_dir){
  int nsample = (traj_lessthan - traj_start)/traj_inc;
  if(type == 1){
    type1234Data type1_mom111(1,Lt,nsample);    
    type1234Data type1_mom_1_1_1(1,Lt,nsample);
    type1234Data type1_mom_111(1,Lt,nsample);
    type1234Data type1_mom1_1_1(1,Lt,nsample);

#pragma omp parallel for
    for(int i=0;i<nsample;i++){
      const int traj = traj_start + i*traj_inc;
      type1_mom111.parse(typeFile(traj,type,tsep_k_pi,tsep_pipi,data_dir,{1,1,1}),i);
      type1_mom_1_1_1.parse(typeFile(traj,type,tsep_k_pi,tsep_pipi,data_dir,{-1,-1,-1}),i);
      type1_mom_111.parse(typeFile(traj,type,tsep_k_pi,tsep_pipi,data_dir,{-1,1,1}),i);
      type1_mom1_1_1.parse(typeFile(traj,type,tsep_k_pi,tsep_pipi,data_dir,{1,-1,-1}),i);
    }
    return type1_mom111 * (1.0/8.0) + type1_mom_1_1_1 * (1.0/8.0) + type1_mom_111 * (3.0/8.0) + type1_mom1_1_1 * (3.0/8.0); //Project onto A2
  }else if(type == 2 || type == 3 || type == 4){
    type1234Data typedata(type,Lt,nsample);
#pragma omp parallel for
    for(int i=0;i<nsample;i++){
      const int traj = traj_start + i*traj_inc;
      typedata.parse(typeFile(traj,type,tsep_k_pi,tsep_pipi,data_dir),i);
    }
#ifdef DAIQIAN_COMPATIBILITY_MODE
    if(type == 2 || type == 3) typedata = typedata * 0.5; //correct for missing coefficient
#endif
    return typedata;
  }else{
    error_exit(std::cout << "readType invalid type " << type << std::endl);
  }
}


NumericTensor<rawDataDistributionD,1> readA2projectedBubble(const int traj_start, const int traj_inc, const int traj_lessthan, const int tsep_pipi, const int Lt, const std::string &data_dir){
  bubbleDataAllMomenta raw_bubble_data;
  readBubble(raw_bubble_data, data_dir, tsep_pipi, Lt, traj_start, traj_inc, traj_lessthan);

  NumericTensor<rawDataDistributionD,1> out({Lt});
  for(int t=0;t<Lt;t++) out({t}) = ( raw_bubble_data({1,1,1})(t)  + raw_bubble_data({-1,-1,-1})(t)
				     + raw_bubble_data({-1,1,1})(t) + raw_bubble_data({1,-1,-1})(t)
				     + raw_bubble_data({1,-1,1})(t) + raw_bubble_data({-1,1,-1})(t)
				     + raw_bubble_data({1,1,-1})(t) + raw_bubble_data({-1,-1,1})(t) )/8.;
  return out;
}


#endif
