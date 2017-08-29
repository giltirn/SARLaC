#ifndef _FIT_KTOPIPI_READ_DATA_H
#define _FIT_KTOPIPI_READ_DATA_H

#include<algorithm>
#include<expression_parse.h>

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


#ifdef HAVE_HDF5
void writeBubble(HDF5writer &writer, const NumericTensor<rawDataDistributionD,1> &value, const std::string &tag){
  boost::timer::auto_cpu_timer total_time("writeBubble(HDF5writer &, const NumericTensor<rawDataDistributionD,1> &, const std::string &)  %w s\n");
  writer.enter(tag);
  int size = value.size(0);
  int nsample = value({0}).size();
  
  write(writer, size, "size");
  write(writer, nsample, "nsample");

  std::vector<double> data(size * nsample);
  for(int i=0;i<size;i++)
    for(int s=0;s<nsample;s++)
      data[s + nsample*i] = value(&i).sample(s);
  write(writer, data, "data");
  writer.leave();
}
void readBubble(HDF5reader &reader, NumericTensor<rawDataDistributionD,1> &value, const std::string &tag){
  boost::timer::auto_cpu_timer total_time("readBubble(HDF5writer &, NumericTensor<rawDataDistributionD,1> &, const std::string &)  %w s\n");
  reader.enter(tag);
  int size;
  int nsample;
  
  read(reader, size, "size");
  read(reader, nsample, "nsample");

  std::vector<double> data(size * nsample);
  read(reader, data, "data");

  value.resize(&size, rawDataDistributionD(nsample));
  
  for(int i=0;i<size;i++)
    for(int s=0;s<nsample;s++)
      value(&i).sample(s) = data[s + nsample*i];

  reader.leave();
}
#endif

NumericTensor<rawDataDistributionD,1> getA2projectedBubble(const Args &args, const CMDline &cmdline){
  NumericTensor<rawDataDistributionD,1> bubble;
  if(cmdline.load_data_checkpoint){
#ifdef HAVE_HDF5
    std::ostringstream file; file << cmdline.load_data_checkpoint_stub << "_bubble.hdf5";
    std::cout << "Loading checkpoint data for bubble from " << file.str() << std::endl;
    HDF5reader rd(file.str());
    readBubble(rd,bubble,"bubble");
#else
    error_exit(std::cout << "Checkpointing of data requires HDF5\n");
#endif
  }else{
    bubble = readA2projectedBubble(args.traj_start,args.traj_inc,args.traj_lessthan,args.tsep_pipi,args.Lt,args.data_dir);
  }
  if(cmdline.save_data_checkpoint){
#ifdef HAVE_HDF5
    std::ostringstream file; file << cmdline.save_data_checkpoint_stub << "_bubble.hdf5";
    std::cout << "Saving checkpoint data for bubble to " << file.str() << std::endl;
    HDF5writer wr(file.str());
    writeBubble(wr,bubble,"bubble");
#else
    error_exit(std::cout << "Checkpointing of data requires HDF5\n");
#endif
  }
  return bubble;
}

void readUKfitVectorEntry(jackknifeDistribution<double> &into, const std::string &filename, const int idx){
  XMLreader rd(filename);
  std::cout << "readUKfitVectorEntry reading file " << filename << std::endl;
  UKvalenceDistributionContainer<jackknifeDistribution<double> > con;
  read(rd,con,"data_in_file");
  std::cout << "Read " << con.Nentries << " distributions from file " << filename << std::endl;  
  into = con.list[idx];
}


template<typename FitFuncPolicies>
void readFrozenParams(fitter<FitFuncPolicies> &fitter, const int Q, const CMDline &cmdline, const int nsample){
  if(!cmdline.load_freeze_data) return;

  typedef typename FitFuncPolicies::jackknifeFitParameters jackknifeFitParameters;
  jackknifeFitParameters values(nsample);
  
  std::vector<int> freeze;
  FreezeParams fparams;
  {
    std::ifstream fi(cmdline.freeze_data.c_str());
    fi >> fparams;
    assert(!fi.fail());
  }
  for(int i=0;i<fparams.sources.size();i++){
    if(fparams.sources[i].Qlist.size() > 0 && std::find(fparams.sources[i].Qlist.begin(), fparams.sources[i].Qlist.end(), Q) == fparams.sources[i].Qlist.end()) continue;
    
    std::cout << "readFrozenParams loading freeze data for parameter " << fparams.sources[i].param_idx << std::endl;
    freeze.push_back(fparams.sources[i].param_idx);

    jackknifeDistribution<double> fval;
    if(fparams.sources[i].reader == UKfitXMLvectorReader) readUKfitVectorEntry(fval, fparams.sources[i].filename, fparams.sources[i].input_idx);
    else error_exit(std::cout << "readFrozenParams unknown reader " << fparams.sources[i].reader << std::endl);

    if(fval.size() != nsample) error_exit(std::cout << "readFrozenParams read jackknife of size " << fval.size() << ", expected " << nsample << std::endl);

    if(fparams.sources[i].operation != ""){
      expressionAST AST = mathExpressionParse(fparams.sources[i].operation);
      if(!AST.containsSymbol("x")) error_exit(std::cout << "readFrozenParams expect math expression to be a function of 'x', got \"" << fparams.sources[i].operation << "\"\n");
      for(int s=0;s<nsample;s++){
	AST.assignSymbol("x",fval.sample(s));
	fval.sample(s) = AST.value();
      }
    }

    std::cout << "readFrozenParams read " << fval << std::endl;
    
    for(int s=0;s<nsample;s++)
      values.sample(s)(fparams.sources[i].param_idx) = fval.sample(s);    
  }

  fitter.freeze(freeze, values);
}


#endif
