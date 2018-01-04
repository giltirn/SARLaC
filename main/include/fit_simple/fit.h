#ifndef _FIT_SIMPLE_FIT_H_
#define _FIT_SIMPLE_FIT_H_

inline int getMassParamIdx(const FitFuncType type){
  switch(type){
  case FCosh:
  case FSinh:
  case FExp:
    return 1;
  default:
    error_exit(std::cout << "getMassParamIdx unknown fit func " << type << std::endl);
  }
}
template<typename FitFunc>
FitFunc* getFitFunc(const Args &args){ assert(0); }

template<>
FitCosh* getFitFunc<FitCosh>(const Args &args){ return new FitCosh(args.Lt); }

template<>
FitSinh* getFitFunc<FitSinh>(const Args &args){ return new FitSinh(args.Lt); }

template<>
FitExp* getFitFunc<FitExp>(const Args &args){ return new FitExp; }

template<typename FitFunc, template<typename> class CostFunctionPolicy>
void fitSpecFFcorr(const jackknifeCorrelationFunctionD &data_j, const doubleJackknifeCorrelationFunctionD &data_dj, const Args &args);

template<typename FitFunc>
inline void fitSpecFF(const jackknifeCorrelationFunctionD &data_j, const doubleJackknifeCorrelationFunctionD &data_dj, const Args &args, const CMDline &cmdline){
  return args.correlated ? fitSpecFFcorr<FitFunc,correlatedFitPolicy>(data_j,data_dj,args,cmdline) : fitSpecFFcorr<FitFunc,uncorrelatedFitPolicy>(data_j, data_dj,args,cmdline);
};
inline void fitResampled(const jackknifeCorrelationFunctionD &data_j, const doubleJackknifeCorrelationFunctionD &data_dj, const Args &args, const CMDline &cmdline){
  switch(args.fitfunc){
  case FCosh:
    return fitSpecFF<FitCosh>(data_j, data_dj,args,cmdline);
  case FSinh:
    return fitSpecFF<FitSinh>(data_j, data_dj,args,cmdline);
  case FExp:
    return fitSpecFF<FitExp>(data_j, data_dj,args,cmdline);    
  default:
    error_exit(std::cout << "fit: Invalid fitfunc " << args.fitfunc << std::endl);
  };
}
inline void fit(const rawDataCorrelationFunctionD &data, const Args &args, const CMDline &cmdline){
  doubleJackknifeCorrelationFunctionD data_dj(args.Lt, 
					      [&](const int t){
						return typename doubleJackknifeCorrelationFunctionD::ElementType(t,  doubleJackknifeDistributionD(data.value(t)));
					      }
					      );
  jackknifeCorrelationFunctionD data_j(args.Lt, 
				       [&](const int t){
					 return typename jackknifeCorrelationFunctionD::ElementType(t,  jackknifeDistributionD(data.value(t))); 
				       }
				       );

  if(cmdline.save_combined_data){
#ifdef HAVE_HDF5
    std::cout << "Writing double-jackknife data to " << cmdline.save_combined_data_file << std::endl;
    HDF5writer writer(cmdline.save_combined_data_file);
    write(writer, data_dj, "data");
#else
    error_exit("fitSpecFFcorr: Saving amplitude data requires HDF5\n");
#endif
  }
  fitResampled(data_j, data_dj, args, cmdline);
}
  
#endif
