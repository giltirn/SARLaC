#ifndef _FIT_PIPI_GND_EXC_SIGMA_GPARITY_FITFUNC_H
#define _FIT_PIPI_GND_EXC_SIGMA_GPARITY_FITFUNC_H

//Define the mappings between the sub-fit parameters and the full set of parameters plus index the parameters
typedef std::unordered_map<std::string, std::string> SubFitFuncParameterMap;
typedef std::unordered_map<std::string,size_t> ParamTagIdxMap;
typedef taggedValueContainer<double,std::string> Params;

inline std::string getStub(const Operator op){
  switch(op){
  case Operator::PiPiGnd:
    return "pipi_gnd";
  case Operator::PiPiExc:
    return "pipi_exc";
  case Operator::Sigma:
    return "sigma";
  }
}

//nstate applies only to "MultiState" variants
void setupParameterMaps(std::map< std::pair<Operator,Operator>, SubFitFuncParameterMap > &subfit_pmaps,
			ParamTagIdxMap &param_map,
			Params &guess,
			const std::vector<Operator> &incl_ops,
			const FitFuncType fitfunc, const double Ascale, const int nstate){
#define DEF(NM, GUESS) { param_map[NM] = p++;  guesses.push_back(GUESS); }  

  std::vector<double> guesses;
  int nops = incl_ops.size();

  //Define the full set of parameters and default guesses  
  if(fitfunc == FitFuncType::FSimGenOneState){ 
    int p=0;
    for(int i=0;i<nops;i++) DEF("A" + getStub(incl_ops[i]) + "_0", 1e6 / sqrt(Ascale)); 
    DEF("E0",0.35);
    for(int i=0;i<nops;i++)
      for(int j=i;j<nops;j++)
	DEF("C" + getStub(incl_ops[i]) + "_" + getStub(incl_ops[j]), 0.);
      

    //Map internal sub-fit parameters to outer params
    for(int i=0;i<nops;i++){
      auto istub = getStub(incl_ops[i]);
      for(int j=i;j<nops;j++){
	auto jstub = getStub(incl_ops[j]);
	subfit_pmaps[{incl_ops[i],incl_ops[j]}] =
	  {  {"Asrc","A" + istub + "_0"}, {"Asnk", "A" + jstub + "_0"}, {"E", "E0"},	     
	     {"Csys", "C" + istub + "_" + jstub} };
      }
    }


  }else if(fitfunc == FitFuncType::FSimGenTwoState){ 
    int p=0;
    for(int i=0;i<nops;i++) DEF("A" + getStub(incl_ops[i]) + "_0", 1e6 / sqrt(Ascale)); 
    DEF("E0",0.35);
    for(int i=0;i<nops;i++) DEF("A" + getStub(incl_ops[i]) + "_1", 1e6 / sqrt(Ascale));
    DEF("E1", 0.7);
    for(int i=0;i<nops;i++)
      for(int j=i;j<nops;j++)
	DEF("C" + getStub(incl_ops[i]) + "_" + getStub(incl_ops[j]), 0.);
      

    //Map internal sub-fit parameters to outer params
    for(int i=0;i<nops;i++){
      auto istub = getStub(incl_ops[i]);
      for(int j=i;j<nops;j++){
	auto jstub = getStub(incl_ops[j]);
	subfit_pmaps[{incl_ops[i],incl_ops[j]}] =
	  {  {"Asrc0","A" + istub + "_0"}, {"Asnk0", "A" + jstub + "_0"}, {"E0", "E0"},
	     {"Asrc1","A" + istub + "_1"}, {"Asnk1", "A" + jstub + "_1"}, {"E1", "E1"},
	     {"Csys", "C" + istub + "_" + jstub} };
      }
    }

  }else if(fitfunc == FitFuncType::FSimGenThreeState){ 
    int p=0;
    for(int i=0;i<nops;i++) DEF("A" + getStub(incl_ops[i]) + "_0", 1e6 / sqrt(Ascale)); 
    DEF("E0",0.35);
    for(int i=0;i<nops;i++) DEF("A" + getStub(incl_ops[i]) + "_1", 1e6 / sqrt(Ascale));
    DEF("E1", 0.7);
    for(int i=0;i<nops;i++) DEF("A" + getStub(incl_ops[i]) + "_2", 1e6 / sqrt(Ascale));
    DEF("E2", 0.9);
    for(int i=0;i<nops;i++)
      for(int j=i;j<nops;j++)
	DEF("C" + getStub(incl_ops[i]) + "_" + getStub(incl_ops[j]), 0.);
      
    //Map internal sub-fit parameters to outer params
    for(int i=0;i<nops;i++){
      auto istub = getStub(incl_ops[i]);
      for(int j=i;j<nops;j++){
	auto jstub = getStub(incl_ops[j]);
	subfit_pmaps[{incl_ops[i],incl_ops[j]}] =
	  {  {"Asrc0","A" + istub + "_0"}, {"Asnk0", "A" + jstub + "_0"}, {"E0", "E0"},
	     {"Asrc1","A" + istub + "_1"}, {"Asnk1", "A" + jstub + "_1"}, {"E1", "E1"},
	     {"Asrc2","A" + istub + "_2"}, {"Asnk2", "A" + jstub + "_2"}, {"E2", "E2"},
	     {"Csys", "C" + istub + "_" + jstub} };
      }
    }  
  }else if(fitfunc == FitFuncType::FSimGenThreeStateLogEdiff){ 
    int p=0;
    for(int i=0;i<nops;i++) DEF("A" + getStub(incl_ops[i]) + "_0", 1e6 / sqrt(Ascale)); 
    DEF("logE0",-1.05);
    for(int i=0;i<nops;i++) DEF("A" + getStub(incl_ops[i]) + "_1", 1e6 / sqrt(Ascale));
    DEF("logE1mE0", -1.05);
    for(int i=0;i<nops;i++) DEF("A" + getStub(incl_ops[i]) + "_2", 1e6 / sqrt(Ascale));
    DEF("logE2mE1", -1.6);
    for(int i=0;i<nops;i++)
      for(int j=i;j<nops;j++)
	DEF("C" + getStub(incl_ops[i]) + "_" + getStub(incl_ops[j]), 0.);
      
    //Map internal sub-fit parameters to outer params
    for(int i=0;i<nops;i++){
      auto istub = getStub(incl_ops[i]);
      for(int j=i;j<nops;j++){
	auto jstub = getStub(incl_ops[j]);
	subfit_pmaps[{incl_ops[i],incl_ops[j]}] =
	  {  {"Asrc0","A" + istub + "_0"}, {"Asnk0", "A" + jstub + "_0"}, {"logE0", "logE0"},
	     {"Asrc1","A" + istub + "_1"}, {"Asnk1", "A" + jstub + "_1"}, {"logE1mE0", "logE1mE0"},
	     {"Asrc2","A" + istub + "_2"}, {"Asnk2", "A" + jstub + "_2"}, {"logE2mE1", "logE2mE1"},
	     {"Csys", "C" + istub + "_" + jstub} };
      }
    }
  }else if(fitfunc == FitFuncType::FSimGenMultiState || fitfunc == FitFuncType::FSimGenMultiStateSub){ 
    //Define the full set of outer parameters
    int p=0;
    for(int state=0;state<nstate;state++){
      std::string statestr = anyToStr(state);
      for(int i=0;i<nops;i++) DEF("A" + getStub(incl_ops[i]) + "_" + statestr, 1e6 / sqrt(Ascale)); 
      DEF("E" + statestr,0.3 * (state+1));
    }
    for(int i=0;i<nops;i++)
      for(int j=i;j<nops;j++)
	DEF("C" + getStub(incl_ops[i]) + "_" + getStub(incl_ops[j]), 0.);
      
    //Map internal sub-fit parameters to outer params
    for(int i=0;i<nops;i++){
      auto istub = getStub(incl_ops[i]);
      for(int j=i;j<nops;j++){
	auto jstub = getStub(incl_ops[j]);
	SubFitFuncParameterMap &op_params = subfit_pmaps[{incl_ops[i],incl_ops[j]}];
	for(int state=0;state<nstate;state++){
	  op_params[stringize("Asrc%d",state)] = stringize("A%s_%d",istub.c_str(),state);
	  op_params[stringize("Asnk%d",state)] = stringize("A%s_%d",jstub.c_str(),state);
	  op_params[stringize("E%d",state)] = stringize("E%d",state);
	}
	op_params["Csys"] = "C" + istub + "_" + jstub;
      }
    }
  }else if(fitfunc == FitFuncType::FSimGenMultiStateLogEdiff){ 
    //Define the full set of outer parameters
    int p=0;
    for(int state=0;state<nstate;state++){
      std::string statestr = anyToStr(state);
      std::string statem1str = anyToStr(state-1);

      for(int i=0;i<nops;i++) DEF("A" + getStub(incl_ops[i]) + "_" + statestr, 1e6 / sqrt(Ascale)); 
      if(state == 0){ DEF("logE" + statestr, -1.05); }
      else{ DEF("logE" + statestr + "mE" + statem1str, -1.05); }
    }
    for(int i=0;i<nops;i++)
      for(int j=i;j<nops;j++)
	DEF("C" + getStub(incl_ops[i]) + "_" + getStub(incl_ops[j]), 0.);
      
    //Map internal sub-fit parameters to outer params
    for(int i=0;i<nops;i++){
      auto istub = getStub(incl_ops[i]);
      for(int j=i;j<nops;j++){
	auto jstub = getStub(incl_ops[j]);
	SubFitFuncParameterMap &op_params = subfit_pmaps[{incl_ops[i],incl_ops[j]}];
	for(int state=0;state<nstate;state++){
	  op_params[stringize("Asrc%d",state)] = stringize("A%s_%d",istub.c_str(),state);
	  op_params[stringize("Asnk%d",state)] = stringize("A%s_%d",jstub.c_str(),state);
	  if(state == 0) op_params[stringize("logE%d",state)] = stringize("logE%d",state);
	  else op_params[stringize("logE%dmE%d",state,state-1)] = stringize("logE%dmE%d",state,state-1);
	}
	op_params["Csys"] = "C" + istub + "_" + jstub;
      }
    }
  }else if(fitfunc == FitFuncType::FSimGenMultiStateCparam){ 
    //Define the full set of outer parameters
    int p=0;
    for(int state=0;state<nstate;state++){
      std::string statestr = anyToStr(state);
      for(int i=0;i<nops;i++) DEF("A" + getStub(incl_ops[i]) + "_" + statestr, 1e6 / sqrt(Ascale)); 
      DEF("E" + statestr,0.3 * (state+1));
    }
    for(int i=0;i<nops;i++)
      DEF("C" + getStub(incl_ops[i]), 0.);
      
    //Map internal sub-fit parameters to outer params
    for(int i=0;i<nops;i++){
      auto istub = getStub(incl_ops[i]);
      for(int j=i;j<nops;j++){
	auto jstub = getStub(incl_ops[j]);
	SubFitFuncParameterMap &op_params = subfit_pmaps[{incl_ops[i],incl_ops[j]}];
	for(int state=0;state<nstate;state++){
	  op_params[stringize("Asrc%d",state)] = stringize("A%s_%d",istub.c_str(),state);
	  op_params[stringize("Asnk%d",state)] = stringize("A%s_%d",jstub.c_str(),state);
	  op_params[stringize("E%d",state)] = stringize("E%d",state);
	}
	op_params["Csrc"] = "C" + istub;
	op_params["Csnk"] = "C" + jstub;
      }
    }
  }else{
    assert(0);
  }


  //Setup the default guesses
  guess.resize(&param_map); assert(guess.size() == guesses.size());
  for(int i=0;i<guess.size();i++) guess(i) = guesses[i];


#undef DEF
}

void loadGuess(Params &guess, const std::string &guess_file){
  std::map<std::string,double> gmap;
  parse(gmap, guess_file);
  guess.mapImport(gmap);
  std::cout << "Loaded guess: " << guess << std::endl;
}

void saveGuess(const Params &guess, const std::string &guess_file){
  std::map<std::string,double> gmap;
  guess.mapExport(gmap);
  std::ofstream of(guess_file);
  parser_tools::parser_output_print(of, gmap);
}



#endif
