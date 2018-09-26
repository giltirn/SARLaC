#ifndef _FIT_PIPI_GND_EXC_SIGMA_GPARITY_FITFUNC_H
#define _FIT_PIPI_GND_EXC_SIGMA_GPARITY_FITFUNC_H

//Define the mappings between the sub-fit parameters and the full set of parameters plus index the parameters
typedef std::unordered_map<std::string, std::string> SubFitFuncParameterMap;
typedef std::unordered_map<std::string,size_t> ParamTagIdxMap;
typedef taggedValueContainer<double,std::string> Params;

enum Operator { PiPiGnd, PiPiExc, Sigma };

std::ostream & operator<<(std::ostream &os, Operator op){
  switch(op){
  case PiPiGnd:
    os << "PiPiGnd"; return os;
  case PiPiExc:
    os << "PiPiExc"; return os;
  case Sigma:
    os << "Sigma"; return os;  
  };
  assert(0);
}

void setupParameterMaps(std::map< std::pair<Operator,Operator>, SubFitFuncParameterMap > &subfit_pmaps,
			ParamTagIdxMap &param_map,
			Params &guess,
			FitFuncType fitfunc, const double Ascale){
#define DEF(NM, GUESS) { param_map[NM] = p++;  guesses.push_back(GUESS); }  
  static const std::vector<std::string> op_stubs = { "pipi_gnd", "pipi_exc", "sigma" };
  static const std::vector<Operator> ops = {PiPiGnd, PiPiExc, Sigma};

  std::vector<double> guesses;
  
  if(fitfunc == FitFuncType::FSimGenTwoState){ 
    //Define the full set of parameters and default guesses

    int p=0;
    for(int i=0;i<3;i++) DEF("A" + op_stubs[i] + "_0", 1e6 / sqrt(Ascale)); 
    DEF("Epipi",0.35);
    for(int i=0;i<3;i++) DEF("A" + op_stubs[i] + "_1", 1e6 / sqrt(Ascale));
    DEF("Eexc", 0.7);
    for(int i=0;i<3;i++)
      for(int j=i;j<3;j++)
	DEF("C" + op_stubs[i] + "_" + op_stubs[j], 0.);
      

    //Map internal sub-fit parameters to outer params
    for(int i=0;i<3;i++){
      for(int j=i;j<3;j++){
	subfit_pmaps[{ops[i],ops[j]}] =
	  {  {"Asrc0","A" + op_stubs[i] + "_0"}, {"Asnk0", "A" + op_stubs[j] + "_0"}, {"E0", "Epipi"},
	     {"Asrc1","A" + op_stubs[i] + "_1"}, {"Asnk1", "A" + op_stubs[j] + "_1"}, {"E1", "Eexc"},
	     {"Csys", "C" + op_stubs[i] + "_" + op_stubs[j]} };
      }
    }

  }else if(fitfunc == FitFuncType::FSimGenThreeState){ 
    //Define the full set of parameters and default guesses

    int p=0;
    for(int i=0;i<3;i++) DEF("A" + op_stubs[i] + "_0", 1e6 / sqrt(Ascale)); 
    DEF("E0",0.35);
    for(int i=0;i<3;i++) DEF("A" + op_stubs[i] + "_1", 1e6 / sqrt(Ascale));
    DEF("E1", 0.7);
    for(int i=0;i<3;i++) DEF("A" + op_stubs[i] + "_2", 1e6 / sqrt(Ascale));
    DEF("E2", 0.9);
    for(int i=0;i<3;i++)
      for(int j=i;j<3;j++)
	DEF("C" + op_stubs[i] + "_" + op_stubs[j], 0.);
      
    //Map internal sub-fit parameters to outer params
    for(int i=0;i<3;i++){
      for(int j=i;j<3;j++){
	subfit_pmaps[{ops[i],ops[j]}] =
	  {  {"Asrc0","A" + op_stubs[i] + "_0"}, {"Asnk0", "A" + op_stubs[j] + "_0"}, {"E0", "E0"},
	     {"Asrc1","A" + op_stubs[i] + "_1"}, {"Asnk1", "A" + op_stubs[j] + "_1"}, {"E1", "E1"},
	     {"Asrc2","A" + op_stubs[i] + "_2"}, {"Asnk2", "A" + op_stubs[j] + "_2"}, {"E2", "E2"},
	     {"Csys", "C" + op_stubs[i] + "_" + op_stubs[j]} };
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
