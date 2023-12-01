#ifndef _FIT_PIPI_GND_EXC_GPARITY_FITFUNC_H
#define _FIT_PIPI_GND_EXC_GPARITY_FITFUNC_H

//Define the mappings between the sub-fit parameters and the full set of parameters plus index the parameters
typedef std::unordered_map<std::string, std::string> SubFitFuncParameterMap;
typedef std::unordered_map<std::string,size_t> ParamTagIdxMap;
typedef taggedValueContainer<double,std::string> Params;

void setupParameterMaps(std::map< std::pair<std::string,std::string>, SubFitFuncParameterMap > &subfit_pmaps,
			ParamTagIdxMap &param_map,
			Params &guess,
			FitFuncType fitfunc, const double Ascale){
#define DEF(NM, GUESS) { param_map[NM] = p++;  guesses.push_back(GUESS); }  
  static const std::vector<std::string> op = { "gnd", "exc" };

  std::vector<double> guesses;
  
  if(fitfunc == FitFuncType::FSimGenOneState){
    //Define the full set of parameters and default guesses
    int p=0;
    DEF("Apipi_gnd", 1e6 / sqrt(Ascale));
    DEF("Apipi_exc", 1e6 / sqrt(Ascale));
    DEF("Epipi", 0.35);
    DEF("Cpipi_gnd_gnd", 0);
    DEF("Cpipi_exc_exc", 0);
    DEF("Cpipi_gnd_exc", 0);

    //Map internal sub-fit parameters to outer params
    for(int i=0;i<2;i++){
      for(int j=i;j<2;j++){
	subfit_pmaps[{op[i],op[j]}] =
	  {  {"Asrc","Apipi_" + op[i]}, {"Asnk", "Apipi_" + op[j]}, {"E", "Epipi"},      
	     {"Csys", "Cpipi_" + op[i] + "_" + op[j]} };
      }
    }
  }else if(fitfunc == FitFuncType::FSimGenTwoState){
    //Define the full set of parameters and default guesses
    int p=0;
    for(int i=0;i<2;i++) DEF("Apipi_" + op[i] + "_0", 1e6 / sqrt(Ascale)); 
    DEF("Epipi",0.35);
    for(int i=0;i<2;i++) DEF("Apipi_" + op[i] + "_1", 1e6 / sqrt(Ascale));
    DEF("Eexc", 0.7);
    for(int i=0;i<2;i++)
      for(int j=i;j<2;j++)
	DEF("Cpipi_" + op[i] + "_" + op[j], 0.);

    //Map internal sub-fit parameters to outer params
    for(int i=0;i<2;i++){
      for(int j=i;j<2;j++){
	subfit_pmaps[{op[i],op[j]}] =
	  {  {"Asrc0","Apipi_" + op[i] + "_0"}, {"Asnk0", "Apipi_" + op[j] + "_0"}, {"E0", "Epipi"},
	     {"Asrc1","Apipi_" + op[i] + "_1"}, {"Asnk1", "Apipi_" + op[j] + "_1"}, {"E1", "Eexc"},
	     {"Csys", "Cpipi_" + op[i] + "_" + op[j]} };
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
