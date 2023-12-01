#pragma once

template<typename T>
void parseOrTemplate(T &args, const std::string &params_file, const std::string &template_file){
  if(params_file == "TEMPLATE"){
    std::ofstream of(template_file);
    (std::cout << "Outputting data generation template to '" << template_file << "' and exiting\n").flush();
    of << args;
    of.close();
    exit(0);
  } 
  parse(args, params_file);
}

double estimatePvalue(const double q2, const std::vector<double> &q2s){
  double c=0, n=q2s.size();
  for(double v : q2s)
    if(v > q2) c+=1.;
  return c/n;
}

