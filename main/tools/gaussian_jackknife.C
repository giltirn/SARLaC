//Generate a jackknife distribution with a chosen central value and error using a random gaussian

#include<common.h>
#include<utils.h>
#include<random.h>

using namespace SARLaC;


jackknifeDistributionD randomJackknife(double cen, double err, double err_tol, int N, int iter_max = 1000, int iter = 0){
  double dev = err / sqrt(double(N-1));
  jackknifeDistributionD out(N);
  for(int i=0;i<N;i++)
    out.sample(i) = gaussianRandom<double>(cen, dev);

  //Shift slightly to correct mean
  double delta = cen - out.mean();
  for(int i=0;i<N;i++)
    out.sample(i) += delta;

  //Make it recursive to get the error correct to within tol
  if( abs(  (out.standardError() - err) / err  ) > err_tol ){
    if(iter >= iter_max){
      throw "Error: Could not generate distribution with desired tolerance";
    }
    return randomJackknife(cen, err, err_tol, N, iter_max, iter+1);
  }
  return out;
}

void writeBootXML(const jackknifeDistributionD &out, const std::string &tag, const std::string &file){
  std::ofstream ff(file);
  ff << std::setprecision(16);
  ff << "<?xml version=\"1.0\"?>\n\n";
  ff << "<data_in_file>\n";
  ff << "  <Nentries>1</Nentries>\n";
  ff << "  <list>\n";
  ff << "    <elem>\n";
  ff << "      <SampleType>SuperJackBoot</SampleType>\n";
  ff << "      <Nmeas>" << out.size() << "</Nmeas>\n";
  ff << "      <Nensembles>1</Nensembles>\n";
  ff << "      <Ensembles>\n";
  ff << "        <elem>\n";
  ff << "          <tag>" << tag << "</tag>\n";
  ff << "          <SampleType>Jackknife</SampleType>\n";
  ff << "          <EnsembleSize>" << out.size() << "</EnsembleSize>\n";
  ff << "          <avg>" << out.best() << "</avg>\n";
  ff << "          <values>";
  for(int s=0;s<out.size();s++)
    ff << out.sample(s) << " ";
  ff << "</values>\n";
  ff << "        </elem>\n";
  ff << "      </Ensembles>\n";
  ff << "    </elem>\n";
  ff << "  </list>\n";
  ff << "</data_in_file>\n";
}



int main(const int argc, const char* argv[]){
  assert(argc >= 5);
  int i=1;
  double cen = strToAny<double>(argv[i++]);
  double std_err = strToAny<double>(argv[i++]);
  int nsample = strToAny<int>(argv[i++]);
  std::string outfile = argv[i++];

  double err_tol = 1e-3;
  int ntry = 10000;
  
  int seed = 4967;

  bool superjack_bootxml = false;
  std::string superjack_bootxml_tag;
  
  while(i<argc){
    std::string arg(argv[i]);
    if(arg == "-cen_tol"){ 
      std::cout << "Warning: cen_tol is deprecated" << std::endl;
      i+=2;
    }else if(arg == "-err_tol"){ 
      err_tol = strToAny<double>(argv[i+1]);
      i+=2;
    }else if(arg == "-ntry"){ 
      ntry = strToAny<int>(argv[i+1]);
      i+=2;
    }else if(arg == "-seed"){ 
      seed = strToAny<int>(argv[i+1]);
      i+=2;
    }else if(arg == "-superjack_bootxml"){
      superjack_bootxml = true;
      superjack_bootxml_tag = argv[i+1];
      i+=2;
    }else{
      error_exit(std::cout << "Unrecognized argument: " << arg << std::endl);
    }
  }
  
  RNG.initialize(seed);

  jackknifeDistributionD out(nsample, cen);
  if(std_err != 0.){
    try{
      out = randomJackknife(cen, std_err, err_tol, nsample, ntry);
    }catch(const char* e){
      std::cerr << "Error: " << e << std::endl;
      return 1;
    }
  }
  std::cout << out << std::endl;

  if(superjack_bootxml)
    writeBootXML(out, superjack_bootxml_tag, outfile);
  else 
    writeParamsStandard(out, outfile);
  return 0;
}
