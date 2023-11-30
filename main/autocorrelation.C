//Compute the autocorrelation of a jackknife distribution
#include<distribution.h>
#include<fit/fit_wrapper/fit_wrapper_freeze.h>
#include<common.h>
#include<plot.h>

using namespace CPSfit;

std::vector<int> parseIndices(const std::string &ind){
  std::vector<int> indices;
  std::istringstream buffer(ind);
  std::string line;
  while(std::getline(buffer, line, ',')){
    int i;
    std::stringstream ss; ss << line; ss >> i;
    indices.push_back(i);
  }
  return indices;
}

std::vector<double> autoCorrelation(const jackknifeDistributionD &A){
  const int nsample = A.size();
  double Abar = A.mean();
  double Astddev = A.standardError()/sqrt(nsample - 1.);

  const int maxDelta = nsample-1;

  std::vector<double> out(maxDelta+1);
  
  for(int Delta = 0; Delta <= maxDelta; Delta++){
    double v = 0.;
    int N = 0;
    for(int i=0;i<nsample;i++){
      if(i + Delta < nsample){
	v += (A.sample(i) - Abar)*(A.sample(i+Delta) - Abar);
      }
    }
    v = v/Astddev/Astddev/nsample;
    
    out[Delta] = v;
  }
  return out;
}

std::vector<double> integratedAutoCorrelation(const std::vector<double> &C){
  int N = C.size();
  std::vector<double> out(N);
  for(int cut=1; cut < N; cut++){
    double v = 0.5;
    for(int i=1;i<=cut;i++)
      v += C[i];
    out[cut-1] = v;
  }
  return out;
}

std::vector<bootstrapDistributionD> autoCorrelation(const jackknifeDistributionD &A, const int bin_size){
  typename bootstrapDistributionD::initType boot_params;
  boot_params.boots = 500;
  boot_params.confidence = 68;

  const int nsample = A.size();
  const int maxDelta = nsample-1;

  double Abar = A.mean();
  double Astddev = A.standardError()/sqrt(nsample - 1.);

  std::vector<bootstrapDistributionD> out;
  
  std::cout << "Computing autocorrelation function:\n";

  for(int Delta = 0; Delta <= maxDelta; Delta++){
    int ndelta = 0;
    for(int i=0;i<nsample;i++)      
      if(i + Delta < nsample)
	++ndelta;

    if(ndelta < 2*bin_size) continue; //only one binned sample, no error possible

    //Construct  (Y_i - Ybar)(Y_{i+Delta} - Ybar)/sigma^2   for each sample i
    rawDataDistributionD C_samples(ndelta);

    int nn=0;
    for(int i=0;i<nsample;i++)      
      if(i + Delta < nsample){
	double v = (A.sample(i) - Abar)*(A.sample(i+Delta) - Abar)/Astddev/Astddev;
	C_samples.sample(nn++) = v;
      }  
    
    //These are many approximations to the autocorrelation function. We wish to average them. Do this under a boostrap to get a width also
    //Increase bin size until this error stops growing

    rawDataDistributionD C_samples_binned = C_samples.bin(bin_size,true);
    bootstrapDistributionD C_boot(C_samples_binned, boot_params);

    std::cout << Delta << " " << C_boot << std::endl;
    out.push_back(C_boot);
  }
  return out;
}


//tau(cut) = 1/2 + \sum_1^cut corr(t)
std::vector<bootstrapDistributionD> integratedAutoCorrelation(const std::vector<bootstrapDistributionD> &C, int stop = -1){
  int N = stop != -1 ? stop+1 : C.size();
  bootstrapDistributionD sum = C[0];
  for(int i=0;i<sum.size();i++) sum.sample(i) = 0.5;
  sum.best() = 0.5;

  std::vector<bootstrapDistributionD> out(N);
  out[0] = sum;

  for(int t=1;t<N;t++){
    sum = sum + C[t];
    out[t] = sum;
  }
  return out;
}


int main(const int argc, const char** argv){
  assert(argc >= 5);
  std::string file = argv[1];
  std::vector<int> idx = parseIndices(argv[2]);
  int bin_size = strToAny<int>(argv[3]);
  std::string plot_stub = argv[4];

  int stop_sep = -1;
  int traj_inc = 1; //scale trajectory separations by meas frequency

  {
    int i=5;
    while(i<argc){      
      std::string arg = argv[i];
      if(arg == "-stop"){
	stop_sep = strToAny<int>(argv[i+1]);
	i+=2;
      }else if(arg == "-traj_inc"){
	traj_inc = strToAny<int>(argv[i+1]);
	i+=2;
      }else{
	error_exit(std::cout << "Unknown argument: \"" << arg << "\"");
      }
    }
  }

  jackknifeDistributionD A;
  readHDF5file(A,file,idx);

  std::vector<bootstrapDistributionD> C = autoCorrelation(A, bin_size);
  std::vector<bootstrapDistributionD> tau_int = integratedAutoCorrelation(C, stop_sep);
  
  //Scale by traj_inc
  for(int i=1;i<tau_int.size();i++) tau_int[i] = tau_int[i] * traj_inc;

  std::cout << "Integrated autocorrelation length:\n";
  for(int i=0;i<tau_int.size();i++)
    std::cout << i*traj_inc << " " << tau_int[i] << std::endl;


  {
    MatPlotLibScriptGenerate plot;
    
    struct Accessor{
      const std::vector<bootstrapDistributionD> &vec;
      int traj_inc;
      Accessor(const std::vector<bootstrapDistributionD> &vec, int traj_inc): vec(vec), traj_inc(traj_inc){}
      
      inline double x(const int i) const{ return i*traj_inc; }
      inline double upper(const int i) const{ return vec[i].confidenceRegion().second; }
      inline double lower(const int i) const{ return vec[i].confidenceRegion().first; }
      inline int size() const{ return vec.size(); }
    };
    
    Accessor acc(tau_int, traj_inc);
    plot.errorBand(acc);
    plot.setXlabel(R"($\Delta_{\rm cut}$)");
    plot.setYlabel(R"($\tau_{\rm int}(\Delta_{\rm cut}$))");
    plot.write(plot_stub +".py", plot_stub+".pdf");
  }


  return 0;
}

