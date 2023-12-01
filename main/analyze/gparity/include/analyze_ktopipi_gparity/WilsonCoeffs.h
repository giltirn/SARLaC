#ifndef _ANALYZE_KTOPIPI_WILSON_COEFFS_H_
#define _ANALYZE_KTOPIPI_WILSON_COEFFS_H_

class WilsonCoeffs{
  typedef std::complex<double> complexD;
  std::vector<double> y;
  std::vector<double> z;
  complexD tau;
 public:
  WilsonCoeffs(): y(10),z(10){}
  WilsonCoeffs(const std::string &file): WilsonCoeffs(){ read(file); }
  WilsonCoeffs(const double mu, const PerturbativeVariables &var, const std::complex<double> &_tau): WilsonCoeffs(){ compute(mu,var,_tau); }

  //Expect tau as a pair of doubles on the *first* line
  //After for each Q_i, expect line format   i  z_i  y_i
  void read(const std::string &file, const bool verbose = true){
    std::ifstream ff(file.c_str());
    if (!ff.good()){ std::cout << "WilsonCoeffs::read failed to open file " << file << '\n'; std::cout.flush(); exit(-1); }

    double tr, ti;
    ff >> tr >> ti;
    assert(!ff.bad());

    tau = complexD(tr,ti);

    if(verbose) std::cout << "Read tau = (" << tr << ", " << ti << ")\n";

    for(int i=1;i<=10;i++){
      int ii; ff >> ii; assert(!ff.bad());
      if(ii!=i) error_exit(std::cout << "WilsonCoeffs::read found coeffs for Q_" << ii << " where expected those for Q_" << i << '\n');

      double zi,yi;
      ff >> z[i-1] >> y[i-1];
      assert(!ff.bad());
      if(verbose) std::cout << "Read z[" << i << "] = " << z[i-1] << " and y[" << i << "] = " << y[i-1] << std::endl;
    }
    ff.close();
  }
  
  void compute(const double mu, const PerturbativeVariables &var, const std::complex<double> &_tau){
    NumericVector<double> zz, yy;
    DeltaS1WilsonCoeffs::compute3fWilsonCoeffs(yy,zz,mu,var);
    for(int i=0;i<10;i++){
      y[i] = yy[i];
      z[i] = zz[i];
    }
    tau = _tau;
  }       
  void write(const std::string &filename) const{
    std::ofstream ff(filename.c_str());
    ff << tau.real() << " " << tau.imag() << std::endl;
    for(int i=0;i<10;i++)
      ff << i+1 << " " << z[i] << " " << y[i] << std::endl;
  }
  void writeLatex(const std::string &filename) const{
    std::ofstream ff(filename.c_str());
    ff << "$\\tau = " << tau.real() << " + " << tau.imag() << " i$" << std::endl;
    ff << "\\begin{table}[h]\n\\begin{tabular}{cc}\n\\hline\\hline\n $i$ & $y_i$ & $z_i$ \\\\ \\hline \\\\" << std::endl;
    for(int i=0;i<10;i++)
      ff << i+1 << " & " << y[i] << " & " << z[i] << "\\\\" << std::endl;
    ff << "\\end{tabular}\n\\end{table}" << std::endl;
  }

  complexD operator()(const int i) const{ return z[i] + tau * y[i]; }
  double getZ(const int i) const{ return z[i]; }
  double getY(const int i) const{ return y[i]; }
};

#endif
