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
  
  complexD operator()(const int i) const{ return z[i] + tau * y[i]; }
};

#endif
