#ifndef _PIPI_STATIONARY_PLOT_H
#define _PIPI_STATIONARY_PLOT_H

void plotDeterminantTest(const std::string &file_stub, const ResampledData<jackknifeCorrelationFunctionD> &data, const std::vector<Operator> &ops, const int Lt){
  struct plotdata{
    std::vector<jackknifeDistributionD> v;
    std::vector<int> t;
    
    double x(const int i) const{ return t[i]; }
    double y(const int i) const{ return v[i].best(); }
    double dxm(const int i) const{ return 0; }
    double dxp(const int i) const{ return 0; }
    double dym(const int i) const{ return v[i].standardError(); }
    double dyp(const int i) const{ return v[i].standardError(); }
    int size() const{ return v.size(); }
  };
  plotdata pdata;
  
  int n = ops.size();
  for(int t=0;t<Lt;t++){
    NumericSquareMatrix<jackknifeDistributionD> C(n);
    for(int i=0;i<n;i++)
      for(int j=i;j<n;j++)
	C(i,j) = C(j,i) = data.correlator(ops[i],ops[j]).value(t);
    jackknifeDistributionD det = determinant(C);
    for(int i=0;i<n;i++) det = det / C(i,i);
    
    pdata.t.push_back(t);
    pdata.v.push_back(det);
  }
  
  MatPlotLibScriptGenerate plot;
  plot.plotData(pdata);
  plot.write(file_stub+".py", file_stub+".pdf");
}


#endif
