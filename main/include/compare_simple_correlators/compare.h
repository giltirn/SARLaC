#ifndef _COMPARE_SIMPLE_CORRELATORS_COMPARE_H_H
#define _COMPARE_SIMPLE_CORRELATORS_COMPARE_H_H

void compareRelativeDifferences(const jackknifeCorrelationFunctionD &A, const jackknifeCorrelationFunctionD &B, const std::string plot_stub = "reldiff"){
  assert(A.size() == B.size());
  const int sz = A.size();
  jackknifeCorrelationFunctionD reldiff_j(sz,
					  [&](const int i){
					    assert(A.coord(i) == A.coord(i));
					    jackknifeDistributionD reldiff = 2.*(B.value(i) - A.value(i))/(B.value(i) + A.value(i));
					    return typename jackknifeCorrelationFunctionD::ElementType(A.coord(i), std::move(reldiff));
					  }
					  );
  
  std::cout << "Jackknife relative differences:\n";
  for(int t=0;t<sz;t++){    
    std::cout << reldiff_j.coord(t) << " " << reldiff_j.value(t) << std::endl;
  }

  MatPlotLibScriptGenerate plot;
  typedef DataSeriesAccessor<jackknifeCorrelationFunctionD, ScalarCoordinateAccessor<double>,  DistributionPlotAccessor<jackknifeDistributionD> > accessor;
  accessor a(reldiff_j);
  plot.plotData(a);
  plot.write(plot_stub+".py",plot_stub+".eps");
}

#endif
