#ifndef _GLOBALFIT_DATA_PARAMS_H_
#define _GLOBALFIT_DATA_PARAMS_H_

GENERATE_ENUM_AND_PARSER(DataType, (mpi2)(mK2)(mOmega) );

#define DATAPARAMS_MEMBERS			\
  (int, lattice)				\
  (DataType, type)				\
  (double, mx)					\
  (double, my)					\
  (double, ml)					\
  (double, mh)


struct DataParams{
  GENERATE_MEMBERS(DATAPARAMS_MEMBERS)

  DataParams(int lattice, DataType type, double mx, double my, double ml, double mh): lattice(lattice), type(type), mx(mx), my(my), ml(ml), mh(mh){}

  DataParams() = default;
  DataParams(const DataParams &r) = default;

  DataParams(const DataParams &base, const double mres): DataParams(base){
    mx += mres;
    my += mres;
    ml += mres;
    mh += mres;
  }
};

GENERATE_PARSER(DataParams, DATAPARAMS_MEMBERS)


superJackknifeDistribution<DataParams> operator+(const DataParams &base, const superJackknifeDistribution<double> &mres){
  superJackknifeDistribution<DataParams> out(mres.getLayout(), DataParams(base, mres.best()));
  for(int s=0;s<out.size();s++) out.sample(s) = DataParams(base, mres.sample(s));
  return out;
}

typedef correlationFunction<superJackknifeDistribution<DataParams> , superJackknifeDistribution<double> > DataSeriesType;

#endif
