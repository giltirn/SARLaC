#ifndef _FIT_MPI_GPARITY_AMA_READ_DATA_H
#define _FIT_MPI_GPARITY_AMA_READ_DATA_H

typedef NumericSquareMatrix<rawDataDistribution<double> > rawDataDistributionMatrix;
typedef NumericVector<rawDataDistribution<double> > rawDataDistributionVector;

struct configData{
  rawDataDistributionMatrix *all_data;
  int conf;

  configData(rawDataDistributionMatrix &ad, const int c): all_data(&ad), conf(c){}
  inline int getLt() const{ return all_data->size(); }
};
namespace configData_grammar{
  namespace ascii = boost::spirit::x3::ascii;
  namespace x3 = boost::spirit::x3;
  using boost::fusion::at_c;

  struct set_entry{
    template <typename Context> void operator()(Context const& ctx) const{
      const auto &attr = x3::_attr(ctx);
      configData &con = x3::_val(ctx);
      distribution<double> &d = con.all_data->operator()(at_c<0>(attr),at_c<1>(attr));
      d.sample(con.conf) = at_c<2>(attr); //just the real part
    }
  };
  
  auto const line = x3::lexeme[x3::int_ > ' ' > x3::int_ > ' ' > x3::double_ > ' ' > x3::double_];
  x3::rule<struct line_rule_, configData> const line_rule = "line_rule";
  auto const line_rule_def = line[set_entry()];
  BOOST_SPIRIT_DEFINE(line_rule);
};

std::istream & operator>>(std::istream &is, configData &s){
#if 0
  namespace ascii = boost::spirit::x3::ascii;
  namespace x3 = boost::spirit::x3;
  
  is >> std::noskipws;

  boost::spirit::istream_iterator f(is), l;
  const int nelems = s.getLt()*s.getLt();
  bool r = x3::phrase_parse(f, l, x3::repeat(nelems)[configData_grammar::line_rule], ascii::space, s);
  
  if(!r){
    throw std::runtime_error("Parsing ConfigData failed\n");
  }
#else

  const int nelems = s.getLt()*s.getLt();
  int i,j;
  double imag;
  for(int e=0;e<nelems;e++){
    if(!(is >> i >> j)) error_exit(std::cout << "configData failed to read indices for config" << s.conf << "\n");
    double &to = s.all_data->operator()(i,j).sample(s.conf);
    if(!(is >> to >> imag)) error_exit(std::cout << "configData failed to read indices for config " << s.conf << "\n");
  }

#endif
  return is;
}

void read(configData &into, const std::string &fmt, const int traj){
  const std::string file = subsIdx(fmt,traj);
  std::ifstream is(file.c_str());
  if(!is.good()){ std::cout << "Could not open file \"" << file << "\"\n"; std::cout.flush(); exit(-1); }
  is >> into;
  if(is.fail() || is.bad()){ std::cout << "Error reading file \"" << file << "\"\n"; std::cout.flush(); exit(-1); }
  is.close();
}
void read(rawDataDistributionMatrix &exact, rawDataDistributionMatrix &sloppy, const SloppyExact &args, const int traj, const int sample_idx){
  if(args.include_data){
    configData e(exact, sample_idx);
    configData s(sloppy, sample_idx);    
    read(e, args.exact_fmt, traj);
    read(s, args.sloppy_fmt, traj);
  }
}





#endif
