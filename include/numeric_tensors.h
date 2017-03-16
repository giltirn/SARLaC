#ifndef _NUMERIC_TENSORS_H
#define _NUMERIC_TENSORS_H

template<typename Numeric>
class NumericVector{
  std::vector<Numeric> v;
public:
  NumericVector():v(){}
  explicit NumericVector(const int n):v(n){}
  NumericVector(const int n, const Numeric &def):v(n,def){}

  int size() const{ return v.size(); }
  
  void resize(const int n){
    v.resize(n);
  }
  void zero(){ for(int i=0;i<v.size();i++) v[i] = 0.; }

  std::string print() const{
    std::ostringstream os;
    os << v[0];
    for(int i=1;i<v.size();i++)
      os << " " << v[i];
    return os.str();
  }

  Numeric & operator()(const int i){ return v[i]; }
  const Numeric & operator()(const int i) const { return v[i]; }

  Numeric & operator[](const int i){ return v[i]; }
  const Numeric & operator[](const int i) const { return v[i]; }
  
  NumericVector &operator+=(const NumericVector &r){
    for(int i=0;i<v.size();i++) v[i] += r.v[i];
  }
};

template<typename Numeric>
class SVDinvertPolicy;



template<typename Numeric, typename InvertPolicy = SVDinvertPolicy<Numeric> >
class NumericMatrix: public InvertPolicy{ //square matrix
  std::vector<std::vector<Numeric> > m;
public:
  NumericMatrix():m(){}
  explicit NumericMatrix(const int n): m(n, std::vector<Numeric>(n)){}

  int size() const{ return m.size(); }
  
  void resize(const int n){
    m.resize(n);
    for(int i=0;i<n;i++)
      m[i].resize(n);      
  }
  void zero(){
    const int n = m.size();
    for(int i=0;i<n;i++)
      for(int j=0;j<n;j++)
	m[i][j] = 0.;
  }
  void invert(const NumericMatrix<Numeric> &what){
    this->InvertPolicy::invert(m,what.m);
  }

  std::string print() const{
    std::ostringstream os;
    for(int i=0;i<m.size();i++){
      os << m[i][0];
      for(int j=1;j<m[i].size();j++)
	os << " " << m[i][j];
      os << std::endl;
    }
    return os.str();
  }
  
  Numeric & operator()(const int i, const int j){ return m[i][j]; }
  const Numeric & operator()(const int i, const int j) const { return m[i][j]; }
};

template<typename Numeric>
class SVDinvertPolicy{
 protected:
  inline static void invert(std::vector<std::vector<Numeric> > &inv_m, const std::vector<std::vector<Numeric> > &m){
    assert(m.size() == m[0].size()); //square matrix
    inv_m.resize(m.size(), std::vector<Numeric>(m.size()));
    svd_inverse(inv_m, m);
  }
};

#endif
