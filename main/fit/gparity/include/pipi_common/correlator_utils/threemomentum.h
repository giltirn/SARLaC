#ifndef _PIPI_THREEMOMENTUM_H_
#define _PIPI_THREEMOMENTUM_H_

#include <boost/functional/hash.hpp>

#include<config.h>
#include<utils/macros.h>

SARLAC_START_NAMESPACE

typedef std::array<int,3> threeMomentum;
typedef std::pair<threeMomentum, threeMomentum> sinkSourceMomenta;

std::ostream & operator<<(std::ostream &os, const threeMomentum &mom){
  os << "(" << mom[0] << ", " << mom[1] << ", " << mom[2] << ")";
  return os;
}
std::ostream & operator<<(std::ostream &os, const sinkSourceMomenta &mom){
  os << "Snk:" << mom.first << " Src:" << mom.second;
  return os;
}


inline threeMomentum operator-(const threeMomentum &p){
  return threeMomentum({-p[0],-p[1],-p[2]});
}
inline threeMomentum operator*(const threeMomentum &p, const int r){
  return threeMomentum({p[0]*r,p[1]*r,p[2]*r});
}
inline threeMomentum operator-(const threeMomentum &a, const threeMomentum &b){
  return threeMomentum({a[0]-b[0],a[1]-b[1],a[2]-b[2]});
}
inline threeMomentum operator+(const threeMomentum &a, const threeMomentum &b){
  return threeMomentum({a[0]+b[0],a[1]+b[1],a[2]+b[2]});
}



inline std::string momStr(const threeMomentum &p){
  std::ostringstream os;
  for(int i=0;i<3;i++)
    os << (p[i] < 0 ? "_" : "") << abs(p[i]);
  return os.str();
}
inline sinkSourceMomenta momComb(const int snkx, const int snky, const int snkz,
				 const int srcx, const int srcy, const int srcz){
  return sinkSourceMomenta( {snkx,snky,snkz}, {srcx,srcy,srcz} );
}
inline sinkSourceMomenta momComb(const threeMomentum &snk, const threeMomentum &src){
  return sinkSourceMomenta(snk,src);
}

//(abc) (acb) (bac) (bca) (cab) (cba)
inline threeMomentum axisPerm(const int i, const threeMomentum &p){
  const int a=p[0], b=p[1], c=p[2];

  switch(i){
  case 0:
    return p;
  case 1:
    return threeMomentum({a,c,b});
  case 2:
    return threeMomentum({b,a,c});
  case 3:
    return threeMomentum({b,c,a});
  case 4:
    return threeMomentum({c,a,b});
  case 5:
    return threeMomentum({c,b,a});
  default:
    assert(0);
  }
}

inline threeMomentum cyclicPermute(const int n, const threeMomentum &p){
  return threeMomentum({ p[n % 3], p[ (1+n) % 3], p[ (2+n) % 3 ]});
}

SARLAC_END_NAMESPACE

namespace boost{
  inline std::size_t hash_value(const SARLaC::threeMomentum &p){
  std::size_t seed = 0;
  boost::hash_combine(seed, p[0]);
  boost::hash_combine(seed, p[1]);
  boost::hash_combine(seed, p[2]);
  return seed;
}
};

#endif
