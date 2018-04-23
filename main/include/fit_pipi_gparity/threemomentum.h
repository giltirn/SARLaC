#ifndef _PIPI_THREEMOMENTUM_H_
#define _PIPI_THREEMOMENTUM_H_

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


#endif
