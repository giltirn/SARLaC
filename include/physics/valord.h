#ifndef _VALORD_H
#define _VALORD_H

#include<map>
#include<cmath>
#include<iostream>

#include<config.h>
#include<utils/macros.h>
#include<utils/utils/string.h>

CPSFIT_START_NAMESPACE

//A class that contains a power series in 2 variables 'as' and 'a' which are taken as coefficients of alpha_s and alpha_EM in perturbative expansion
struct ValOrd{
  int minas;
  int maxas;
  int maxa;
  int mina;

  std::map<std::pair<int,int>, double> v;

  void zero(){
    v.clear();
    v[{0,0}] = 0;
    minas = maxas = maxa = mina = 0;
  }

  ValOrd(){
    zero();
  }
  ValOrd(const int powas, const int powa, const int val): ValOrd(){
    this->operator()(powas,powa) = val;
  }

  double value() const{
    double out = 0.;
    for(auto it = v.begin(); it != v.end(); ++it) out += it->second;
    return out;
  }

  double & operator()(const int powas, const int powa){
    minas = std::min(powas, minas);
    maxas = std::max(powas, maxas);
    mina = std::min(powa, mina);
    maxa = std::max(powa, maxa);
    return v[{powas,powa}];
  }
  double operator()(const int powas, const int powa) const{
    auto it = v.find({powas,powa});
    return it == v.end() ? 0. : it->second;
  }

  ValOrd &operator=(const double r){
    zero();
    this->operator()(0,0) = r;
    return *this;
  }
   

};


inline ValOrd operator+(const ValOrd &a, const ValOrd &b){
  ValOrd out(a);
  for(int powas = std::min(a.minas,b.minas); powas <= std::max(a.maxas,b.maxas); powas++){
    for(int powa = std::min(a.mina,b.mina); powa <= std::max(a.maxa,b.maxa); powa++){
      out(powas, powa) += b(powas, powa);
    }
  }
  return out;
}
inline ValOrd operator+(const ValOrd &a, const double b){
  ValOrd out(a);
  out(0,0) += b;
  return out;
}
inline ValOrd operator+(const double a, const ValOrd &b){
  ValOrd out(b);
  out(0,0) += a;
  return out;
}

inline ValOrd operator-(const ValOrd &a){
  ValOrd out(a);
  for(std::map<std::pair<int,int>, double>::iterator it = out.v.begin(); it != out.v.end(); it++)
    it->second = -it->second;
  return out;
}


inline ValOrd operator-(const ValOrd &a, const ValOrd &b){
  ValOrd out(a);
  for(int powas = std::min(a.minas,b.minas); powas <= std::max(a.maxas,b.maxas); powas++){
    for(int powa = std::min(a.mina,b.mina); powa <= std::max(a.maxa,b.maxa); powa++){
      out(powas, powa) -= b(powas, powa);
    }
  }
  return out;
}
inline ValOrd operator-(const ValOrd &a, const double b){
  ValOrd out(a);
  out(0,0) -= b;
  return out;
}
inline ValOrd operator-(const double a, const ValOrd &b){
  ValOrd out(-b);
  out(0,0) += a;
  return out;
}




inline ValOrd operator*(const ValOrd &a, const ValOrd &b){
  ValOrd out;
  for(auto ait = a.v.begin(); ait != a.v.end(); ++ait){
    int powasl = ait->first.first;
    int powal = ait->first.second;

    for(auto bit = b.v.begin(); bit != b.v.end(); ++bit){
      int powasr = bit->first.first;
      int powar = bit->first.second;

      int powasout = powasl + powasr;
      int powaout = powal + powar;
      out(powasout, powaout) += ait->second * bit->second;
    }
  }
  return out;
}
inline ValOrd operator*(const ValOrd &a, const double b){
  ValOrd out(a);
  for(std::map<std::pair<int,int>, double>::iterator it = out.v.begin(); it != out.v.end(); it++)
    it->second *= b;
  return out;
}
inline ValOrd operator*(const double a, const ValOrd &b){
  ValOrd out(b);
  for(std::map<std::pair<int,int>, double>::iterator it = out.v.begin(); it != out.v.end(); it++)
    it->second *= a;
  return out;
}


inline ValOrd operator/(const ValOrd &a, const double b){
  ValOrd out(a);
  for(std::map<std::pair<int,int>, double>::iterator it = out.v.begin(); it != out.v.end(); it++)
    it->second /= b;
  return out;
}

ValOrd _AS(1,0,1.);
ValOrd _AE(0,1,1.);
ValOrd _ASM1(-1,0,1.);



std::ostream & operator<<(std::ostream &os, const ValOrd &v){
  bool first = true;
  for(auto it=v.v.begin(); it != v.v.end(); it++){
    if(fabs(it->second) > 1e-10){
      int powas = it->first.first;
      int powa = it->first.second;

      if(!first && it->second > 0.)
	os << "+";
      os << it->second << (powas != 0 ? stringize("as^%d", powas): "") << (powa != 0 ? stringize("a^%d", powa): "");

      first = false;
    }
  }
  if(first) os << 0;

  return os;
}
    
//Keep only terms linear in ae, as and O(0).   ae^2/a counts as ae^2.   Also remove ae as
ValOrd truncateO1(const ValOrd &v){
  ValOrd out(v);
  std::map<std::pair<int,int>, double>::iterator it = out.v.begin();
  while(it != out.v.end()){
    int powas = it->first.first;
    int powa = it->first.second;
      
    //int powtot = powas + powa;
    //if(powtot > 1) it = out.v.erase(it);

    if(powas > 1 || powa > 1 || (powas == 1 && powa==1) ) it = out.v.erase(it);
    else ++it;
  }
  return out;
}

//Expand out in ae.  ae^0 term is kept at all orders of as
//ae^1 term is truncated at as^0

//Keep only terms linear in ae and O(0)
ValOrd truncateO1e(const ValOrd &v){
  ValOrd out(v);
  std::map<std::pair<int,int>, double>::iterator it = out.v.begin();
  while(it != out.v.end()){
    int powas = it->first.first;
    int powa = it->first.second;
    if( (powa ==1 && powas > 0) || powa > 1) it = out.v.erase(it);    
    else ++it;
  }
  return out;
}


CPSFIT_END_NAMESPACE

#endif
