#ifndef _CPSFIT_DUAL_NUMBER_H_
#define _CPSFIT_DUAL_NUMBER_H_

#include<iostream>
#include<cmath>

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

//A naive implementation for dual numbers (cf https://en.wikipedia.org/wiki/Automatic_differentiation)

struct dual{
  double x;
  double xp;
  
  dual() = default;
  explicit dual(double x, double xp = 0.): x(x), xp(xp){}
  dual(std::initializer_list<double> xxp){ 
    auto it = xxp.begin();
    x = *it;
    ++it;
    xp = *it;
  }
};

inline std::ostream & operator<<(std::ostream &os, const dual &a){
  os << "(" << a.x << ", " << a.xp << ")";
  return os;
}


inline dual operator+(const dual &a, const dual &b){
  return dual(a.x+b.x, a.xp+b.xp);
}
inline dual operator+(const double &a, const dual &b){
  return dual(a+b.x, b.xp);
}
inline dual operator+(const dual &a, const double &b){
  return dual(a.x+b, a.xp);
}

inline dual operator-(const dual &a, const dual &b){
  return dual(a.x-b.x, a.xp-b.xp);
}
inline dual operator-(const double &a, const dual &b){
  return dual(a-b.x, -b.xp);
}
inline dual operator-(const dual &a, const double &b){
  return dual(a.x-b, a.xp);
}


inline dual operator*(const dual &a, const dual &b){
  return dual(a.x * b.x, a.x * b.xp + a.xp * b.x);
}
inline dual operator*(const double &a, const dual &b){
  return dual(a * b.x, a * b.xp);
}
inline dual operator*(const dual &a, const double &b){
  return dual(a.x * b, a.xp * b);
}


inline dual operator/(const dual &a, const dual &b){
  return dual(a.x/b.x, (a.xp * b.x - a.x * b.xp)/b.x/b.x);
}
inline dual operator/(const double &a, const dual &b){
  return dual(a/b.x, - a * b.xp/b.x/b.x);
}
inline dual operator/(const dual &a, const double &b){
  return dual(a.x/b, a.xp/b);
}

using ::exp;

inline dual exp(const dual &a){
  return dual(::exp(a.x), ::exp(a.x)*a.xp );
}

using ::log;

inline dual log(const dual &a){
  return dual(::log(a.x), a.xp/a.x );
}

using ::sin;

inline dual sin(const dual &a){
  return dual(::sin(a.x), a.xp*::cos(a.x));
}

using ::cos;

inline dual cos(const dual &a){
  return dual(::cos(a.x), -a.xp*::sin(a.x));
}

using ::pow;

inline dual pow(const dual &a, const double &b){
  return dual( ::pow(a.x,b), b*::pow(a.x, b-1.)*a.xp );
}

using ::sqrt;

inline dual sqrt(const dual &a){
  return dual( ::sqrt(a.x), 0.5/::sqrt(a.x)*a.xp );
}

inline dual operator-(const dual &a){
  return dual( -a.x, -a.xp );
}

using ::sinh;
using ::cosh;

inline dual sinh(const dual &a){
  return dual(::sinh(a.x), a.xp*::cosh(a.x));
}
inline dual cosh(const dual &a){
  return dual(::cosh(a.x), a.xp*::sinh(a.x));
}

using::asinh;
using::acosh;

inline dual asinh(const dual &a){
  return dual(::asinh(a.x), a.xp/::sqrt(a.x*a.x + 1.)  );
}
inline dual acosh(const dual &a){
  return dual(::acosh(a.x), a.xp/::sqrt(a.x*a.x - 1.)  );
}



CPSFIT_END_NAMESPACE

#endif
