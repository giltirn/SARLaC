#ifndef _GSL_VECTOR_WRAPPER_H_
#define _GSL_VECTOR_WRAPPER_H_

//A wrapper class around GSL vectors. Note, this does not check vector dimensions during operations, so be careful!

#include<iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

class GSLvector{
  gsl_vector * v;
  
public:
  inline explicit GSLvector(const int d): v(gsl_vector_calloc(d)){  } //Initialized to zero
  inline explicit GSLvector(const int d, const double init): v(gsl_vector_alloc(d)){ gsl_vector_set_all(v,init); }

  inline const int dim() const{ return v->size; }

  inline GSLvector(const GSLvector &r){ 
    v = gsl_vector_alloc(r.dim());
    gsl_vector_memcpy (v,r.v);
  }

  inline const double & operator[](const int i) const{ return *gsl_vector_const_ptr(v,i); }
  inline double & operator[](const int i) { return *gsl_vector_ptr(v,i); }

  inline GSLvector & operator+=(const GSLvector &r){
    gsl_vector_add(v,r.v); return *this;
  }
  inline GSLvector & operator-=(const GSLvector &r){
    gsl_vector_sub(v,r.v); return *this;
  }
  
  //a[i] = a[i]*b[i]
  inline GSLvector & outer_prod(const GSLvector &r){
    gsl_vector_mul(v,r.v); return *this;
  }

  //a[i] = a[i]/b[i]
  inline GSLvector & outer_div(const GSLvector &r){
    gsl_vector_div(v,r.v); return *this;
  }  

  inline GSLvector & operator*=(const double x){
    gsl_vector_scale(v,x); return *this;
  }

  //a[i] = a[i] + x
  inline GSLvector & operator+=(const double x){
    gsl_vector_add_constant(v,x); return *this;
  }

  //Compute the Euclidean norm ||x||_2 = \sqrt {\sum x_i^2} of the vector x. 
  inline double norm() const{ 
    return gsl_blas_dnrm2 (v);
  }
  //Squared norm
  inline double norm2() const{
    double nrm = norm();
    return nrm*nrm;
  }

  inline ~GSLvector(){
    gsl_vector_free(v);
  }

  friend inline double dot(const GSLvector &a, const GSLvector &b);
};
inline std::ostream & operator<<(std::ostream &os, const GSLvector &v){
  os << "(" << v[0] << "," << v[1] << "," << v[2] << ")"; return os;
}


inline double dot(const GSLvector &a, const GSLvector &b){
  double out;
  gsl_blas_ddot (a.v,b.v,&out);
  return out;
}

inline GSLvector operator-(const GSLvector &a, const GSLvector &b){
  GSLvector out(a); out -= b; return out;
}

inline GSLvector operator*(const double x, const GSLvector &v){
  GSLvector out(v); out *= x; return out;
}
inline GSLvector operator*(const GSLvector &v,const double x){
  GSLvector out(v); out *= x; return out;
}


CPSFIT_END_NAMESPACE

#endif
