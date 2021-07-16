#ifndef TENSOR_H
#define TENSOR_H "$Header:$"
//----------------------------  tensor_base.h  ---------------------------
//    $Id: tensor_base.h,v 1.43 2003/01/03 22:21:26 deal Exp $
//    Version: $Name: Version-4-0-0 $
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tensor_base.h  ---------------------------
//----------------------------  tensor.h  ---------------------------
//    $Id: tensor.h,v 1.47 2003/01/03 22:21:26 deal Exp $
//    Version: $Name: Version-4-0-0 $
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tensor.h  ---------------------------
//
//modified from deal.II/base/include/base/tensor_base.h and tensor.h
//  by John Trangenstein, August 2009
//**********************************************************************
// Copyright 2009 John A. Trangenstein
//
// This software is made available for research and instructional use 
// only. 
// You may copy and use this software without charge for these 
// non-commercial purposes, provided that the copyright notice and 
// associated text is reproduced on all copies.  
// For all other uses (including distribution of modified versions), 
// please contact the author at
//   John A. Trangenstein
//   Department of Mathematics
//   Duke University
//   Durham, NC 27708-0320
//   USA
// or
//   johnt@math.duke.edu
// 
// This software is made available "as is" without any assurance that it
// is completely correct, or that it will work for your purposes.  
// Use the software at your own risk.
//**********************************************************************

#include <ostream>
#include <Vector.H>
#include <Point.H>
#include <Tracer.H>
#include <Types.H>

template <int rank, int dim> class Tensor;
//compare to DEAL Tensor
template <int dim> class Tensor<1,dim> {
  private:
    Point<dim> values;
    inline void unrollRecursion(LAVector &v,int &i) const
    {
      for (unsigned j=0;j<dim;j++) v[i++]=values[j];
    }
    friend class Tensor<2,dim>;
    template <int d> friend 
      ostream& operator<<(ostream&,const Tensor<1,d>&);
  public:
    typedef double array_type[(dim!=0) ? dim : INT_MAX];
    inline explicit Tensor(bool initialize=true) {
      if (initialize) for (int i=0;i<dim;++i) values[i]=0.;
    }
    inline Tensor(const array_type &t) {
      for (int i=0;i<dim;++i) values[i]=t[i];
    }
    inline explicit Tensor(const Point<dim> &p) : values(p) { }
    inline Tensor(const Tensor<1,dim> &t) : values(t.values) { }
    const Point<dim>& getValues() const { return values; }
    Point<dim>& getValues() { return values; }
    inline double operator[](int i) const { return values[i]; }

    inline double& operator[](int i) { return values[i]; }

    inline Tensor<1,dim>& operator=(const Point<dim> &p) { 
      values=p; return *this; 
    }
    inline Tensor<1,dim>& operator=(const Tensor<1,dim> &t) { 
      values=t.values; return *this; 
    }
    inline bool operator==(const Tensor<1,dim> &t) const { 
      return values==t.values; 
    }
    inline bool operator!=(const Tensor<1,dim> &t) const {
      return !(*this==t);
    }
    inline Tensor<1,dim>& operator+=(const Tensor<1,dim> &t) {
      values+=t.values; return *this;
    }
    inline Tensor<1,dim>& operator-=(const Tensor<1,dim> &t) {
      values-=t.values; return *this;
    }
    inline Tensor<1,dim>& operator*=(double x) {values*=x; return *this; }
    inline Tensor<1,dim>& operator/=(double x) {values/=x; return *this; }
    inline double operator*(const Tensor<1,dim> &t) const { 
      return values*t.values; 
    }
    inline Tensor<1,dim> operator+(const Tensor<1,dim> &t) const {
      Tensor<1,dim> t2(*this); return t2+=t; 
    }
    inline Tensor<1,dim> operator-(const Tensor<1,dim> &t) const {
      Tensor<1,dim> t2(*this); return t2-=t; 
    }
    inline Tensor<1,dim> operator*(double x) const {
      Tensor<1,dim> t2(*this); return t2*=x; 
    }
    inline Tensor<1,dim> operator/(double x) const {
      Tensor<1,dim> t2(*this); return t2/=x; 
    }
    inline Tensor<1,dim> operator-() const {
      Tensor<1,dim> t2(-values); return t2;
    }
    void initialize(double r) { values.initialize(r); } 
    void normalize() {
      double len=sqrt(values*values);
      CHECK_POSITIVE(len);
      values/=len;
    }
    inline void unroll(LAVector &v) const {
      CHECK_SAME(v.size(),dim);
      int i=0; unrollRecursion(v,i);
    }
};

template <int dim> inline
ostream& operator<<(ostream &os,const Tensor<1,dim> &t) {
  os << "[" << t.values << "]";
  return os;
}

template <int dim> inline void product(double &t,
const Tensor<1,dim> &t1,const Tensor<1,dim> &t2) {
  t=0.; for (int i=0;i<dim;i++) t+=t1[i]*t2[i];
}

template <int dim> inline void product(Tensor<1,dim> &t2,
const Tensor<1,dim> &t,double x) {
  t2=t; t2*=x;
}

template <int dim> inline void product(Tensor<1,dim> &t2,double x,
const Tensor<1,dim> &t) {
  t2=t; t2*=x;
}

template <int dim> inline void quotient(Tensor<1,dim> &t2,
const Tensor<1,dim> &t,double x) {
  t2=t; t2/=x;
}

//compare to DEAL Tensor
template <int rank, int dim> class Tensor {
  private:
    Tensor<rank-1,dim> subtensor[dim];
    void unrollRecursion(LAVector &v,int &si) const {
      for (unsigned i=0;i<dim;++i) operator[](i).unrollRecursion(v,si);
    }
    friend class Tensor<rank+1,dim>;
  public:
//  the following basically says that array_type is double[dim]...[dim]
    typedef typename Tensor<rank-1,dim>::array_type array_type[dim];

    inline explicit Tensor(bool initialize=true) {
      if (initialize) for (int i=0;i<dim;++i) subtensor[i].initialize(0.);
    }
    Tensor(const array_type &a) {
      for (int i=0;i<dim;++i) {
        subtensor[i]=Tensor<rank-1,dim>(a[i]);
      }
    }
    Tensor<rank-1,dim>& operator [](int i) {
      CHECK_BOUNDS(i,0,dim);
      return subtensor[i];
    }
    const Tensor<rank-1,dim>& operator[](int i) const {
      CHECK_BOUNDS(i,0,dim);
      return subtensor[i];
    }

    Tensor& operator=(const Tensor<rank,dim> &t) {
      for (int i=0;i<dim;++i) subtensor[i]=t.subtensor[i];
      return *this;
    }
    bool operator==(const Tensor<rank,dim> &p) const {
      for (int i=0;i<dim;++i) {
        if (subtensor[i]!=p.subtensor[i]) return false;
      }
      return true;
    }
    bool operator!=(const Tensor<rank,dim> &p) const {
      return !((*this) == p);
    }
    Tensor<rank,dim>& operator+=(const Tensor<rank,dim> &p) {
      for (int i=0;i<dim;++i) subtensor[i]+=p.subtensor[i];
      return *this;
    }
    Tensor<rank,dim>& operator-=(const Tensor<rank,dim> &p) {
      for (int i=0;i<dim;++i) subtensor[i]-=p.subtensor[i];
      return *this;
    }
    Tensor<rank,dim>& operator*=(double s) {
      for (int i=0;i<dim;++i) subtensor[i]*=s; return *this;
    }
    Tensor<rank,dim>& operator/=(double s) {
      for (int i=0;i<dim;++i) subtensor[i]/=s; return *this;
    }
    Tensor<rank,dim> operator+(const Tensor<rank,dim> &t) const {
      Tensor<rank,dim> tmp(*this); tmp+=t; return tmp;
    }
    Tensor<rank,dim> operator-(const Tensor<rank,dim> &t) const {
      Tensor<rank,dim> tmp(*this); tmp-=t; return tmp;
    }
    Tensor<rank,dim> operator-() const {
      Tensor<rank,dim> tmp(*this); 
      for (int i=0;i<dim;++i) tmp.subtensor[i]=-subtensor[i];
      return tmp;
    }
    void unroll(LAVector &v) const {
      CHECK_SAME(v.size(),pow(static_cast<double>(dim),rank));
      int i=0; unrollRecursion(v,i);
    }
    void initialize(double r) {
      for (int i=0; i<dim; ++i) subtensor[i].initialize(r);
    }
};

template <int rank, int dim> inline ostream& 
operator<<(ostream &os,const Tensor<rank,dim> &p) {
  os << "[";
  for (int i=0; i<dim-1; ++i) os << p[i] << ' ';
  os << p[dim-1] << "]";
  return os;
}

template <int rank, int dim> inline void product(Tensor<rank,dim> &tt,
const Tensor<rank,dim> &t,double s) {
  tt=t; tt*=s;
}

template <int rank, int dim> inline void product(Tensor<rank,dim> &tt,
double s,const Tensor<rank,dim> &t) {
  tt=t; tt*=s;
}

template <int rank, int dim> inline void quotient(Tensor<rank,dim> &tt,
const Tensor<rank,dim> &t,double s) {
  tt=t; tt/=s;
}

template <int dim> inline void product(Tensor<1,dim> &t,
const Tensor<2,dim> &t2,const Tensor<1,dim> &t1) {
  t.initialize(0.);
  for (int i=0;i<dim;++i) {
    t[i]=t2[i]*t1;
  }
}

//not in DEAL
template <int dim> inline void transposeProduct(Tensor<1,dim> &t,
const Tensor<2,dim> &t2,const Tensor<1,dim> &t1) {
  t.initialize(0.);
  for (int i=0;i<dim;++i) {
    for (int j=0;j<dim;++j) t[i]+=t2[j][i]*t1[j];
  }
}

template <int dim> inline void product(Tensor<1,dim> &t,
const Tensor<1,dim> &t1, const Tensor<2,dim> &t2) {
  t.initialize(0.);
  for (int i=0;i<dim;++i) {
    for (int j=0;j<dim;++j) t[i]+=t1[j]*t2[j][i];
  }
}

template <int dim> inline void product(Tensor<2,dim> &t,
const Tensor<2,dim> &t1,const Tensor<2,dim> &t2) {
  t.initialize(0.);
  for (int i=0;i<dim;++i) {
    for (int j=0;j<dim;++j) {
      for (int k=0;k<dim;++k) {
        t[i][j]+=t1[i][k]*t2[k][j];
      }
    }
  }
}

template <int dim> inline void product(Tensor<2,dim> &t,
const Tensor<2,dim> &t1,int i1,
const Tensor<2,dim> &t2,int i2) {
  CHECK_BOUNDS(i1,1,3);
  CHECK_BOUNDS(i2,1,3);
  t.initialize(0.);
  if (i1==1) {
    if (i2==1) {
      for (int i=0;i<dim;++i) {
        for (int j=0;j<dim;++j) {
          for (int k=0;k<dim;++k) {
            t[i][j]+=t1[k][i]*t2[k][j];
          }
        }
      }
    } else {
      for (int i=0;i<dim;++i) {
        for (int j=0;j<dim;++j) {
          for (int k=0;k<dim;++k) {
            t[i][j]+=t1[k][i]*t2[j][k];
          }
        }
      }
    }
  } else {
    if (i2==1) {
      for (int i=0;i<dim;++i) {
        for (int j=0;j<dim;++j) {
          for (int k=0;k<dim;++k) {
            t[i][j]+=t1[i][k]*t2[k][j];
          }
        }
      }
    } else {
      for (int i=0;i<dim;++i) {
        for (int j=0;j<dim;++j) {
          for (int k=0;k<dim;++k) {
            t[i][j]+=t1[i][k]*t2[j][k];
          }
        }
      }
    }
  }
}

template <int dim> inline void product(Tensor<2,dim> &t,
const Tensor<3,dim> &t3,int i,const Tensor<1,dim> &t1) {
  CHECK_BOUNDS(i,1,4);
  t.initialize(0.);
  switch (i) {
    case 1:
      for (int i=0;i<dim;++i) {
        for (int j=0;j<dim;++j) {
          for (int k=0;k<dim;++k) {
            t[i][j]+=t3[k][i][j]*t1[k];
          }
        }
      }
      break;
    case 2:
      for (int i=0;i<dim;++i) {
        for (int j=0;j<dim;++j) {
          for (int k=0;k<dim;++k) {
            t[i][j]+=t3[i][k][j]*t1[k];
          }
        }
      }
      break;
    case 3:
      for (int i=0;i<dim;++i) {
        for (int j=0;j<dim;++j) {
          for (int k=0;k<dim;++k) t[i][j]+=t3[i][j][k]*t1[k];
        }
      }
      break;
    default:
      break;
  };
}

template <int dim> inline void product(Tensor<3,dim> &t,
const Tensor<3,dim> &t3,const Tensor<2,dim> &t2) {
  t.initialize(0.);
  for (int i=0;i<dim;++i) {
    for (int j=0;j<dim;++j) {
      for (int k=0;k<dim;++k) {
        for (int l=0;l<dim;++l) {
          t[i][j][k]+=t3[i][j][l]*t2[l][k];
        }
      }
    }
  }
}

template <int dim> inline void product(Tensor<3,dim> &t,
const Tensor<2,dim> &t2,const Tensor<3,dim> &t3) {
  t.initialize(0.);
  for (int i=0;i<dim;++i) {
    for (int j=0;j<dim;++j) {
      for (int k=0;k<dim;++k) {
        for (int l=0;l<dim;++l) {
          t[i][j][k]+=t2[i][l]*t3[l][j][k];
        }
      }
    }
  }
}

template <int dim> inline void product(Tensor<4,dim> &t,
const Tensor<3,dim> &t1,const Tensor<3,dim> &t2){
  t.initialize(0.);
  for (int i=0;i<dim;++i) {
    for (int j=0;j<dim;++j) {
      for (int k=0;k<dim;++k) {
        for (int l=0;l<dim;++l) {
          for (int m=0;m<dim;++m) {
            t[i][j][k][l]+=t1[i][j][m]*t2[m][k][l];
          }
        }
      }
    }
  }
}

template <int dim> inline void outerProduct(Tensor<2,dim> &t,
const Tensor<1,dim> &t1,const Tensor<1,dim> &t2) {
  for (int i=0;i<dim;++i) {
    for (int j=0;j<dim;++j) t[i][j]=t1[i]*t2[j];
  }
}

template <int dim> inline void outerProduct(Tensor<3,dim> &t,
const Tensor<1,dim> &t1,const Tensor<2,dim> &t2) {
  for (int i=0;i<dim;++i) {
    for (int j=0;j<dim;++j) {
      for (int k=0;k<dim;++k) t[i][j][k]=t1[i]*t2[j][k];
    }
  }
}

template <int dim> inline void outerProduct(Tensor<3,dim> &t,
const Tensor<2,dim> &t1,const Tensor<1,dim> &t2) {
  for (int i=0;i<dim;++i) {
    for (int j=0;j<dim;++j) {
      for (int k=0;k<dim;++k) {
        t[i][j][k]=t1[i][j]*t2[k];
      }
    }
  }
}

template <int dim> inline void outerProduct(Tensor<1,dim> &t,double s,
const Tensor<1,dim> &t1) {
  product(t,s,t1);
}

template <int dim> inline void outerProduct(Tensor<1,dim> &t,
const Tensor<1,dim> t1,double s) {
  product(t,t1,s);
}

// rotate counter-clockwise by 90 degrees:
inline void crossProduct(Tensor<1,2> &t,const Tensor<1,2> &t1) {
  t[0]=t1[1];
  t[1]=-t1[0];
}

inline void crossProduct(Tensor<1,3> &t,const Tensor<1,3> &t1,
const Tensor<1,3> &t2) {
  t[0]=t1[1]*t2[2]-t1[2]*t2[1];
  t[1]=t1[2]*t2[0]-t1[0]*t2[2];
  t[2]=t1[0]*t2[1]-t1[1]*t2[0];
}

template <int rank> inline double determinant(const Tensor<rank,1> &t) {
  return determinant(t[0]);
}

inline double determinant(const Tensor<1,1> &t) { return t[0]; }

inline double determinant(const Tensor<2,2> &t) {
  return ((t[0][0]*t[1][1])-(t[1][0]*t[0][1]));
}

inline double determinant (const Tensor<2,3> &t) {
// get this using Maple:
// with(linalg);
// a := matrix(3,3);
// x := det(a);
// readlib(C);
// C(x, optimized);
  return ( t[0][0]*t[1][1]*t[2][2]
          -t[0][0]*t[1][2]*t[2][1]
          -t[1][0]*t[0][1]*t[2][2]
          +t[1][0]*t[0][2]*t[2][1]
          +t[2][0]*t[0][1]*t[1][2]
          -t[2][0]*t[0][2]*t[1][1] );
}

//d det(A)/d A = A^{-T} det(A)
//not in DEAL
inline Tensor<2,1> determinantDerivative(const Tensor<2,1> &array) {
  Tensor<2,1> derivative;
  derivative.initialize(1.);
  return derivative;
}
inline Tensor<2,2> determinantDerivative(const Tensor<2,2> &array) {
  Tensor<2,2> derivative;
  derivative[0][0]=array[1][1];
  derivative[1][0]=-array[0][1];
  derivative[0][1]=-array[1][0];
  derivative[1][1]=array[0][0];
  return derivative;
}
inline Tensor<2,3> determinantDerivative(const Tensor<2,3> &array) {
  Tensor<2,3> derivative;
  derivative[0][0]=array[1][1]*array[2][2]-array[1][2]*array[2][1];
  derivative[1][1]=array[2][2]*array[0][0]-array[2][0]*array[0][2];
  derivative[2][2]=array[0][0]*array[1][1]-array[0][1]*array[1][0];
  derivative[1][0]=array[0][2]*array[2][1]-array[0][1]*array[2][2];
  derivative[2][0]=array[0][1]*array[1][2]-array[0][2]*array[1][1];
  derivative[0][1]=array[1][2]*array[2][0]-array[1][0]*array[2][2];
  derivative[2][1]=array[1][0]*array[0][2]-array[1][2]*array[0][0];
  derivative[0][2]=array[2][1]*array[1][0]-array[2][0]*array[1][1];
  derivative[1][2]=array[2][0]*array[0][1]-array[2][1]*array[0][0];
  return derivative;
}

template <int dim> inline double trace(const Tensor<2,dim> &t) {
  double tr=0.;
  for (int i=0;i<dim;++i) tr+=t[i][i];
  return tr;
}
//not in DEAL
inline Tensor<2,3> traceDerivative(const Tensor<2,3> &array) {
  Tensor<2,3> derivative;
  derivative.initialize(0.);
  for (int i=0;i<3;i++) derivative[i][i]=1.;
  return derivative;
}

//not in DEAL
inline double secondInvariant(const Tensor<2,3> &t) {
  return t[1][1]*t[2][2]+t[2][2]*t[3][3]+t[3][3]*t[1][1]
        -t[1][2]*t[2][1]-t[2][3]*t[3][2]-t[3][1]*t[1][3];
}
//d II_2(A) / d A = I trace(A) - A^{T}
inline Tensor<2,3> secondInvariantDerivative(const Tensor<2,3> &array) {
  Tensor<2,3> derivative;
  derivative[0][0]=array[1][1]+array[2][2];
  derivative[1][1]=array[2][2]+array[0][0];
  derivative[2][2]=array[0][0]+array[1][1];
  derivative[1][0]=-array[0][1];
  derivative[2][0]=-array[0][2];
  derivative[0][1]=-array[1][0];
  derivative[2][1]=-array[1][2];
  derivative[0][2]=-array[2][0];
  derivative[1][2]=-array[2][1];
  return derivative;
}

template <int dim> inline Tensor<2,dim> invert(const Tensor<2,dim> &t) {
  CHECK_BOUNDS(dim,1,4);
  Tensor<2,dim> inverse;
  switch (dim) {
    case 1:
      inverse[0][0]=1.0/t[0][0];
      break;
    case 2: {
      const double det_inv=1.0/determinant(t);
      inverse[0][0]=t[1][1]*det_inv;
      inverse[0][1]=-t[0][1]*det_inv;
      inverse[1][0]=-t[1][0]*det_inv;
      inverse[1][1]=t[0][0]*det_inv;
      break;
    };
    case 3: {
      const double t4=t[0][0]*t[1][1],
                 t6=t[0][0]*t[1][2],
                 t8=t[0][1]*t[1][0],
                t00=t[0][2]*t[1][0],
                t01=t[0][1]*t[2][0],
                t04=t[0][2]*t[2][0],
                t07=1.0/(t4 *t[2][2]-t6 *t[2][1]-t8 *t[2][2]
                        +t00*t[2][1]+t01*t[1][2]-t04*t[1][1]);
      inverse[0][0]=(t[1][1]*t[2][2]-t[1][2]*t[2][1])*t07;
      inverse[0][1]=-(t[0][1]*t[2][2]-t[0][2]*t[2][1])*t07;
      inverse[0][2]=-(-t[0][1]*t[1][2]+t[0][2]*t[1][1])*t07;
      inverse[1][0]=-(t[1][0]*t[2][2]-t[1][2]*t[2][0])*t07;
      inverse[1][1]=(t[0][0]*t[2][2]-t04)*t07;
      inverse[1][2]=-(t6-t00)*t07;
      inverse[2][0]=-(-t[1][0]*t[2][1]+t[1][1]*t[2][0])*t07;
      inverse[2][1]=-(t[0][0]*t[2][1]-t01)*t07;
      inverse[2][2]=(t4-t8)*t07;
      break;
    };
    default:
      break;
  };    
  return inverse;
}

//inverse_transpose(m) * n * det(m)
//not in DEAL
template <int dim> inline void adjointProduct(Tensor<1,dim> &product,
const Tensor<2,dim> &m,const Tensor<1,dim> &n) {
  CHECK_BOUNDS(dim,1,4);
  switch (dim) {
    case 1:
      product[0]=n[0];
      break;
    case 2: {
      product[0]=m[1][1]*n[0]-m[1][0]*n[1];
      product[1]=m[0][0]*n[1]-m[0][1]*n[0];
      break;
    };
    case 3: {
      product[0]=(m[1][1]*m[2][2]-m[2][1]*m[1][2])*n[0]
                -(m[1][0]*m[2][2]-m[2][0]*m[1][2])*n[1]
                +(m[1][0]*m[2][1]-m[2][0]*m[1][1])*n[2];
      product[1]=-(m[0][1]*m[2][2]-m[2][1]*m[0][2])*n[0]
                 +(m[0][0]*m[2][2]-m[2][0]*m[0][2])*n[1]
                 -(m[0][0]*m[2][1]-m[2][0]*m[0][1])*n[2];
      product[2]=(m[0][1]*m[1][2]-m[1][1]*m[0][2])*n[0]
                -(m[0][0]*m[1][2]-m[1][0]*m[0][2])*n[1]
                +(m[0][0]*m[1][1]-m[1][0]*m[0][1])*n[2];
      break;
    };
    default:
      break;
  };    
}

template <int dim> inline void transpose(Tensor<2,dim> &tt,
const Tensor<2,dim> &t) {
  for (int i=0;i<dim;++i) {
    for (int j=0;j<dim;++j) tt[j][i]=t[i][j];
  }
}

inline void transpose(Tensor<2,1> &tt,const Tensor<2,1> &t) { 
  tt=t;
}

inline void transpose(Tensor<2,2> &tt,const Tensor<2,2> &t) {
  tt[0][0]=t[0][0]; tt[0][1]=t[1][0];
  tt[1][0]=t[0][1]; tt[1][1]=t[1][1];
}

inline void transpose(Tensor<2,3> &tt,const Tensor<2,3> &t) {
  tt[0][0]=t[0][0]; tt[0][1]=t[1][0]; tt[0][2]=t[2][0];
  tt[1][0]=t[0][1]; tt[1][1]=t[1][1]; tt[1][2]=t[2][1];
  tt[2][0]=t[0][2]; tt[2][1]=t[1][2]; tt[2][2]=t[2][2];
}

#endif
