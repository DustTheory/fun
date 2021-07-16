// "$Header:$"
//----------------------  tensor_product_polynomials.cc  ------------
//$Id: tensor_product_polynomials.cc,v 1.20 2003/04/21 16:11:54 wolf Exp $
//Version: $Name: Version-4-0-0 $
//
//Copyright (C) 2000, 2001, 2002, 2003 by the deal.II authors
//
//This file is subject to QPL and may not be  distributed
//without copyright and license information. Please refer
//to the file deal.II/doc/license.html for the  text  and
//further information on this license.
//
//----------------------  tensor_product_polynomials.cc  ------------
//
//modified from deal.II/base/source/tensor_product_polynomials.cc
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

#include <ShapeFunction.H>
#include <Array.H>
#include <Shape.H>
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<int dim> void ShapeFunction<dim>::printOn(ostream &os) const {
  os << "ShapeFunction: dim = " << dim << endl;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<int dim> void TensorProductPolynomials<dim>::computeValues(
const Point<dim> &p,NumPtr<REAL> &v) const {
  v.allocate(n_tensor_pols);
  if (dim==1) {
    polynomials[0]->values(p[0],v);
  } else {
    NumPtr<NumPtr<REAL> > lpv(dim);
    for (int d=0;d<dim;d++) {
      lpv[d].allocate(n_pols[d]);
      polynomials[d]->values(p[d],lpv[d]);
    }
    int iv=0;
    if (dim==2) {
      for (int j=0;j<n_pols[1];j++) {
        REAL &lpvj=lpv[1][j];
        for (int i=0;i<n_pols[0];i++) {
          v[iv++]=lpv[0][i]*lpvj;
        }
      }
    } else {
      for (int k=0;k<n_pols[2];k++) {
        REAL &lpvk=lpv[2][k];
        for (int j=0;j<n_pols[1];j++) {
          REAL lpvjk=lpv[1][j]*lpvk;
          for (int i=0;i<n_pols[0];i++) {
            v[iv++]=lpv[0][i]*lpvjk;
          }
        }
      }
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<int dim> void TensorProductPolynomials<dim>::computeGrads(
const Point<dim> &p,NumPtr<Tensor<1,dim> > &g) const {
  g.allocate(n_tensor_pols);
  if (dim==1) {
    NumPtr<REAL> s(n_pols[0]);
    polynomials[0]->slopes(p[0],s);
    for (int i=0;i<n_pols[0];i++) g[i][0]=s[i];
  } else {
    NumPtr<NumPtr<REAL> > lpv(dim);
    NumPtr<NumPtr<REAL> > lps(dim);
    for (int d=0;d<dim;d++) {
      lpv[d].allocate(n_pols[d]);
      lps[d].allocate(n_pols[d]);
      polynomials[d]->values(p[d],lpv[d]);
      polynomials[d]->slopes(p[d],lps[d]);
    }
    int iv=0;
    if (dim==2) {
      for (int j=0;j<n_pols[1];j++) {
        REAL &lpvj=lpv[1][j];
        REAL &lpsj=lps[1][j];
        for (int i=0;i<n_pols[0];i++) {
          Tensor<1,dim> &gv=g[iv++];
          gv[0]=lps[0][i]*lpvj;
          gv[1]=lpv[0][i]*lpsj;
        }
      }
    } else {
      for (int k=0;k<n_pols[2];k++) {
        REAL &lpvk=lpv[2][k];
        REAL &lpsk=lps[2][k];
        for (int j=0;j<n_pols[1];j++) {
          REAL &lpvj=lpv[1][j];
          REAL &lpsj=lps[1][j];
          REAL lpvv=lpvj*lpvk;
          REAL lpsv=lpsj*lpvk;
          REAL lpvs=lpvj*lpsk;
          for (int i=0;i<n_pols[0];i++) {
            Tensor<1,dim> &gv=g[iv++];
            REAL &lpvi=lpv[0][i];
            gv[0]=lps[0][i]*lpvv;
            gv[1]=lpvi*lpsv;
            gv[2]=lpvi*lpvs;
          }
        }
      }
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<int dim> void TensorProductPolynomials<dim>::computeGradGrads(
const Point<dim> &p,NumPtr<Tensor<2,dim> > &gg) const {
  gg.allocate(n_tensor_pols);
  if (dim==1) {
    NumPtr<REAL> s2(n_pols[0]);
    polynomials[0]->slope2s(p[0],s2);
    for (int i=0;i<n_pols[0];i++) gg[i][0][0]=s2[i];
  } else {
    NumPtr<NumPtr<REAL> > lpv(dim);
    NumPtr<NumPtr<REAL> > lps(dim);
    NumPtr<NumPtr<REAL> > lps2(dim);
    for (int d=0;d<dim;d++) {
      lpv[d].allocate(n_pols[d]);
      lps[d].allocate(n_pols[d]);
      lps2[d].allocate(n_pols[d]);
      polynomials[d]->values(p[d],lpv[d]);
      polynomials[d]->slopes(p[d],lps[d]);
      polynomials[d]->slope2s(p[d],lps2[d]);
    }
    int iv=0;
    if (dim==2) {
      for (int j=0;j<n_pols[1];j++) {
        REAL &lpvj=lpv[1][j];
        REAL &lpsj=lps[1][j];
        REAL &lps2j=lps2[1][j];
        for (int i=0;i<n_pols[0];i++) {
          Tensor<2,dim> &ggv=gg[iv++];
          ggv[0][0]=lps2[0][i]*lpvj;
          ggv[1][0]=lps[0][i]*lpsj;
          ggv[1][1]=lpv[0][i]*lps2j;
          ggv[0][1]=ggv[1][0];
        }
      }
    } else {
      for (int k=0;k<n_pols[2];k++) {
        REAL &lpvk=lpv[2][k];
        REAL &lpsk=lps[2][k];
        REAL &lps2k=lps2[2][k];
        for (int j=0;j<n_pols[1];j++) {
          REAL &lpvj=lpv[1][j];
          REAL &lpsj=lps[1][j];
          REAL &lps2j=lps2[1][j];
          REAL lpvv=lpvj*lpvk;
          REAL lpsv=lpsj*lpvk;
          REAL lpvs=lpvj*lpsk;
          REAL lps2v=lps2j*lpvk;
          REAL lpss=lpsj*lpsk;
          REAL lpvs2=lpvj*lps2k;
          for (int i=0;i<n_pols[0];i++) {
            Tensor<2,dim> &ggv=gg[iv++];
            ggv[0][0]=lps2[0][i]*lpvv;
            ggv[1][0]=lps[0][i]*lpsv;
            ggv[2][0]=lps[0][i]*lpvs;
            ggv[1][1]=lpv[0][i]*lps2v;
            ggv[2][1]=lpv[0][i]*lpss;
            ggv[2][2]=lpv[0][i]*lpvs2;
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
          }
        }
      }
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<int dim> void TensorProductPolynomials<dim>::printOn(
ostream &os) const {
  os << "TensorProductPolynomials: n_tensor_pols = " << n_tensor_pols
     << "\n\t n_pols =";
  for (int d=0;d<dim;d++) cout << " " << n_pols[d];
  cout << endl;
  for (int d=0;d<dim;d++) {
    os << "\tpolynomials[" << d << "]:" << endl;
    polynomials[d]->printOn(os);
  }
  ShapeFunction<dim>::printOn(os);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
AnisotropicPolynomials::AnisotropicPolynomials(
const NumPtr<const Polynomial*> &pols,const int (&np)[SPACEDIM]) : 
polynomials(SPACEDIM) {
  CHECK_SAME(pols.getNumber(),SPACEDIM);
  polynomials.copyFrom(pols);
  n_tensor_pols=1;
  for (int d=0;d<SPACEDIM;d++) {
    n_pols[d]=np[d];
    n_tensor_pols*=np[d];
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
void AnisotropicPolynomials::computeValues(const Point<SPACEDIM> &p,
NumPtr<REAL> &v) const {
  v.allocate(n_tensor_pols);
  NumPtr<NumPtr<REAL> > lpv(SPACEDIM);
  for (int d=0;d<SPACEDIM;d++) {
    lpv[d].allocate(n_pols[d]);
    polynomials[d]->values(p[d],lpv[d]);
  }
  int iv=0;
# if (SPACEDIM==2)
  for (int j=0;j<n_pols[1];j++) {
    REAL &lpvj=lpv[1][j];
    for (int i=0;i<n_pols[0];i++) {
      v[iv++]=lpv[0][i]*lpvj;
    }
  }
# else
  for (int k=0;k<n_pols[2];k++) {
    REAL &lpvk=lpv[2][k];
    for (int j=0;j<n_pols[1];j++) {
      REAL lpvjk=lpv[1][j]*lpvk;
      for (int i=0;i<n_pols[0];i++) {
        v[iv++]=lpv[0][i]*lpvjk;
      }
    }
  }
# endif
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
void AnisotropicPolynomials::computeGrads(const Point<SPACEDIM> &p,
NumPtr<Tensor<1,SPACEDIM> > &g) const {
  g.allocate(n_tensor_pols);
  NumPtr<NumPtr<REAL> > lpv(SPACEDIM);
  NumPtr<NumPtr<REAL> > lps(SPACEDIM);
  for (int d=0;d<SPACEDIM;d++) {
    lpv[d].allocate(n_pols[d]);
    lps[d].allocate(n_pols[d]);
    polynomials[d]->values(p[d],lpv[d]);
    polynomials[d]->slopes(p[d],lps[d]);
  }
  int iv=0;
#if (SPACEDIM==2)
  for (int j=0;j<n_pols[1];j++) {
    REAL &lpvj=lpv[1][j];
    REAL &lpsj=lps[1][j];
    for (int i=0;i<n_pols[0];i++) {
      Tensor<1,SPACEDIM> &gv=g[iv++];
      gv[0]=lps[0][i]*lpvj;
      gv[1]=lpv[0][i]*lpsj;
    }
  }
#else
  for (int k=0;k<n_pols[2];k++) {
    REAL &lpvk=lpv[2][k];
    REAL &lpsk=lps[2][k];
    for (int j=0;j<n_pols[1];j++) {
      REAL &lpvj=lpv[1][j];
      REAL &lpsj=lps[1][j];
      REAL lpvv=lpvj*lpvk;
      REAL lpsv=lpsj*lpvk;
      REAL lpvs=lpvj*lpsk;
      for (int i=0;i<n_pols[0];i++) {
        Tensor<1,SPACEDIM> &gv=g[iv++];
        REAL &lpvi=lpv[0][k];
        gv[0]=lps[0][i]*lpvv;
        gv[1]=lpvi*lpsv;
        gv[2]=lpvi*lpvs;
      }
    }
  }
#endif
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
void AnisotropicPolynomials::computeGradGrads(const Point<SPACEDIM> &p,
NumPtr<Tensor<2,SPACEDIM> > &gg) const {
  gg.allocate(n_tensor_pols);
  NumPtr<NumPtr<REAL> > lpv(SPACEDIM);
  NumPtr<NumPtr<REAL> > lps(SPACEDIM);
  NumPtr<NumPtr<REAL> > lps2(SPACEDIM);
  for (int d=0;d<SPACEDIM;d++) {
    lpv[d].allocate(n_pols[d]);
    lps[d].allocate(n_pols[d]);
    lps2[d].allocate(n_pols[d]);
    polynomials[d]->values(p[d],lpv[d]);
    polynomials[d]->slopes(p[d],lps[d]);
    polynomials[d]->slope2s(p[d],lps2[d]);
  }
  int iv=0;
#if (SPACEDIM==2)
  for (int j=0;j<n_pols[1];j++) {
    REAL &lpvj=lpv[1][j];
    REAL &lpsj=lps[1][j];
    REAL &lps2j=lps2[1][j];
    for (int i=0;i<n_pols[0];i++) {
      Tensor<2,SPACEDIM> &ggv=gg[iv++];
      ggv[0][0]=lps2[0][i]*lpvj;
      ggv[1][0]=lps[0][i]*lpsj;
      ggv[1][1]=lpv[0][i]*lps2j;
      ggv[0][1]=ggv[1][0];
    }
  }
#else
  for (int k=0;k<n_pols[2];k++) {
    REAL &lpvk=lpv[2][k];
    REAL &lpsk=lps[2][k];
    REAL &lps2k=lps2[2][k];
    for (int j=0;j<n_pols[1];j++) {
      REAL &lpvj=lpv[1][j];
      REAL &lpsj=lps[1][j];
      REAL &lps2j=lps2[1][j];
      REAL lpvv=lpvj*lpvk;
      REAL lpsv=lpsj*lpvk;
      REAL lpvs=lpvj*lpsk;
      REAL lps2v=lps2j*lpvk;
      REAL lpss=lpsj*lpsk;
      REAL lpvs2=lpvj*lps2k;
      for (int i=0;i<n_pols[0];i++) {
        Tensor<2,SPACEDIM> &ggv=gg[iv++];
        ggv[0][0]=lps2[0][i]*lpvv;
        ggv[1][0]=lps[0][i]*lpsv;
        ggv[2][0]=lps[0][i]*lpvs;
        ggv[1][1]=lpv[0][i]*lps2v;
        ggv[2][1]=lpv[0][i]*lpss;
        ggv[2][2]=lpv[0][i]*lpvs2;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    }
  }
#endif
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
void AnisotropicPolynomials::printOn(ostream &os) const {
  os << "\tAnisotropicPolynomials: n_tensor_pols = " << n_tensor_pols
     << endl;
  for (int d=0;d<SPACEDIM;d++) {
    os << "\tn_pols[" << d << "] = " << n_pols[d] << endl;
  }
  for (int i=0;i<polynomials.getNumber();i++) {
    os << "\tpolynomials[" << i << "] : " << endl;
    polynomials[i]->printOn(os);
  }
  ShapeFunction<SPACEDIM>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
void TrianglePolynomials::printOn(ostream &os) const {
  os << "TrianglePolynomials: order = " << order << endl;
  ShapeFunction<2>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
//these may be a bit more linearly independent
template<int dim> void HomogeneousPolynomials<dim>::computeValues(
const Point<dim> &p,NumPtr<REAL> &v) const {
  v.allocate(n_pols);
  NumPtr<NumPtr<REAL> > lpv(dim);
  for (int d=0;d<dim;d++) {
    NumPtr<REAL> &lpvd=lpv[d];
    REAL pd=p[d];
    lpvd.allocate(order+1);
    REAL pp=1.;
    for (int m=0;m<=order;m++) {
      lpvd[m]=pp;
      pp*=pd;
    }
  }
  NumPtr<REAL> &lpv0=lpv[0];
  NumPtr<REAL> &lpv1=lpv[1];
  int iv=0;
  if (dim==2) {
    for (int i=0;i<=order;i++) {
      v[iv++]=lpv0[i]*lpv1[order-i];
    }
  } else {
    NumPtr<REAL> &lpv2=lpv[2];
    for (int j=0;j<=order;j++) {
      REAL lpvj=lpv1[j];
      for (int i=0;i+j<=order;i++) {
        v[iv++]=lpv0[i]*lpvj*lpv2[order-i-j];
      }
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<int dim> void HomogeneousPolynomials<dim>::computeGrads(
const Point<dim> &p,NumPtr<Tensor<1,dim> > &g) const {
  g.allocate(n_pols);
  NumPtr<NumPtr<REAL> > lpv(dim);
  NumPtr<NumPtr<REAL> > lps(dim);
  for (int d=0;d<dim;d++) {
    NumPtr<REAL> &lpvd=lpv[d];
    NumPtr<REAL> &lpsd=lps[d];
    REAL pd=p[d];
    lpvd.allocate(order+1);
    lpsd.allocate(order+1);
    REAL pp=1.;
    for (int m=0;m<=order;m++) {
      lpvd[m]=pp;
      pp*=pd;
    }
    lpsd[0]=0.;
    for (int m=1;m<=order;m++) {
      lpsd[m]=static_cast<REAL>(m)*lpvd[m-1];
    }
  }
  NumPtr<REAL> &lpv0=lpv[0];
  NumPtr<REAL> &lps0=lps[0];
  NumPtr<REAL> &lpv1=lpv[1];
  NumPtr<REAL> &lps1=lps[1];
  int iv=0;
  if (dim==2) {
    for (int i=0;i<=order;i++) {
      Tensor<1,dim> &gv=g[iv++];
      gv[0]=lps0[i]*lpv1[order-i];
      gv[1]=lpv0[i]*lps1[order-i];
    }
  } else {
    NumPtr<REAL> &lpv2=lpv[2];
    NumPtr<REAL> &lps2=lps[2];
    for (int j=0;j<=order;j++) {
      REAL &lpvj=lpv1[j];
      REAL &lpsj=lps1[j];
      for (int i=0;i+j<=order;i++) {
        Tensor<1,dim> &gv=g[iv++];
        gv[0]=lps0[i]*lpvj*lpv2[order-i-j];
        gv[1]=lpv0[i]*lpsj*lpv2[order-i-j];
        gv[2]=lpv0[i]*lpvj*lps2[order-i-j];
      }
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<int dim> void HomogeneousPolynomials<dim>::computeGradGrads(
const Point<dim> &p,NumPtr<Tensor<2,dim> > &gg) const {
  gg.allocate(n_pols);
  NumPtr<NumPtr<REAL> > lpv(dim);
  NumPtr<NumPtr<REAL> > lps(dim);
  NumPtr<NumPtr<REAL> > lpss(dim);
  for (int d=0;d<dim;d++) {
    NumPtr<REAL> &lpvd=lpv[d];
    NumPtr<REAL> &lpsd=lps[d];
    NumPtr<REAL> &lpssd=lpss[d];
    REAL pd=p[d];
    lpvd.allocate(order+1);
    lpsd.allocate(order+1);
    lpssd.allocate(order+1);
    REAL pp=1.;
    for (int m=0;m<=order;m++) {
      lpvd[m]=pp;
      pp*=pd;
    }
    lpsd[0]=0.;
    for (int m=1;m<=order;m++) {
      lpsd[m]=static_cast<REAL>(m)*lpvd[m-1];
    }
    lpssd[0]=0.;
    for (int m=1;m<=order;m++) {
      lpssd[m]=static_cast<REAL>(m)*lpsd[m-1];
    }
  }
  NumPtr<REAL> &lpv0=lpv[0];
  NumPtr<REAL> &lps0=lps[0];
  NumPtr<REAL> &lpss0=lpss[0];
  NumPtr<REAL> &lpv1=lpv[1];
  NumPtr<REAL> &lps1=lps[1];
  NumPtr<REAL> &lpss1=lpss[1];
  int iv=0;
  if (dim==2) {
    for (int i=0;i<=order;i++) {
      Tensor<2,dim> &ggv=gg[iv++];
      ggv[0][0]=lpss0[i]*lpv1[order-i];
      ggv[1][0]=lps0[i]*lps1[order-i];
      ggv[1][1]=lpv0[i]*lpss1[order-i];
      ggv[0][1]=ggv[1][0];
    }
  } else {
    NumPtr<REAL> &lpv2=lpv[2];
    NumPtr<REAL> &lps2=lps[2];
    NumPtr<REAL> &lpss2=lpss[2];
    for (int j=0;j<=order;j++) {
      REAL &lpvj=lpv1[j];
      REAL &lpsj=lps1[j];
      REAL &lpssj=lpss1[j];
      for (int i=0;i+j<=order;i++) {
        Tensor<2,dim> &ggv=gg[iv++];
        ggv[0][0]=lpss0[i]*lpvj*lpv2[order-i-j];
        ggv[1][0]=lps0[i]*lpsj*lpv2[order-i-j];
        ggv[2][0]=lps0[i]*lpvj*lps2[order-i-j];
        ggv[1][1]=lpv0[i]*lpssj*lpv2[order-i-j];
        ggv[2][1]=lpv0[i]*lpsj*lps2[order-i-j];
        ggv[2][2]=lpv0[i]*lpvj*lpss2[order-i-j];
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<int dim> void HomogeneousPolynomials<dim>::printOn(ostream &os)
const {
  os << "HomogeneousPolynomials: order = " << order << endl;
  os << "\tn_pols = " << n_pols << endl;
  ShapeFunction<dim>::printOn(os);
}
#endif
#if (SPACEDIM>1)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
HomogeneousLagrangePolynomials2::HomogeneousLagrangePolynomials2(
int degree) : order(degree),n_pols(degree+1),zeros(degree+1) {
  if (order>=0) zeros[0]=0.;
  if (order>0) {
    REAL ro=1./static_cast<REAL>(order);
    REAL z=ro;
    for (int i=1;i<order;i++,z+=ro) zeros[i]=z;
    zeros[order]=1.;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void hnumden(const NumPtr<REAL> &zeros,int n,REAL arg1,REAL arg2,
REAL z,REAL &prod) {
  prod=1.;
  REAL den=1.;
  for (int m=0;m<n;m++) {
    const REAL &zm=zeros[m];
    prod*=arg1-zm*(arg1+arg2);
    den*=z-zm;
  }
  prod/=den;
}
void hnumdens(const NumPtr<REAL> &zeros,int n,REAL arg1,REAL arg2,
REAL z,REAL &prod,REAL &dprodda1,REAL &dprodda2) {
  prod=1.;
  dprodda1=0.;
  dprodda2=0.;
  REAL den=1.;
  for (int m=0;m<n;m++) {
    const REAL &zm=zeros[m];
    prod*=arg1-zm*(arg1+arg2);
    den*=z-zm;
    REAL p=1.;
    for (int l=0;l<n;l++) if (l!=m) p*=arg1-zeros[l]*(arg1+arg2);
    dprodda1+=p*(1.-zm);
    dprodda2-=p*zm;
  }
  prod/=den;
  dprodda1/=den;
  dprodda2/=den;
}
void hnumdens2(const NumPtr<REAL> &zeros,int n,REAL arg1,REAL arg2,
REAL z,REAL &prod,REAL &dprodda1,REAL &dprodda2,REAL &d2prodda1da1,
REAL &d2prodda1da2,REAL &d2prodda2da2) {
  prod=1.;
  dprodda1=0.;
  dprodda2=0.;
  d2prodda1da1=0.;
  d2prodda1da2=0.;
  d2prodda2da2=0.;
  REAL den=1.;
  for (int m=0;m<n;m++) {
    const REAL &zm=zeros[m];
    prod*=arg1-zm*(arg1+arg2);
    den*=z-zm;
    REAL p=1.;
    for (int l=0;l<n;l++) {
      if (l==m) continue;
      const REAL &zl=zeros[l];
      p*=arg1-zl*(arg1+arg2);
      REAL p2=1.;
      for (int k=0;k<n;k++) {
        if (k!=m && k!=l) p2*=arg1-zeros[k]*(arg1+arg2);
      }
      d2prodda1da1+=p2*(1.-zm)*(1.-zl);
      d2prodda1da2-=p2*(1.-zm)*zl;
      d2prodda2da2+=p2*zm*zl;
    }
    dprodda1+=p*(1.-zm);
    dprodda2-=p*zm;
  }
  prod/=den;
  dprodda1/=den;
  dprodda2/=den;
  d2prodda1da1/=den;
  d2prodda1da2/=den;
  d2prodda2da2/=den;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void HomogeneousLagrangePolynomials2::computeValues(const Point<2> &p,
NumPtr<REAL> &v) const {
  v.allocate(order+1);
  for (int i=0;i<=order;i++) {
    REAL pi,pnmi;
    hnumden(zeros,i,p[0],p[1],zeros[i],pi);
    hnumden(zeros,order-i,p[1],p[0],zeros[order-i],pnmi);
    v[i]=pi*pnmi;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void HomogeneousLagrangePolynomials2::computeGrads(const Point<2> &p,
NumPtr<Tensor<1,2> > &g) const {
  g.allocate(order+1);
  for (int i=0;i<=order;i++) {
    REAL pi,dpida1,dpida2;
    REAL pnmi,dpnmida1,dpnmida2;
    hnumdens(zeros,i,p[0],p[1],zeros[i],pi,dpida1,dpida2);
    hnumdens(zeros,order-i,p[1],p[0],zeros[order-i],pnmi,dpnmida1,
      dpnmida2);
    g[i][0]=dpida1*pnmi+pi*dpnmida2;
    g[i][1]=dpida2*pnmi+pi*dpnmida1;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void HomogeneousLagrangePolynomials2::computeGradGrads(
const Point<2> &p,NumPtr<Tensor<2,2> > &gg) const {
  gg.allocate(order+1);
  for (int i=0;i<=order;i++) {
    REAL pi,dpida1,dpida2,dpida1da1,dpida1da2,dpida2da2;
    hnumdens2(zeros,i,p[0],p[1],zeros[i],pi,dpida1,dpida2,dpida1da1,
      dpida1da2,dpida2da2);
    int j=order-i;
    REAL pj,dpjda1,dpjda2,dpjda1da1,dpjda1da2,dpjda2da2;
    hnumdens2(zeros,j,p[1],p[0],zeros[j],pj,dpjda1,
      dpjda2,dpjda1da1,dpjda1da2,dpjda2da2);
    gg[i][0][0]=dpida1da1*pj+2.*dpida1*dpjda2+pi*dpjda2da2;
    gg[i][0][1]=dpida1da2*pj+dpida1*dpjda1+dpida2*dpjda2
               + pi*dpjda1da2;
    gg[i][1][1]=dpida2da2*pj+2.*dpida2*dpjda1+pi*dpjda1da1;
    gg[i][1][0]=gg[i][0][1];
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void HomogeneousLagrangePolynomials2::printOn(ostream &os) const {
  os << "HomogeneousLagrangePolynomials2 : order = " << order
     << "\n\tn_pols = " << n_pols << endl;
  for (int i=0;i<zeros.getNumber();i++) {
    cout << "\tzeros[" << i << "] = " << zeros[i] << endl;
  }
  ShapeFunction<2>::printOn(os);
}
#endif
#if (SPACEDIM==3)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
HomogeneousLagrangePolynomials3::HomogeneousLagrangePolynomials3(
int degree) : order(degree),n_pols(((degree+1)*(degree+2))/2),
zeros(degree+1) {
  if (order>=0) zeros[0]=0.;
  if (order>0) {
    REAL ro=1./static_cast<REAL>(order);
    REAL z=ro;
    for (int i=1;i<order;i++,z+=ro) zeros[i]=z;
    zeros[order]=1.;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void HomogeneousLagrangePolynomials3::computeValues(const Point<3> &p,
NumPtr<REAL> &v) const {
  v.allocate(n_pols);
  NumPtr<NumPtr<REAL> > prod(3);
  for (int d=0;d<3;d++) {
    NumPtr<REAL> &pd=prod[d];
    REAL arg2=p[(d+1)%3]+p[(d+2)%3];
    pd.allocate(order+1);
    for (int i=0;i<=order;i++) {
      hnumden(zeros,i,p[d],arg2,zeros[i],pd[i]);
    }
  }
  int iv=0;
  NumPtr<REAL> &p0=prod[0];
  NumPtr<REAL> &p1=prod[1];
  NumPtr<REAL> &p2=prod[2];
  for (int j=0;j<=order;j++) {
    REAL &p1j=p1[j];
    for (int i=0;i+j<=order;i++) {
      v[iv++]=p0[i]*p1j*p2[order-i-j];
    }
  }
  for (int d=0;d<3;d++) {
    prod[d].cleanup();
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void HomogeneousLagrangePolynomials3::computeGrads(const Point<3> &p,
NumPtr<Tensor<1,3> > &g) const {
  g.allocate(n_pols);
  NumPtr<NumPtr<REAL> > prod(3);
  NumPtr<NumPtr<REAL> > dprodda1(3);
  NumPtr<NumPtr<REAL> > dprodda2(3);
  for (int d=0;d<3;d++) {
    NumPtr<REAL> &pd=prod[d];
    NumPtr<REAL> &dpda1=dprodda1[d];
    NumPtr<REAL> &dpda2=dprodda2[d];
    pd.allocate(order+1);
    dpda1.allocate(order+1);
    dpda2.allocate(order+1);
    REAL arg2=p[(d+1)%3]+p[(d+2)%3];
    for (int i=0;i<=order;i++) {
      hnumdens(zeros,i,p[d],arg2,zeros[i],pd[i],dpda1[i],dpda2[i]);
    }
  }
  NumPtr<REAL> &p0=prod[0];
  NumPtr<REAL> &p1=prod[1];
  NumPtr<REAL> &p2=prod[2];
  NumPtr<REAL> &dprod0da1=dprodda1[0];
  NumPtr<REAL> &dprod1da1=dprodda1[1];
  NumPtr<REAL> &dprod2da1=dprodda1[2];
  NumPtr<REAL> &dprod0da2=dprodda2[0];
  NumPtr<REAL> &dprod1da2=dprodda2[1];
  NumPtr<REAL> &dprod2da2=dprodda2[2];
  int iv=0;
  for (int j=0;j<=order;j++) {
    REAL &p1j=p1[j];
    REAL &dp1jda1=dprod1da1[j];
    REAL &dp1jda2=dprod1da2[j];
    for (int i=0;i+j<=order;i++) {
      int k=order-i-j;
      REAL &p0i=p0[i];
      REAL &p2k=p2[k];
      REAL &dp0ida1=dprod0da1[i];
      REAL &dp2kda1=dprod2da1[k];
      REAL &dp0ida2=dprod0da2[i];
      REAL &dp2kda2=dprod2da2[k];
      Tensor<1,3> &gv=g[iv++];
      gv[0]=dp0ida1*p1j*p2k+p0i*dp1jda2*p2k+p0i*p1j*dp2kda2;
      gv[1]=dp0ida2*p1j*p2k+p0i*dp1jda1*p2k+p0i*p1j*dp2kda2;
      gv[2]=dp0ida2*p1j*p2k+p0i*dp1jda2*p2k+p0i*p1j*dp2kda1;
    }
  }
  for (int d=0;d<3;d++) {
    prod[d].cleanup();
    dprodda1[d].cleanup();
    dprodda2[d].cleanup();
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void HomogeneousLagrangePolynomials3::computeGradGrads(
const Point<3> &p,NumPtr<Tensor<2,3> > &gg) const {
  gg.allocate(n_pols);
  NumPtr<NumPtr<REAL> > prod(3);
  NumPtr<NumPtr<REAL> > dprodda1(3);
  NumPtr<NumPtr<REAL> > dprodda2(3);
  NumPtr<NumPtr<REAL> > d2prodda1da1(3);
  NumPtr<NumPtr<REAL> > d2prodda1da2(3);
  NumPtr<NumPtr<REAL> > d2prodda2da2(3);
  for (int d=0;d<3;d++) {
    NumPtr<REAL> &pd=prod[d];
    NumPtr<REAL> &dpda1=dprodda1[d];
    NumPtr<REAL> &dpda2=dprodda2[d];
    NumPtr<REAL> &d2pda1da1=d2prodda1da1[d];
    NumPtr<REAL> &d2pda1da2=d2prodda1da2[d];
    NumPtr<REAL> &d2pda2da2=d2prodda2da2[d];
    pd.allocate(order+1);
    dpda1.allocate(order+1);
    dpda2.allocate(order+1);
    d2pda1da1.allocate(order+1);
    d2pda1da2.allocate(order+1);
    d2pda2da2.allocate(order+1);
    REAL arg2=p[(d+1)%3]+p[(d+2)%3];
    for (int i=0;i<=order;i++) {
      hnumdens2(zeros,i,p[d],arg2,zeros[i],pd[i],dpda1[i],dpda2[i],
        d2pda1da1[i],d2pda1da2[i],d2pda2da2[i]);
    }
  }
  NumPtr<REAL> &p0=prod[0];
  NumPtr<REAL> &p1=prod[1];
  NumPtr<REAL> &p2=prod[2];
  NumPtr<REAL> &dprod0da1=dprodda1[0];
  NumPtr<REAL> &dprod1da1=dprodda1[1];
  NumPtr<REAL> &dprod2da1=dprodda1[2];
  NumPtr<REAL> &dprod0da2=dprodda2[0];
  NumPtr<REAL> &dprod1da2=dprodda2[1];
  NumPtr<REAL> &dprod2da2=dprodda2[2];
  NumPtr<REAL> &d2prod0da1da1=d2prodda1da1[0];
  NumPtr<REAL> &d2prod1da1da1=d2prodda1da1[1];
  NumPtr<REAL> &d2prod2da1da1=d2prodda1da1[2];
  NumPtr<REAL> &d2prod0da1da2=d2prodda1da2[0];
  NumPtr<REAL> &d2prod1da1da2=d2prodda1da2[1];
  NumPtr<REAL> &d2prod2da1da2=d2prodda1da2[2];
  NumPtr<REAL> &d2prod0da2da2=d2prodda2da2[0];
  NumPtr<REAL> &d2prod1da2da2=d2prodda2da2[1];
  NumPtr<REAL> &d2prod2da2da2=d2prodda2da2[2];
  int iv=0;
  for (int j=0;j<=order;j++) {
    REAL &p1j=p1[j];
    REAL &dp1jda1=dprod1da1[j];
    REAL &dp1jda2=dprod1da2[j];
    REAL &d2p1jda1da1=d2prod1da1da1[j];
    REAL &d2p1jda1da2=d2prod1da1da2[j];
    REAL &d2p1jda2da2=d2prod1da2da2[j];
    for (int i=0;i+j<=order;i++) {
      int k=order-i-j;
      REAL &p0i=p0[i];
      REAL &p2k=p2[k];
      REAL &dp0ida1=dprod0da1[i];
      REAL &dp2kda1=dprod2da1[k];
      REAL &dp0ida2=dprod0da2[i];
      REAL &dp2kda2=dprod2da2[k];
      REAL &d2p0ida1da1=d2prod0da1da1[i];
      REAL &d2p2kda1da1=d2prod2da1da1[k];
      REAL &d2p0ida1da2=d2prod0da1da2[i];
      REAL &d2p2kda1da2=d2prod2da1da2[k];
      REAL &d2p0ida2da2=d2prod0da2da2[i];
      REAL &d2p2kda2da2=d2prod2da2da2[k];
      Tensor<2,3> &ggv=gg[iv++];
      ggv[0][0]=
        d2p0ida1da1*p1j*p2k+p0i*d2p1jda2da2*p2k+p0i*p1j*d2p2kda2da2
        +2.*(dp0ida1*dp1jda2*p2k+p0i*dp1jda2*dp2kda2+dp0ida1*p1j*dp2kda2);
      ggv[0][1]=
        d2p0ida1da2*p1j*p2k+p0i*d2p1jda1da2*p2k+p0i*p1j*d2p2kda2da2
        +(dp0ida1*dp1jda1+dp0ida2*dp1jda2)*p2k
        +(dp1jda2*dp2kda2+dp1jda1*dp2kda2)*p0i
        +(dp2kda2*dp0ida2+dp2kda2*dp0ida1)*p1j;
      ggv[0][2]=
        d2p0ida1da2*p1j*p2k+p0i*d2p1jda2da2*p2k+p0i*p1j*d2p2kda1da2
        +(dp0ida1*dp1jda2+dp0ida2*dp1jda2)*p2k
        +(dp1jda2*dp2kda1+dp1jda2*dp2kda2)*p0i
        +(dp2kda2*dp0ida2+dp2kda1*dp0ida1)*p1j;
      ggv[1][1]=
        d2p0ida2da2*p1j*p2k+p0i*d2p1jda1da1*p2k+p0i*p1j*d2p2kda2da2
        +2.*(dp0ida2*dp1jda1*p2k+p0i*dp1jda1*dp2kda2+dp0ida2*p1j*dp2kda2);
      ggv[1][2]=
        d2p0ida2da2*p1j*p2k+p0i*d2p1jda1da2*p2k+p0i*p1j*d2p2kda1da2
        +(dp0ida2*dp1jda2+dp0ida2*dp1jda1)*p2k
        +(dp1jda1*dp2kda1+dp1jda2*dp2kda2)*p0i
        +(dp2kda2*dp0ida2+dp2kda1*dp0ida2)*p1j;
      ggv[2][2]=
        d2p0ida2da2*p1j*p2k+p0i*d2p1jda2da2*p2k+p0i*p1j*d2p2kda1da1
        +2.*(dp0ida2*dp1jda2*p2k+p0i*dp1jda2*dp2kda1+dp0ida2*p1j*dp2kda1);
      ggv[1][0]=ggv[0][1];
      ggv[2][0]=ggv[0][2];
      ggv[2][1]=ggv[1][2];
    }
  }
  for (int d=0;d<3;d++) {
    prod[d].cleanup();
    dprodda1[d].cleanup();
    dprodda2[d].cleanup();
    d2prodda1da1[d].cleanup();
    d2prodda1da2[d].cleanup();
    d2prodda2da2[d].cleanup();
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void HomogeneousLagrangePolynomials3::printOn(ostream &os) const {
  os << "HomogeneousLagrangePolynomials3 : order = " << order
     << "\n\tn_pols = " << n_pols << endl;
  for (int i=0;i<zeros.getNumber();i++) {
    cout << "\tzeros[" << i << "] = " << zeros[i] << endl;
  }
  ShapeFunction<3>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void TetrahedronPolynomials::printOn(ostream &os) const {
  os << "TetrahedronPolynomials: order = " << order << endl;
  ShapeFunction<3>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void PrismPolynomials::computeValues(const Point<3> &p,NumPtr<REAL> &v) 
const {
  int ntp=tp->getNumber();
  int nlp=order2+1;
  v.allocate(ntp*nlp);
  NumPtr<REAL> tpv;
  tp->computeValues(Point<2>(p[0],p[1]),tpv);
  NumPtr<REAL> lpv(nlp);
  lp->values(p[2],lpv);
  int iv=0;
  for (int k=0;k<nlp;k++) {
    REAL &lpvk=lpv[k];
    for (int ij=0;ij<ntp;ij++) {
      v[iv++]=tpv[ij]*lpvk;
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void PrismPolynomials::computeGrads(const Point<3> &p,
NumPtr<Tensor<1,3> > &g) const {
  int ntp=tp->getNumber();
  int nlp=order2+1;
  g.allocate(ntp*nlp);
  NumPtr<REAL> tpv;
  NumPtr<Tensor<1,2> > tps;
  NumPtr<REAL> lpv(nlp);
  NumPtr<REAL> lps(nlp);
  Point<2> p2(p[0],p[1]);
  tp->computeValues(p2,tpv);
  tp->computeGrads(p2,tps);
  lp->values(p[2],lpv);
  lp->slopes(p[2],lps);
  int iv=0;
  for (int k=0;k<nlp;k++) {
    REAL lpvk=lpv[k];
    REAL lpsk=lps[k];
    for (int ij=0;ij<ntp;ij++) {
      Tensor<1,3> &gv=g[iv++];
      Tensor<1,2> &tpsij=tps[ij];
      gv[0]=tpsij[0]*lpvk;
      gv[1]=tpsij[1]*lpvk;
      gv[2]=tpv[ij]*lpsk;
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void PrismPolynomials::computeGradGrads(const Point<3> &p,
NumPtr<Tensor<2,3> > &gg) const {
  int ntp=tp->getNumber();
  int nlp=order2+1;
  gg.allocate(ntp*nlp);
  NumPtr<REAL> tpv;
  NumPtr<Tensor<1,2> > tps;
  NumPtr<Tensor<2,2> > tps2;
  NumPtr<REAL> lpv(nlp);
  NumPtr<REAL> lps(nlp);
  NumPtr<REAL> lps2(nlp);
  Point<2> p2(p[0],p[1]);
  tp->computeValues(p2,tpv);
  tp->computeGrads(p2,tps);
  tp->computeGradGrads(p2,tps2);
  lp->values(p[2],lpv);
  lp->slopes(p[2],lps);
  lp->slope2s(p[2],lps2);
  int iv=0;
  for (int k=0;k<nlp;k++) {
    REAL &lpvk=lpv[k];
    REAL &lpsk=lps[k];
    REAL &lps2k=lps2[k];
    for (int ij=0;ij<ntp;ij++) {
      Tensor<2,3> &ggv=gg[iv++];
      Tensor<1,2> &tpsij=tps[ij];
      Tensor<2,2> &tps2ij=tps2[ij];
      ggv[0][0]=tps2ij[0][0]*lpvk;
      ggv[1][0]=tps2ij[1][0]*lpvk;
      ggv[1][1]=tps2ij[1][1]*lpvk;
      ggv[2][0]=tpsij[0]*lpsk;
      ggv[2][1]=tpsij[1]*lpsk;
      ggv[2][2]=tpv[ij]*lps2k;
      ggv[0][1]=ggv[1][0];
      ggv[0][2]=ggv[2][0];
      ggv[1][2]=ggv[2][1];
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void PrismPolynomials::printOn(ostream &os) const {
  os << "PrismPolynomials: order01 = " << order01
     << "\n\torder2 = " << order2 << endl;
  os << "\tlp:" << endl;
  lp->printOn(os);
  os << "\ttp:" << endl;
  tp->printOn(os);
  ShapeFunction<3>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void C0LagrangeIntervalPolynomials::printOn(ostream &os) const {
  os << "C0LagrangeIntervalPolynomials:" << endl;
  TensorProductPolynomials<1>::printOn(os);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline void setZeros(int order,NumPtr<REAL> &zeros) {
  if (zeros.getNumber()>0) zeros[0]=1.;
  if (order>0) {
    zeros[1]=0.;
    REAL ro=1./static_cast<REAL>(order);
    for (int i=1;i<order;i++) zeros[i+1]=static_cast<REAL>(i)*ro;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
C0LagrangeTrianglePolynomials::C0LagrangeTrianglePolynomials(int ord) :
TrianglePolynomials(ord),zeros(ord+1) {
  setZeros(order,zeros);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
void numden(const NumPtr<REAL> &zeros,int n,REAL bc,REAL sup,
REAL &num) {
  num=1.;
  REAL den=1.;
  for (int m=1;m<=n;m++) {
    const REAL &z=zeros[m];
    num*=bc-z;
    den*=sup-z;
  }
  num /= den;
}
void numdens(const NumPtr<REAL> &zeros,int n,REAL bc,REAL sup,
REAL &num,REAL &dnum) {
  num=1.;
  dnum=0.;
  REAL den=1.;
  for (int m=1;m<=n;m++) {
    const REAL &z=zeros[m];
    num*=bc-z;
    den*=sup-z;
    REAL prod=1.;
    for (int l=1;l<=n;l++) if (l!=m) prod*=bc-zeros[l];
    dnum+=prod;
  }
  num /= den;
  dnum /= den;
}
void numdens2(const NumPtr<REAL> &zeros,int n,REAL bc,REAL sup,
REAL &num,REAL &dnum,REAL &d2num) {
  num=1.;
  dnum=0.;
  d2num=0.;
  REAL den=1.;
  for (int m=1;m<=n;m++) {
    const REAL &z=zeros[m];
    num*=bc-z;
    den*=sup-z;
    REAL prod=1.;
    for (int l=1;l<=n;l++) {
      if (l==m) continue;
      prod*=bc-zeros[l];
      REAL prod2=1.;
      for (int k=1;k<=n;k++) {
        if (k!=m && k!=l) prod2*=bc-zeros[k];
      }
      d2num+=prod2;
    }
    dnum+=prod;
  }
  num /= den;
  dnum /= den;
  d2num /= den;
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
void C0LagrangeTrianglePolynomials::computeValues(const Point<2> &p,
NumPtr<REAL> &v) const {
  int n_polys=getNumber();
  v.allocate(n_polys);
  Point<3> bc(p[0],p[1],1.-(p[0]+p[1]));
  NumPtr<NumPtr<REAL> > pnum(3);
  for (int d=0;d<3;d++) {
    NumPtr<REAL> &pnumd=pnum[d];
    pnumd.allocate(order+1);
    for (int i=0;i<=order;i++) {
      REAL sup=zeros[(i==order ? 0 : i+1)];
      numden(zeros,i,bc[d],sup,pnumd[i]);
    }
  }

  int iv=0;
  NumPtr<REAL> &pnum0=pnum[0];
  NumPtr<REAL> &pnum1=pnum[1];
  NumPtr<REAL> &pnum2=pnum[2];
  for (int j=0;j<=order;j++) {
    REAL &num1=pnum1[j];
    for (int i=0;i+j<=order;i++) {
      int k=order-i-j;
      v[iv++]=num1*pnum0[i]*pnum2[k];
    }
  }

  for (int d=0;d<3;d++) {
    pnum[d].cleanup();
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
void C0LagrangeTrianglePolynomials::computeGrads(const Point<2> &p,
NumPtr<Tensor<1,2> > &g) const {
  int n_polys=getNumber();
  g.allocate(n_polys);
  Point<3> bc(p[0],p[1],1.-(p[0]+p[1]));
  NumPtr<NumPtr<REAL> > pnum(3);
  NumPtr<NumPtr<REAL> > pdnum(3);
  for (int d=0;d<3;d++) {
    NumPtr<REAL> &pnumd=pnum[d];
    NumPtr<REAL> &pdnumd=pdnum[d];
    pnumd.allocate(order+1);
    pdnumd.allocate(order+1);
    for (int i=0;i<=order;i++) {
      REAL sup=zeros[(i==order ? 0 : i+1)];
      numdens(zeros,i,bc[d],sup,pnumd[i],pdnumd[i]);
    }
  }

  NumPtr<REAL> &pnum0=pnum[0];
  NumPtr<REAL> &pnum1=pnum[1];
  NumPtr<REAL> &pnum2=pnum[2];
  NumPtr<REAL> &pdnum0=pdnum[0];
  NumPtr<REAL> &pdnum1=pdnum[1];
  NumPtr<REAL> &pdnum2=pdnum[2];
  int iv=0;
  for (int j=0;j<=order;j++) {
    REAL &num1=pnum1[j];
    REAL &dnum1=pdnum1[j];
    for (int i=0;i+j<=order;i++) {
      int k=order-i-j;
      REAL &num0=pnum0[i];
      REAL &num2=pnum2[k];
      REAL &dnum2=pdnum2[k];
      Tensor<1,2> &gv=g[iv++];
      gv[0]=num1*(pdnum0[i]*num2-num0*dnum2);
      gv[1]=num0*(dnum1*num2-num1*dnum2);
    }
  }

  for (int d=0;d<3;d++) {
    pnum[d].cleanup();
    pdnum[d].cleanup();
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
void C0LagrangeTrianglePolynomials::computeGradGrads(const Point<2> &p,
NumPtr<Tensor<2,2> > &gg) const {
  int n_polys=getNumber();
  gg.allocate(n_polys);
  Point<3> bc(p[0],p[1],1.-(p[0]+p[1]));
  NumPtr<NumPtr<REAL> > pnum(3);
  NumPtr<NumPtr<REAL> > pdnum(3);
  NumPtr<NumPtr<REAL> > pd2num(3);
  for (int d=0;d<3;d++) {
    NumPtr<REAL> &pnumd=pnum[d];
    NumPtr<REAL> &pdnumd=pdnum[d];
    NumPtr<REAL> &pd2numd=pd2num[d];
    pnumd.allocate(order+1);
    pdnumd.allocate(order+1);
    pd2numd.allocate(order+1);
    for (int i=0;i<=order;i++) {
      REAL sup=zeros[(i==order ? 0 : i+1)];
      numdens2(zeros,i,bc[d],sup,pnumd[i],pdnumd[i],pd2numd[i]);
    }
  }

  NumPtr<REAL> &pnum0=pnum[0];
  NumPtr<REAL> &pnum1=pnum[1];
  NumPtr<REAL> &pnum2=pnum[2];
  NumPtr<REAL> &pdnum0=pdnum[0];
  NumPtr<REAL> &pdnum1=pdnum[1];
  NumPtr<REAL> &pdnum2=pdnum[2];
  NumPtr<REAL> &pd2num0=pd2num[0];
  NumPtr<REAL> &pd2num1=pd2num[1];
  NumPtr<REAL> &pd2num2=pd2num[2];
  int iv=0;
  for (int j=0;j<=order;j++) {
    REAL &num1=pnum1[j];
    REAL &dnum1=pdnum1[j];
    REAL &d2num1=pd2num1[j];
    for (int i=0;i+j<=order;i++) {
      int k=order-i-j;
      REAL &num0=pnum0[i];
      REAL &dnum0=pdnum0[i];
      REAL &num2=pnum2[k];
      REAL &dnum2=pdnum2[k];
      REAL &d2num2=pd2num2[k];
      Tensor<2,2> &ggv=gg[iv++];
      ggv[0][0]=num1*(pd2num0[i]*num2-2.*dnum0*dnum2+num0*d2num2);
      ggv[1][0]=dnum0*dnum1*num2+num0*num1*d2num2
               -(num0*dnum1+dnum0*num1)*dnum2;
      ggv[1][1]=num0*(d2num1*num2-2.*dnum1*dnum2+num1*d2num2);
      ggv[0][1]=ggv[1][0];
    }
  }

  for (int d=0;d<3;d++) {
    pnum[d].cleanup();
    pdnum[d].cleanup();
    pd2num[d].cleanup();
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
void C0LagrangeTrianglePolynomials::shapeFunctionsOrder(
const Shape *element,NumPtr<int> &renumber) const {
  CHECK_POINTER(dynamic_cast<const Triangle*>(element));
  element->reorderEquallySpacedLatticePoints(order,renumber);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
void C0LagrangeTrianglePolynomials::printOn(ostream &os) const {
  os << "C0LagrangeTrianglePolynomials: order = " << order << endl;
  os << "\tzeros =";
  for (int i=0;i<zeros.getNumber();i++) os << " " << zeros[i];
  os << endl;
  ShapeFunction<2>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
void C0LagrangeQuadrilateralPolynomials::printOn(ostream &os) const {
  os << "C0LagrangeQuadrilateralPolynomials" << endl;
  TensorProductPolynomials<2>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
C0LagrangeTetrahedronPolynomials::C0LagrangeTetrahedronPolynomials(
int ord) : TetrahedronPolynomials(ord),zeros(ord+1) {
  setZeros(order,zeros);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void C0LagrangeTetrahedronPolynomials::computeValues(const Point<3> &p,
NumPtr<REAL> &v) const {
  int n_polys=getNumber();
  v.allocate(n_polys);
  Point<4> bc(p[0],p[1],p[2],1.-(p[0]+p[1]+p[2]));
  NumPtr<NumPtr<REAL> > pnum(4);
  for (int d=0;d<4;d++) {
    NumPtr<REAL> &pnumd=pnum[d];
    pnumd.allocate(order+1);
    for (int i=0;i<=order;i++) {
      REAL sup=zeros[(i==order ? 0 : i+1)];
      numden(zeros,i,bc[d],sup,pnumd[i]);
    }
  }

  int iv=0;
  NumPtr<REAL> &pnum0=pnum[0];
  NumPtr<REAL> &pnum1=pnum[1];
  NumPtr<REAL> &pnum2=pnum[2];
  NumPtr<REAL> &pnum3=pnum[3];
  for (int k=0;k<=order;k++) {
    REAL &num2=pnum2[k];
    for (int j=0;j+k<=order;j++) {
      REAL &num1=pnum1[j];
      for (int i=0;i+j+k<=order;i++) {
        int l=order-i-j-k;
        v[iv++]=num1*num2*pnum0[i]*pnum3[l];
      }
    }
  }

  for (int d=0;d<4;d++) {
    pnum[d].cleanup();
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void C0LagrangeTetrahedronPolynomials::computeGrads(const Point<3> &p,
NumPtr<Tensor<1,3> > &g) const {
  int n_polys=getNumber();
  g.allocate(n_polys);
  Point<4> bc(p[0],p[1],p[2],1.-(p[0]+p[1]+p[2]));
  NumPtr<NumPtr<REAL> > pnum(4);
  NumPtr<NumPtr<REAL> > pdnum(4);
  for (int d=0;d<4;d++) {
    NumPtr<REAL> &pnumd=pnum[d];
    NumPtr<REAL> &pdnumd=pdnum[d];
    pnumd.allocate(order+1);
    pdnumd.allocate(order+1);
    for (int i=0;i<=order;i++) {
      REAL sup=zeros[(i==order ? 0 : i+1)];
      numdens(zeros,i,bc[d],sup,pnumd[i],pdnumd[i]);
    }
  }

  NumPtr<REAL> &pnum0=pnum[0];
  NumPtr<REAL> &pnum1=pnum[1];
  NumPtr<REAL> &pnum2=pnum[2];
  NumPtr<REAL> &pnum3=pnum[3];
  NumPtr<REAL> &pdnum0=pdnum[0];
  NumPtr<REAL> &pdnum1=pdnum[1];
  NumPtr<REAL> &pdnum2=pdnum[2];
  NumPtr<REAL> &pdnum3=pdnum[3];
  int iv=0;
  for (int k=0;k<=order;k++) {
    REAL &num2=pnum2[k];
    REAL &dnum2=pdnum2[k];
    for (int j=0;j+k<=order;j++) {
      REAL &num1=pnum1[j];
      REAL &dnum1=pdnum1[j];
      for (int i=0;i+j+k<=order;i++) {
        int l=order-i-j-k;
        REAL &num0=pnum0[i];
        REAL &num3=pnum3[l];
        REAL &dnum3=pdnum3[l];
        Tensor<1,3> &gv=g[iv++];
        gv[0]=num1*num2*(pdnum0[i]*num3-num0*dnum3);
        gv[1]=num0*num2*(dnum1*num3-num1*dnum3);
        gv[2]=num0*num1*(dnum2*num3-num2*dnum3);
      }
    }
  }

  for (int d=0;d<4;d++) {
    pnum[d].cleanup();
    pdnum[d].cleanup();
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void C0LagrangeTetrahedronPolynomials::computeGradGrads(
const Point<3> &p,NumPtr<Tensor<2,3> > &gg) const {
  int n_polys=getNumber();
  gg.allocate(n_polys);
  Point<4> bc(p[0],p[1],p[2],1.-(p[0]+p[1]+p[2]));
  NumPtr<NumPtr<REAL> > pnum(4);
  NumPtr<NumPtr<REAL> > pdnum(4);
  NumPtr<NumPtr<REAL> > pd2num(4);
  for (int d=0;d<4;d++) {
    NumPtr<REAL> &pnumd=pnum[d];
    NumPtr<REAL> &pdnumd=pdnum[d];
    NumPtr<REAL> &pd2numd=pd2num[d];
    pnumd.allocate(order+1);
    pdnumd.allocate(order+1);
    pd2numd.allocate(order+1);
    for (int i=0;i<=order;i++) {
      REAL sup=zeros[(i==order ? 0 : i+1)];
      numdens2(zeros,i,bc[d],sup,pnumd[i],pdnumd[i],pd2numd[i]);
    }
  }

  NumPtr<REAL> &pnum0=pnum[0];
  NumPtr<REAL> &pnum1=pnum[1];
  NumPtr<REAL> &pnum2=pnum[2];
  NumPtr<REAL> &pnum3=pnum[3];
  NumPtr<REAL> &pdnum0=pdnum[0];
  NumPtr<REAL> &pdnum1=pdnum[1];
  NumPtr<REAL> &pdnum2=pdnum[2];
  NumPtr<REAL> &pdnum3=pdnum[3];
  NumPtr<REAL> &pd2num0=pd2num[0];
  NumPtr<REAL> &pd2num1=pd2num[1];
  NumPtr<REAL> &pd2num2=pd2num[2];
  NumPtr<REAL> &pd2num3=pd2num[3];
  int iv=0;
  for (int k=0;k<=order;k++) {
    REAL &num2=pnum2[k];
    REAL &dnum2=pdnum2[k];
    REAL &d2num2=pd2num2[k];
    for (int j=0;j+k<=order;j++) {
      REAL &num1=pnum1[j];
      REAL &dnum1=pdnum1[j];
      REAL &d2num1=pd2num1[j];
      for (int i=0;i+j+k<=order;i++) {
        int l=order-i-j-k;
        REAL &num0=pnum0[i];
        REAL &dnum0=pdnum0[i];
        REAL &num3=pnum3[l];
        REAL &dnum3=pdnum3[l];
        REAL &d2num3=pd2num3[l];
        Tensor<2,3> &ggv=gg[iv++];
        ggv[0][0]=(pd2num0[i]*num3-2.*dnum0*dnum3+num0*d2num3)
                 *num1*num2;
        ggv[1][0]=((-dnum1*dnum3+num1*d2num3)*num0
                  +(dnum1*num3-num1*dnum3)*dnum0)*num2;
        ggv[2][0]=((-dnum2*dnum3+num2*d2num3)*num0
                  +(dnum2*num3-num2*dnum3)*dnum0)*num1;
        ggv[1][1]=(d2num1*num3-2.*dnum1*dnum3+num1*d2num3)
                 *num2*num0;
        ggv[2][1]=((-dnum2*dnum3+num2*d2num3)*num1
                  +(dnum2*num3-num2*dnum3)*dnum1)*num0;
        ggv[2][2]=(d2num2*num3-2.*dnum2*dnum3+num2*d2num3)
                 *num0*num1;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    }
  }

  for (int d=0;d<4;d++) {
    pnum[d].cleanup();
    pdnum[d].cleanup();
    pd2num[d].cleanup();
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void C0LagrangeTetrahedronPolynomials::shapeFunctionsOrder(
const Shape *element,NumPtr<int> &renumber) const {
  CHECK_POINTER(dynamic_cast<const Tetrahedron*>(element));
  element->reorderEquallySpacedLatticePoints(order,renumber);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void C0LagrangeTetrahedronPolynomials::printOn(ostream &os) const {
  os << "C0LagrangeTetrahedronPolynomials: order = " << order << endl;
  ShapeFunction<3>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void C0LagrangePrismPolynomials::shapeFunctionsOrder(
const Shape *element,NumPtr<int> &renumber) const {
  CHECK_POINTER(dynamic_cast<const Prism*>(element));
  CHECK_SAME(order01,order2);
  element->reorderEquallySpacedLatticePoints(order01,renumber);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void C0LagrangePrismPolynomials::printOn(ostream &os) const {
  os << "C0LagrangePrismPolynomials:" << endl;
  PrismPolynomials::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void C0LagrangeHexahedronPolynomials::printOn(ostream &os) const {
  os << "C0LagrangeHexahedronPolynomials" << endl;
  TensorProductPolynomials<3>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void HierarchicalIntervalPolynomials::shapeFunctionsOrder(
const Shape *element,NumPtr<int> &renumber) const {
  CHECK_POINTER(dynamic_cast<const Interval*>(element));
  renumber.allocate(getNumber());
  for (int i=0;i<renumber.getNumber();i++) { renumber[i]=i; }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void HierarchicalIntervalPolynomials::computeValues(const Point<1> &p,
NumPtr<REAL> &v) const {
  v.allocate(n_pols);
  HierarchicalPolynomial hp;
  hp.values(p[0],v);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void HierarchicalIntervalPolynomials::computeGrads(const Point<1> &p,
NumPtr<Tensor<1,1> > &g) const {
  g.allocate(n_pols);
  NumPtr<REAL> s(n_pols);
  HierarchicalPolynomial hp;
  hp.slopes(p[0],s);
  for (int i=0;i<n_pols;i++) g[i][0]=s[i];
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void HierarchicalIntervalPolynomials::computeGradGrads(
const Point<1> &p,NumPtr<Tensor<2,1> > &gg) const {
  gg.allocate(n_pols);
  NumPtr<REAL> s2(n_pols);
  HierarchicalPolynomial hp;
  hp.slope2s(p[0],s2);
  for (int i=0;i<n_pols;i++) gg[i][0][0]=s2[i];
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void HierarchicalIntervalPolynomials::printOn(ostream &os) const {
  os << "HierarchicalIntervalPolynomials: n_pols = " << n_pols << endl;
  ShapeFunction<1>::printOn(os);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
/*
following Peano, Comp. & Maths. with Appls 2(1976) 211-224
the following polynomials in the barycentric coordinates (a,b,c) form 
a basis for polynomials of 2 variables of order n:
  1
  if n>0: a, b
  for 1<=k < n: a^k*b, b^k*c, c^k*a
  for 3<=3m <= n: (a*b*c)^m
  for 4<=3m+1<=n: a*(a*b*c)^m, b*(a*b*c)^m
  for 5<=3m+2<=n, for 1 <= k < n-3m: 
    a^k*b*(a*b*c)^m, b^k*c*(a*b*c)^m, c^k*a*(a*b*c)^m
The first 3 entries (1,a,b) are associated with the triangle vertices.
The polynomials a^k*b, b^k*c, c^k*a are associated with the triangle
sides.
The remaining polynomials are associated with the triangle interior.

We have to guarantee that the hierarchical polynomials agree with those
on quadrilaterals along shared lines.
Along side c=0, a+b=1
  B_0 = 1-b = a
  B_1 = b
  B_2 = sqrt(6) b (b-1) = a b {- sqrt(6)} = a b phi_0(b-a)
  B_3 = sqrt(40) b (b-1) (b-1/2) = a b {- sqrt(10) (b-a)} = a b phi_1(b-a)
Here
  phi_k(x) = h_{k+2}(x) / [ x ( 1 - x ) ]
where  h_{k+2}(x) is a hierarchical polynomial.
The forms using phi_k are guaranteed to vanish on the other sides of the
triangle.

The basis functions are partially described in Szabo and Babuska, p. 103
and better described for tetrahedra on p. 244
Assuming that n>0, this leads to the following basis:
  vertex modes ordered 
    c, (v0) 
    a, (v1)
    b  (v2)
  line modes ordered
    for 1 <= k < n: 
      b*c*phi_k(c-b), (v0,v2)
      c*a*phi_k(a-c), (v0,v1)
      a*b*phi_k(b-a), (v1,v2)
  interior modes (order unimportant):
    for 0 <= m <= n-3
      for 0 <= j <= m
        a*b*c * P_{m-j}(b-a) * P_j(2*c-1)
  where P_i(x) is the Legendre polynomial of degree i.
*/
void HierarchicalTrianglePolynomials::computeValues(const Point<2> &p,
NumPtr<REAL> &v) const {
  int n_polys=getNumber();
  v.allocate(n_polys);
  Point<3> bc(p[0],p[1],1.-(p[0]+p[1]));
  int ij=0;
  v[ij++]=bc[2];
  v[ij++]=bc[0];
  v[ij++]=bc[1];
  if (order>=2) {
    HierarchicalSidePolynomial hsp;
    NumPtr<REAL> hspv(order-1);
//  line 0: v0 -> v2
    hsp.values(bc[1],hspv);
    for (int j=1;j<order;j++) {
      v[ij++]=bc[1]*bc[2]*hspv[j-1];
    }
//  line 1: v0 -> v1
    hsp.values(bc[0],hspv);
    for (int i=1;i<order;i++) {
      v[ij++]=bc[2]*bc[0]*hspv[i-1];
    }
//  line 2: v1 -> v2
    hsp.values(bc[1],hspv);
    for (int k=1;k<order;k++) {
      v[ij++]=bc[0]*bc[1]*hspv[k-1];
    }
    if (order>=3) {
      REAL abc=bc[0]*bc[1]*bc[2];
      LegendrePolynomial lp;
      NumPtr<REAL> lpv0(order-2);
      NumPtr<REAL> lpv1(order-2);
      lp.values(bc[0],lpv0);
      lp.values(bc[1],lpv1);
      for (int m=0;m+3<=order;m++) {
        for (int j=0;j<=m;j++) {
          v[ij++]=abc*lpv0[m-j]*lpv1[j];
        }
      }
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
void HierarchicalTrianglePolynomials::computeGrads(const Point<2> &p,
NumPtr<Tensor<1,2> > &g) const {
  int n_polys=getNumber();
  g.allocate(n_polys);
  Point<3> bc(p[0],p[1],1.-(p[0]+p[1]));

  int ij=0;
  Tensor<1,2> &gv0=g[ij++];
  gv0[0]=-1.;
  gv0[1]=-1.;
  Tensor<1,2> &gv1=g[ij++];
  gv1[0]=1.;
  gv1[1]=0.;
  Tensor<1,2> &gv2=g[ij++];
  gv2[0]=0.;
  gv2[1]=1.;
  if (order>=2) {
    HierarchicalSidePolynomial hsp;
    NumPtr<REAL> hspv(order-1);
    NumPtr<REAL> hsps(order-1);
//  line0: v0 -> v2
    hsp.values(bc[1],hspv);
    hsp.slopes(bc[1],hsps);
    for (int j=1;j<order;j++) {
      Tensor<1,2> &gij=g[ij++];
      gij[0]=-bc[1]*hspv[j-1];
      gij[1]=(bc[2]-bc[1])*hspv[j-1]+bc[1]*bc[2]*hsps[j-1];
    }
//  line 1: v0 -> v1
    hsp.values(bc[0],hspv);
    hsp.slopes(bc[0],hsps);
    for (int i=1;i<order;i++) {
      Tensor<1,2> &gij=g[ij++];
      gij[0]=(bc[2]-bc[0])*hspv[i-1]+bc[2]*bc[0]*hsps[i-1];
      gij[1]=-bc[0]*hspv[i-1];
    }
//  line 2: v1 -> v2
    hsp.values(bc[1],hspv);
    hsp.slopes(bc[1],hsps);
    for (int k=1;k<order;k++) {
      Tensor<1,2> &gij=g[ij++];
      gij[0]=bc[1]*hspv[k-1];
      gij[1]=bc[0]*(hspv[k-1]+bc[1]*hsps[k-1]);
    }
    if (order>=3) {
      REAL abc=bc[0]*bc[1]*bc[2];
      Tensor<1,2> gabc;
      gabc[0]=bc[1]*(bc[2]-bc[0]);
      gabc[1]=bc[0]*(bc[2]-bc[1]);
      LegendrePolynomial lp;
      NumPtr<REAL> lpv0(order-2);
      NumPtr<REAL> lps0(order-2);
      NumPtr<REAL> lpv1(order-2);
      NumPtr<REAL> lps1(order-2);
      lp.values(bc[0],lpv0);
      lp.slopes(bc[0],lps0);
      lp.values(bc[1],lpv1);
      lp.slopes(bc[1],lps1);
      for (int m=0;m+3<=order;m++) {
        for (int j=0;j<=m;j++) {
          Tensor<1,2> &gij=g[ij++];
          gij[0]=(gabc[0]*lpv0[m-j]+abc*lps0[m-j])*lpv1[j];
          gij[1]=(gabc[1]*lpv1[j]+abc*lps1[j])*lpv0[m-j];
        }
      }
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
void HierarchicalTrianglePolynomials::computeGradGrads(const Point<2> &p,
NumPtr<Tensor<2,2> > &gg) const {
  int n_polys=getNumber();
  gg.allocate(n_polys);
  Point<3> bc(p[0],p[1],1.-(p[0]+p[1]));

  int ij=0;
  Tensor<2,2> &ggv0=gg[ij++];
  ggv0[0][0]=0.;
  ggv0[1][0]=0.;
  ggv0[0][1]=0.;
  ggv0[1][1]=0.;
  Tensor<2,2> &ggv1=gg[ij++];
  ggv1[0][0]=0.;
  ggv1[1][0]=0.;
  ggv1[0][1]=0.;
  ggv1[1][1]=0.;
  Tensor<2,2> &ggv2=gg[ij++];
  ggv2[0][0]=0.;
  ggv2[1][0]=0.;
  ggv2[0][1]=0.;
  ggv2[1][1]=0.;
  if (order>=2) {
    HierarchicalSidePolynomial hsp;
    NumPtr<REAL> hspv(order-1);
    NumPtr<REAL> hsps(order-1);
    NumPtr<REAL> hspss(order-1);
//  line0: v0 -> v2
    hsp.values(bc[1],hspv);
    hsp.slopes(bc[1],hsps);
    hsp.slope2s(bc[1],hspss);
    for (int j=1;j<order;j++) {
      Tensor<2,2> &ggij=gg[ij++];
      ggij[0][0]=0.;
      ggij[1][0]=-hspv[j-1]-bc[1]*hsps[j-1];
      ggij[1][1]=-2.*hspv[j-1]+2.*(bc[2]-bc[1])*hsps[j-1]
                +bc[1]*bc[2]*hspss[j-1];
      ggij[0][1]=ggij[1][0];
    }
//  line 1: v0 -> v1
    hsp.values(bc[0],hspv);
    hsp.slopes(bc[0],hsps);
    hsp.slope2s(bc[0],hspss);
    for (int i=1;i<order;i++) {
      Tensor<2,2> &ggij=gg[ij++];
      ggij[0][0]=-2.*hspv[i-1]+2.*(bc[2]-bc[0])*hsps[i-1]
                +bc[2]*bc[0]*hspss[i-1];
      ggij[1][0]=-hspv[i-1]-bc[0]*hsps[i-1];
      ggij[1][1]=0.;
      ggij[0][1]=ggij[1][0];
    }
//  line 2: v1 -> v2
    hsp.values(bc[1],hspv);
    hsp.slopes(bc[1],hsps);
    hsp.slope2s(bc[1],hspss);
    for (int k=1;k<order;k++) {
      Tensor<2,2> &ggij=gg[ij++];
      ggij[0][0]=0.;
      ggij[1][0]=hspv[k-1]+bc[1]*hsps[k-1];
      ggij[1][1]=2.*bc[0]*hsps[k-1]+bc[0]*bc[1]*hspss[k-1];
      ggij[0][1]=ggij[1][0];
    }
    if (order>=3) {
      REAL abc=bc[0]*bc[1]*bc[2];
      Tensor<1,2> gabc;
      gabc[0]=bc[1]*(bc[2]-bc[0]);
      gabc[1]=bc[0]*(bc[2]-bc[1]);
      Tensor<2,2> ggabc;
      ggabc[0][0]=-2.*bc[1];
      ggabc[1][0]=bc[2]-bc[0]-bc[1];
      ggabc[1][1]=-2.*bc[0];
      ggabc[0][1]=ggabc[1][0];
      LegendrePolynomial lp;
      NumPtr<REAL> lpv0(order-2);
      NumPtr<REAL> lps0(order-2);
      NumPtr<REAL> lpss0(order-2);
      NumPtr<REAL> lpv1(order-2);
      NumPtr<REAL> lps1(order-2);
      NumPtr<REAL> lpss1(order-2);
      lp.values(bc[0],lpv0);
      lp.slopes(bc[0],lps0);
      lp.slope2s(bc[0],lpss0);
      lp.values(bc[1],lpv1);
      lp.slopes(bc[1],lps1);
      lp.slope2s(bc[1],lpss1);
      for (int m=0;m+3<=order;m++) {
        for (int j=0;j<=m;j++) {
          Tensor<2,2> &ggij=gg[ij++];
          ggij[0][0]=(ggabc[0][0]*lpv0[m-j]+2.*gabc[0]*lps0[m-j]
                     +abc*lpss0[m-j])*lpv1[j];
          ggij[1][0]=ggabc[1][0]*lpv0[m-j]*lpv1[j]
            +gabc[1]*lps0[m-j]*lpv1[j]+gabc[0]*lpv0[m-j]*lps1[j]
            +abc*lps0[m-j]*lps1[j];
          ggij[1][1]=(ggabc[1][1]*lpv1[j]+2.*gabc[1]*lps1[j]
                     +abc*lpss1[j])*lpv0[m-j];
          ggij[0][1]=ggij[1][0];
        }
      }
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
void HierarchicalTrianglePolynomials::shapeFunctionsOrder(
const Shape *element,NumPtr<int> &renumber) const {
  CHECK_POINTER(dynamic_cast<const Triangle*>(element));
  renumber.allocate(getNumber());
  for (int i=0;i<renumber.getNumber();i++) { renumber[i]=i; }
  if (order<2) return;
  int nl=order-1;
//line 0: v0 -> v2
  int a=2;
  if (element->vertexIndex(2)<element->vertexIndex(0)) {
    int s=1;
    for (int j=1;j<order;j++) {
      renumber[a+j]=s*(a+j);
      s=-s;
    }
  }
//line 1: v0 -> v1
  a+=nl;
  if (element->vertexIndex(1)<element->vertexIndex(0)) {
    int s=1;
    for (int i=1;i<order;i++) {
      renumber[a+i]=s*(a+i);
      s=-s;
    }
  }
//line 2: v1 -> v2
  a+=nl;
  if (element->vertexIndex(2)<element->vertexIndex(1)) {
    int s=1;
    for (int j=1;j<order;j++) {
      renumber[a+j]=s*(a+j);
      s=-s;
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
void HierarchicalTrianglePolynomials::nonzeroShapesOnFace(
const Shape1 *e,int f,NumPtr<int> &indices) const {
  indices.allocate(order+1);
  int iv=0;
  for (int v=0;v<2;v++) indices[iv++]=e->elementVertexFromFace(f,v);
  int offset=e->numberVertices()+f*(order-1)-1;
  for (int i=1;i<order;i++) indices[iv++]=offset+i;
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
void HierarchicalTrianglePolynomials::printOn(ostream &os) const {
  os << "HierarchicalTrianglePolynomials: order = " << order << endl;
  ShapeFunction<2>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
void HierarchicalQuadrilateralPolynomials::computeValues(
const Point<2> &p,NumPtr<REAL> &v) const {
  v.allocate(n_pols);
  NumPtr<NumPtr<REAL> > hpv(2);
  HierarchicalPolynomial hp;
  for (int d=0;d<SPACEDIM;d++) {
    hpv[d].allocate(n_pols);
    hp.values(p[d],hpv[d]);
  }
  int iv=0;
  NumPtr<REAL> &hpv0=hpv[0];
  NumPtr<REAL> &hpv1=hpv[1];
//vertices
  for (int j=0;j<2;j++) {
    const REAL &hpvj=hpv1[j];
    for (int i=0;i<2;i++) {
      v[iv++]=hpv0[i]*hpvj;
    }
  }
//side 0
  REAL hpvi=hpv0[0];
  for (int j=2;j<n_pols_1d;j++) {
    v[iv++]=hpvi*hpv1[j];
  }
//side 1
  hpvi=hpv0[1];
  for (int j=2;j<n_pols_1d;j++) {
    v[iv++]=hpvi*hpv1[j];
  }
//side 2
  REAL hpvj=hpv1[0];
  for (int i=2;i<n_pols_1d;i++) {
    v[iv++]=hpv0[i]*hpvj;
  }
//side 3
  hpvj=hpv1[1];
  for (int i=2;i<n_pols_1d;i++) {
    v[iv++]=hpv0[i]*hpvj;
  }
//interior
  for (int j=2;j<n_pols_1d;j++) {
    const REAL &hpvj=hpv1[j];
    for (int i=2;i<n_pols_1d;i++) {
      v[iv++]=hpv0[i]*hpvj;
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
void HierarchicalQuadrilateralPolynomials::computeGrads(
const Point<2> &p,NumPtr<Tensor<1,2> > &g) const {
  g.allocate(n_pols);
  NumPtr<NumPtr<REAL> > hpv(2);
  NumPtr<NumPtr<REAL> > hps(2);
  HierarchicalPolynomial hp;
  for (int d=0;d<2;d++) {
    hpv[d].allocate(n_pols_1d);
    hps[d].allocate(n_pols_1d);
    hp.values(p[d],hpv[d]);
    hp.slopes(p[d],hps[d]);
  }
  int iv=0;
  NumPtr<REAL> &hpv0=hpv[0];
  NumPtr<REAL> &hpv1=hpv[1];
  NumPtr<REAL> &hps0=hps[0];
  NumPtr<REAL> &hps1=hps[1];
//vertices
  for (int j=0;j<2;j++) {
    const REAL &hpvj=hpv1[j];
    const REAL &hpsj=hps1[j];
    for (int i=0;i<2;i++) {
      Tensor<1,2> &gv=g[iv++];
      gv[0]=hps0[i]*hpvj;
      gv[1]=hpv0[i]*hpsj;
    }
  }
//side 0
  REAL hpvi=hpv0[0];
  REAL hpsi=hps0[0];
  for (int j=2;j<n_pols_1d;j++) {
    Tensor<1,2> &gv=g[iv++];
    gv[0]=hpsi*hpv1[j];
    gv[1]=hpvi*hps1[j];
  }
//side 1
  hpvi=hpv0[1];
  hpsi=hps0[1];
  for (int j=2;j<n_pols_1d;j++) {
    Tensor<1,2> &gv=g[iv++];
    gv[0]=hpsi*hpv1[j];
    gv[1]=hpvi*hps1[j];
  }
//side 0
  REAL hpvj=hpv1[0];
  REAL hpsj=hps1[0];
  for (int i=2;i<n_pols_1d;i++) {
    Tensor<1,2> &gv=g[iv++];
    gv[0]=hps0[i]*hpvj;
    gv[1]=hpv0[i]*hpsj;
  }
//side 1
  hpvj=hpv1[1];
  hpsj=hps1[1];
  for (int i=2;i<n_pols_1d;i++) {
    Tensor<1,2> &gv=g[iv++];
    gv[0]=hps0[i]*hpvj;
    gv[1]=hpv0[i]*hpsj;
  }
//interior
  for (int j=2;j<n_pols_1d;j++) {
    const REAL &hpvj=hpv1[j];
    const REAL &hpsj=hps1[j];
    for (int i=2;i<n_pols_1d;i++) {
      Tensor<1,2> &gv=g[iv++];
      gv[0]=hps0[i]*hpvj;
      gv[1]=hpv0[i]*hpsj;
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
void HierarchicalQuadrilateralPolynomials::computeGradGrads(
const Point<2> &p,NumPtr<Tensor<2,2> > &gg) const {
  gg.allocate(n_pols);
  NumPtr<NumPtr<REAL> > hpv(2);
  NumPtr<NumPtr<REAL> > hps(2);
  NumPtr<NumPtr<REAL> > hps2(2);
  HierarchicalPolynomial hp;
  for (int d=0;d<2;d++) {
    hpv[d].allocate(n_pols);
    hps[d].allocate(n_pols);
    hps2[d].allocate(n_pols);
    hp.values(p[d],hpv[d]);
    hp.slopes(p[d],hps[d]);
    hp.slope2s(p[d],hps2[d]);
  }
  int iv=0;
  NumPtr<REAL> &hpv0=hpv[0];
  NumPtr<REAL> &hpv1=hpv[1];
  NumPtr<REAL> &hps0=hps[0];
  NumPtr<REAL> &hps1=hps[1];
  NumPtr<REAL> &hpss0=hps2[0];
  NumPtr<REAL> &hpss1=hps2[1];
//vertices
  for (int j=0;j<2;j++) {
    const REAL &hpvj=hpv1[j];
    const REAL &hpsj=hps1[j];
    const REAL &hpssj=hpss1[j];
    for (int i=0;i<2;i++) {
      Tensor<2,2> &ggv=gg[iv++];
      ggv[0][0]=hpss0[i]*hpvj;
      ggv[1][0]=hps0[i]*hpsj;
      ggv[1][1]=hpv0[i]*hpssj;
      ggv[0][1]=ggv[1][0];
    }
  }
//side 0
  REAL hpvi=hpv0[0];
  REAL hpsi=hps0[0];
  REAL hpssi=hpss0[0];
  for (int j=2;j<n_pols_1d;j++) {
    Tensor<2,2> &ggv=gg[iv++];
    ggv[0][0]=hpssi*hpv1[j];
    ggv[1][0]=hpsi*hps1[j];
    ggv[1][1]=hpvi*hpss1[j];
    ggv[0][1]=ggv[1][0];
  }
//side 1
  hpvi=hpv0[1];
  hpsi=hps0[1];
  hpssi=hpss0[1];
  for (int j=2;j<n_pols_1d;j++) {
    Tensor<2,2> &ggv=gg[iv++];
    ggv[0][0]=hpssi*hpv1[j];
    ggv[1][0]=hpsi*hps1[j];
    ggv[1][1]=hpvi*hpss1[j];
    ggv[0][1]=ggv[1][0];
  }
//side 2
  REAL hpvj=hpv1[0];
  REAL hpsj=hps1[0];
  REAL hpssj=hpss1[0];
  for (int i=2;i<n_pols_1d;i++) {
    Tensor<2,2> &ggv=gg[iv++];
    ggv[0][0]=hpss0[i]*hpvj;
    ggv[1][0]=hps0[i]*hpsj;
    ggv[1][1]=hpv0[i]*hpssj;
    ggv[0][1]=ggv[1][0];
  }
//side 3
  hpvj=hpv1[1];
  hpsj=hps1[1];
  hpssj=hpss1[1];
  for (int i=2;i<n_pols_1d;i++) {
    Tensor<2,2> &ggv=gg[iv++];
    ggv[0][0]=hpss0[i]*hpvj;
    ggv[1][0]=hps0[i]*hpsj;
    ggv[1][1]=hpv0[i]*hpssj;
    ggv[0][1]=ggv[1][0];
  }
//interior
  for (int j=2;j<n_pols_1d;j++) {
    const REAL &hpvj=hpv1[j];
    const REAL &hpsj=hps1[j];
    const REAL &hpssj=hpss1[j];
    for (int i=2;i<n_pols_1d;i++) {
      Tensor<2,2> &ggv=gg[iv++];
      ggv[0][0]=hpss0[i]*hpvj;
      ggv[1][0]=hps0[i]*hpsj;
      ggv[1][1]=hpv0[i]*hpssj;
      ggv[0][1]=ggv[1][0];
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
void HierarchicalQuadrilateralPolynomials::shapeFunctionsOrder(
const Shape *element,NumPtr<int> &renumber) const {
  CHECK_POINTER(dynamic_cast<const Quadrilateral*>(element));
  renumber.allocate(getNumber());
  for (int i=0;i<renumber.getNumber();i++) { renumber[i]=i; }
  if (n_pols_1d<=2) return;
  int order=n_pols_1d-1;
  int nl=order-1;
  int a=3;
//side 0
  if (element->vertexIndex(1)<element->vertexIndex(0)) {
    int s=1;
    for (int j=1;j<order;j++) {
      renumber[a+j]=s*(a+j);
      s=-s;
    }
  }
//side 1
  a+=nl;
  if (element->vertexIndex(3)<element->vertexIndex(1)) {
    int s=1;
    for (int j=1;j<order;j++) {
      renumber[a+j]=s*(a+j);
      s=-s;
    }
  }
//side 2
  a+=nl;
  if (element->vertexIndex(1)<element->vertexIndex(0)) {
    int s=1;
    for (int i=1;i<order;i++) {
      renumber[a+i]=s*(a+i);
      s=-s;
    }
  }
//side 3
  a+=nl;
  if (element->vertexIndex(3)<element->vertexIndex(2)) {
    int s=1;
    for (int i=1;i<order;i++) {
      renumber[a+i]=s*(a+i);
      s=-s;
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
void HierarchicalQuadrilateralPolynomials::nonzeroShapesOnFace(
const Shape1 *e,int f,NumPtr<int> &indices) const {
  indices.allocate(n_pols_1d);
  int order=n_pols_1d-1;
  int iv=0;
  for (int v=0;v<2;v++) indices[iv++]=e->elementVertexFromFace(f,v);
  int offset=e->numberVertices()+f*(order-1)-1;
  for (int i=1;i<order;i++) indices[iv++]=offset+i;
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM>1)
void HierarchicalQuadrilateralPolynomials::printOn(ostream &os) const {
  os << "HierarchicalQuadrilateralPolynomials: npols_1d = " << n_pols_1d
     << "\n\tn_pols = " << n_pols << endl;
  ShapeFunction<2>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
/*
generalizing Peano, Comp. & Maths. with Appls 2(1976) 211-224
the following polynomials in the barycentric coordinates (a,b,c,d) form 
a basis for polynomials of 3 variables of order n:
  1
  if 0<n: a, b, c
  for 1<=k<n: d^k*a, a^k*b, d^k*b, d^k*c, a^k*c, b^k*c
  for 3<=3m<=n: (d*b*c)^m, (d*c*a)^m, (d*a*b)^m, (a*c*b)^m 
  for 4<=3m+1<=n: d*(d*b*c)^m, b*(d*b*c)^m, d*(d*c*a)^m, c*(d*c*a)^m,
                  d*(d*a*b)^m, a*(d*a*b)^m, a*(a*c*b)^m, c*(a*c*b)^m
  for 4<=4p<=n: (d*a*b*c)^p
  for 5<=3m+2<=n, for 1<=k<n-3m: 
    d^k*c*(d*b*c)^m, d^k*b*(d*b*c)^m, b^k*c*(d*b*c)^m,
    d^k*a*(d*c*a)^m, d^k*c*(d*c*a)^m, a^k*c*(d*c*a)^m,
    d^k*b*(d*a*b)^m, d^k*a*(d*a*b)^m, a^k*b*(d*a*b)^m,
    a^k*b*(a*c*b)^m, a^k*c*(a*c*b)^m, b^k*c*(a*c*b)^m
  for 5<=4*p+1<=n: a*(d*a*b*c)^p, b*(d*a*b*c)^p, c*(d*a*b*c)^p
  for 6<=4*p+2<=n, for 1<=k<n-4*p 
    d^k*a*(d*a*b*c)^p, a^k*b*(d*a*b*c)^p, d^k*b*(d*a*b*c)^p,
    d^k*c*(d*a*b*c)^p, a^k*c*(d*a*b*c)^p, b^k*c*(d*a*b*c)^p
  for 7<=4*p+3<=n, for 3<=3*m<=n-4*p:
    (d*b*c)^m*(d*a*b*c)^p, (d*c*a)^m*(d*a*b*c)^p,
    (d*a*b)^m*(d*a*b*c)^p, (a*c*b)^m*(d*a*b*c)^p
  for 8<=4*p+4<=n, for 4<=3*m+1<=n-4*p
    d*(d*b*c)^m*(d*a*b*c)^p, b*(d*b*c)^m*(d*a*b*c)^p, 
    d*(d*c*a)^m*(d*a*b*c)^p, c*(d*c*a)^m*(d*a*b*c)^p,
    d*(d*a*b)^m*(d*a*b*c)^p, a*(d*a*b)^m*(d*a*b*c)^p, 
    a*(a*c*b)^m*(d*a*b*c)^p, c*(a*c*b)^m*(d*a*b*c)^p
  for 9<=4*p+5<=n, for 5<=3m+2<=n-4*p, for 1<=k<n-4*p-3*m: 
    d^k*c*(d*b*c)^m*(d*a*b*c)^p, d^k*b*(d*b*c)^m*(d*a*b*c)^p, 
                                 b^k*c*(d*b*c)^m*(d*a*b*c)^p,
    d^k*a*(d*c*a)^m*(d*a*b*c)^p, d^k*c*(d*c*a)^m*(d*a*b*c)^p, 
                                 a^k*c*(d*c*a)^m*(d*a*b*c)^p,
    d^k*b*(d*a*b)^m*(d*a*b*c)^p, d^k*a*(d*a*b)^m*(d*a*b*c)^p, 
                                 a^k*b*(d*a*b)^m*(d*a*b*c)^p,
    a^k*b*(a*c*b)^m*(d*a*b*c)^p, a^k*c*(a*c*b)^m*(d*a*b*c)^p, 
                                 b^k*c*(a*c*b)^m*(d*a*b*c)^p
The first 3 entries (1,a,b) are associated with the triangle vertices.
The polynomials a^k*b, b^k*c, c^k*a are associated with the triangle
sides.
The remaining polynomials are associated with the triangle interior.

We can replace barycentric coordinate powers a^k, b^k, c^k, d^k with 
the hierarchical polynomials A_k, B_k, C_k, D_k of the same order.  
We do this so that
  1) hierarchical polynomials agree on sides shared by triangles and 
     quadrilaterals
  2) the ordering of the basis functions corresponds to the ordering in
     DoFHandler::distributeDofsOnCell
Assuming that n>0, this leads to the following basis:
  vertex modes ordered as in Shape.H: 
    D_1, (v0)
    A_1, (v1)
    B_1, (v2)
    C_1  (v3)
  line modes ordered as in Shape.H:
    for 1<=k<n: 
      D_k*A_1, (v0,v1)
      A_k*B_1, (v1,v2)
      D_k*B_1, (v0,v2)
      D_k*C_1, (v0,v3)
      A_k*C_1, (v1,v3)
      B_k*C_1  (v2,v3)
  interior side modes ordered as in Shape.H:
    side (v0,v2,v3):
      for 3<=3m<=n: D_m*B*m*C_m
      for 4<=3m+1<=n: D_{m+1}*B_m*C_m, B_{m+1}*D_m*C_m
      for 5<=3m+2<=n, for 1<=k<n-3m: 
        D_{m+k}*B_m*C_{m+1}, D_{m+k}*B_{m+1}*C_m, D_m*B_{m+k}*C_{m+1}
    side (v0,v3,v1):
      for 3<=3m<=n: D_m*A_m*C_m
      for 4<=3m+1<=n: D_{m+1}*A_m*C_m, D_m*A_{m+1}*C_m
      for 5<=3m+2<=n, for 1<=k<n-3m: 
        D_{m+k}*A_{m+1}*C_m, D_{m+k}*A_m*C_{m+1}, D_m*A_{m+k}*C_{m+1}
    side (v0,v1,v2):
      for 3<=3m<=n: D_m*A_m*B_m
      for 4<=3m+1<=n: D_{m+1}*A_m*B_m, D_m*A_{m+1}*B_m 
      for 5<=3m+2<=n, for 1<=k<n-3m:
        D_{m+k}*A_m*B_{m+1}, D_{k+m}*A_{m+1}*B_m, D_m*A_{m+k}*B_{m+1}
    side (v1,v3,v2):
      for 3<=3m<=n: A_m*B_m*C_m
      for 4<=3m+1<=n: A_{m+1}*B_m*C_m, A_m*B_{m+1}*C_m
      for 5<=3m+2<=n, for 1<=k<n-3m:
        A_{m+k}*B_{m+1}*C_m, A_{m+k}*B_m*C_{m+1}, A_m*B_{m+k}*C_{m+1}
  interior modes (order unimportant):
    for 4<=4p<=n: 
      D_p*A_p*B_p*C_p
    for 5<=4*p+1<=n: 
      D_{p+1}*A_p*B_p*C_p, 
      D_p*A_{p+1}*B_p*C_p, 
      D_p*A_p*B_{p+1}*C_p
    for 6<=4*p+2<=n, for 1<=k<n-4*p 
      D_{p+k}*A_{p+1}*B_p*C_p, 
      D_p*A_{p+k}*B_{p+1}*C_p 
      D_{p+k}*A_p*B_{p+1}*C_p, 
      D_{p+k}*A_p*B_p*C_{p+1}, 
      D_p*A_{p+k}*B_p*C_{p+1}, 
      D_p*A_p*B_{p+k}*C_{p+1},
    for 7<=4*p+3<=n, for 3<=3*m<=n-4*p:
      D_{p+m}*A_p*B_{p+m}*C_{p+m},
      D_{p+m}*A_{p+m}*B_p*C_{p+m},
      D_{p+m}*A_{p+m}*B_{p+m}*C_p,
      D_p*A_{p+m}*B_{p+m}*C_{p+m}
    for 8<=4*p+4<=n, for 4<=3*m+1<=n-4*p
      D_{p+m+1}*A_p*B_{p+m}*C_{m+p},
      D_{p+m}*A_p*B_{p+m+1}*C_{m+p},
      D_{p+m+1}*A_{p+m}*B_p*C_{m+p},
      D_{p+m}*A_{p+m+1}*B_p*C_{m+p},
      D_{p+m+1}*A_{p+m}*B_{p+m}*C_p,
      D_{p+m}*A_{p+m+1}*B_{p+m}*C_p,
      D_p*A_{p+m+1}*B_{p+m}*C_{p+m},
      D_p*A_{p+m}*B_{p+m+1}*C_{p+m}
    for 9<=4*p+5<=n, for 5<=3m+2<=n-4*p, for 1<=k<n-4*p-3*m: 
      D_{p+m+k}*A_p*B_{p+m}*C_{p+m+1},
      D_{p+m+k}*A_p*B_{p+m+1}*C_{p+m},
      D_{p+m}*A_p*B_{p+m+k}*C_{p+m+1},
      D_{p+m+k}*A_{p+m+1}*B_p*C_{p+m},
      D_{p+m+k}*A_{p+m}*B_p*C_{p+m+1},
      D_{p+m}*A_{p+m+k}*B_p*C_{p+m+1},
      D_{p+m+k}*A_{p+m}*B_{p+m+1}*C_p,
      D_{p+m+k}*A_{p+m+1}*B_{p+m}*C_p,
      D_{p+m}*A_{p+m+k}*B_{p+m+1}*C_p,
      D_p*A_{p+m+k}*B_{p+m+1}*C_{p+m},
      D_p*A_{p+m+k}*B_{p+m}*C_{p+m+1},
      D_p*A_{p+m}*B_{p+m+k}*C_{p+m+1}
*/
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalTetrahedronPolynomials::computeValues(const Point<3> &p,
NumPtr<REAL> &v) const {
  CHECK_TEST(order<=1);
  v.allocate(4);
//vertices only
  v[0]=1.-(p[0]+p[1]+p[2]);
  v[1]=p[0];
  v[2]=p[1];
  v[3]=p[2];
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalTetrahedronPolynomials::computeGrads(const Point<3> &p,
NumPtr<Tensor<1,3> > &g) const {
  CHECK_TEST(order<=1);
  g.allocate(4);
  int ijk=0;
//vertices only
  Tensor<1,3> &gv0=g[ijk++];
  gv0[0]=-1.;
  gv0[1]=-1.;
  gv0[2]=-1.;
  Tensor<1,3> &gv1=g[ijk++];
  gv1[0]=1.;
  gv1[1]=0.;
  gv1[2]=0.;
  Tensor<1,3> &gv2=g[ijk++];
  gv2[0]=0.;
  gv2[1]=1.;
  gv2[2]=0.;
  Tensor<1,3> &gv3=g[ijk++];
  gv3[0]=0.;
  gv3[1]=0.;
  gv3[2]=1.;
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalTetrahedronPolynomials::computeGradGrads(
const Point<3> &p,NumPtr<Tensor<2,3> > &gg) const {
  CHECK_TEST(order<=1);
  gg.allocate(4);
  for (int ijk=0;ijk<4;ijk++) gg[ijk].initialize(0.);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalTetrahedronPolynomials::computeValuesFor(const Shape *e,
const Point<SPACEDIM> &p,NumPtr<REAL> &v) const {
  int n_polys=getNumber();
  v.allocate(n_polys);
  Point<4> bc(p[0],p[1],p[2],1.-(p[0]+p[1]+p[2]));
  NumPtr<int> ev(e->numberVertices());
  for (int i=0;i<e->numberVertices();i++) ev[i]=e->vertexIndex(i);

  int ijk=0;
//vertices
  v[ijk++]=bc[3];
  v[ijk++]=bc[0];
  v[ijk++]=bc[1];
  v[ijk++]=bc[2];

  if (order<2) return;
  HierarchicalSidePolynomial hsp;
  NumPtr<REAL> hspv0(order-1);
  NumPtr<REAL> hspv1(order-1);
  NumPtr<REAL> hspv2(order-1);
  NumPtr<REAL> hspv3(order-1);
  hsp.values(bc[0],hspv0);
  hsp.values(bc[1],hspv1);
  hsp.values(bc[2],hspv2);
  hsp.values(bc[3],hspv3);
//line 0: v0 -> v1
  NumPtr<REAL> *hspv=(ev[0]<ev[1] ? &hspv0 : &hspv3);
  REAL ab=bc[3]*bc[0];
  for (int i=1;i<order;i++) {
    v[ijk++]=ab*(*hspv)[i-1];
  }
//line 1: v1 -> v2
  hspv=(ev[1]<ev[2] ? &hspv1 : &hspv0);
  ab=bc[0]*bc[1];
  for (int j=1;j<order;j++) {
    v[ijk++]=ab*(*hspv)[j-1];
  }
//line 2: v0 -> v2
  hspv=(ev[0]<ev[2] ? &hspv1 : &hspv3);
  ab=bc[1]*bc[3];
  for (int j=1;j<order;j++) {
    v[ijk++]=ab*(*hspv)[j-1];
  }
//line 3: v0 -> v3
  hspv=(ev[0]<ev[3] ? &hspv2 : &hspv3);
  ab=bc[2]*bc[3];
  for (int k=1;k<order;k++) {
    v[ijk++]=ab*(*hspv)[k-1];
  }
//line 4: v1 -> v3
  hspv=(ev[1]<ev[3] ? &hspv2 : &hspv0);
  ab=bc[0]*bc[2];
  for (int k=1;k<order;k++) {
    v[ijk++]=ab*(*hspv)[k-1];
  }
//line 5: v2 -> v3
  hspv=(ev[2]<ev[3] ? &hspv2 : &hspv1);
  ab=bc[1]*bc[2];
  for (int k=1;k<order;k++) {
    v[ijk++]=ab*(*hspv)[k-1];
  }

  if (order<3) return;
  LegendrePolynomial lp;
  NumPtr<REAL> lpv0(order-2);
  NumPtr<REAL> lpv1(order-2);
  NumPtr<REAL> lpv2(order-2);
  NumPtr<REAL> lpv3(order-2);
  lp.values(bc[0],lpv0);
  lp.values(bc[1],lpv1);
  lp.values(bc[2],lpv2);
  lp.values(bc[3],lpv3);
//face 0: v0, v2, v3
  REAL abc=bc[3]*bc[1]*bc[2];
  NumPtr<REAL> *lpva=0;
  NumPtr<REAL> *lpvb=0;
  if (ev[0]<ev[2]) {
    if (ev[2]<ev[3]) {
      lpva=&lpv1;
      lpvb=&lpv2;
    } else if (ev[0]<ev[3]) {
      lpva=&lpv2;
      lpvb=&lpv1;
    } else {
      lpva=&lpv3;
      lpvb=&lpv1;
    }
  } else if (ev[0]<ev[3]) {
    lpva=&lpv3;
    lpvb=&lpv2;
  } else if (ev[2]<ev[3]) {
    lpva=&lpv2;
    lpvb=&lpv3;
  } else {
    lpva=&lpv1;
    lpvb=&lpv3;
  }
  for (int m=0;m+3<=order;m++) {
    for (int j=0;j<=m;j++) {
      v[ijk++]=abc*(*lpva)[m-j]*(*lpvb)[j];
    }
  }
//face 1: v0, v1, v3
  abc=bc[3]*bc[0]*bc[2];
  if (ev[0]<ev[1]) {
    if (ev[1]<ev[3]) {
      lpva=&lpv0;
      lpvb=&lpv2;
    } else if (ev[0]<ev[3]) {
      lpva=&lpv2;
      lpvb=&lpv0;
    } else {
      lpva=&lpv3;
      lpvb=&lpv0;
    }
  } else if (ev[0]<ev[3]) {
    lpva=&lpv3;
    lpvb=&lpv2;
  } else if (ev[1]<ev[3]) {
    lpva=&lpv2;
    lpvb=&lpv3;
  } else {
    lpva=&lpv0;
    lpvb=&lpv3;
  }
  for (int m=0;m+3<=order;m++) {
    for (int j=0;j<=m;j++) {
      v[ijk++]=abc*(*lpva)[m-j]*(*lpvb)[j];
    }
  }
//face 2: v0, v1, v2
  abc=bc[3]*bc[0]*bc[1];
  if (ev[0]<ev[1]) {
    if (ev[1]<ev[2]) {
      lpva=&lpv0;
      lpvb=&lpv1;
    } else if (ev[0]<ev[2]) {
      lpva=&lpv1;
      lpvb=&lpv0;
    } else {
      lpva=&lpv3;
      lpvb=&lpv0;
    }
  } else if (ev[0]<ev[2]) {
    lpva=&lpv3;
    lpvb=&lpv1;
  } else if (ev[1]<ev[2]) {
    lpva=&lpv1;
    lpvb=&lpv3;
  } else {
    lpva=&lpv0;
    lpvb=&lpv3;
  }
  for (int m=0;m+3<=order;m++) {
    for (int j=0;j<=m;j++) {
      v[ijk++]=abc*(*lpva)[m-j]*(*lpvb)[j];
    }
  }
//face 3: v1, v2, v3
  abc=bc[0]*bc[1]*bc[2];
  if (ev[1]<ev[2]) {
    if (ev[2]<ev[3]) {
      lpva=&lpv1;
      lpvb=&lpv2;
    } else if (ev[1]<ev[3]) {
      lpva=&lpv2;
      lpvb=&lpv1;
    } else {
      lpva=&lpv0;
      lpvb=&lpv1;
    }
  } else if (ev[1]<ev[3]) {
    lpva=&lpv0;
    lpvb=&lpv2;
  } else if (ev[2]<ev[3]) {
    lpva=&lpv2;
    lpvb=&lpv0;
  } else {
    lpva=&lpv1;
    lpvb=&lpv0;
  }
  for (int m=0;m+3<=order;m++) {
    for (int j=0;j<=m;j++) {
      v[ijk++]=abc*(*lpva)[m-j]*(*lpvb)[j];
    }
  }
  if (order<4) return;
//interior
  REAL abcd=bc[0]*bc[1]*bc[2]*bc[3];
  for (int m=0;m+4<=order;m++) {
    for (int k=0;k<=m;k++) {
      for (int j=0;j+k<=m;j++) {
        v[ijk++]=abcd*lpv0[m-j-k]*lpv1[j]*lpv2[k];
      }
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalTetrahedronPolynomials::computeGradsFor(const Shape *e,
const Point<SPACEDIM> &p,NumPtr<Tensor<1,SPACEDIM> > &g) const {
  int n_polys=getNumber();
  g.allocate(n_polys);
  Point<4> bc(p[0],p[1],p[2],1.-(p[0]+p[1]+p[2]));
  NumPtr<int> ev(e->numberVertices());
  for (int i=0;i<e->numberVertices();i++) ev[i]=e->vertexIndex(i);

  int ijk=0;
  Tensor<1,3> &gv0=g[ijk++];
  gv0[0]=-1.;
  gv0[1]=-1.;
  gv0[2]=-1.;
  Tensor<1,3> &gv1=g[ijk++];
  gv1[0]=1.;
  gv1[1]=0.;
  gv1[2]=0.;
  Tensor<1,3> &gv2=g[ijk++];
  gv2[0]=0.;
  gv2[1]=1.;
  gv2[2]=0.;
  Tensor<1,3> &gv3=g[ijk++];
  gv3[0]=0.;
  gv3[1]=0.;
  gv3[2]=1.;

  if (order<2) return;
  HierarchicalSidePolynomial hsp;
  NumPtr<REAL> hspv0(order-1);
  NumPtr<REAL> hspv1(order-1);
  NumPtr<REAL> hspv2(order-1);
  NumPtr<REAL> hspv3(order-1);
  hsp.values(bc[0],hspv0);
  hsp.values(bc[1],hspv1);
  hsp.values(bc[2],hspv2);
  hsp.values(bc[3],hspv3);
  NumPtr<REAL> hsps0(order-1);
  NumPtr<REAL> hsps1(order-1);
  NumPtr<REAL> hsps2(order-1);
  NumPtr<REAL> hsps3(order-1);
  hsp.slopes(bc[0],hsps0);
  hsp.slopes(bc[1],hsps1);
  hsp.slopes(bc[2],hsps2);
  hsp.slopes(bc[3],hsps3);
//line 0: v0 -> v1
  REAL ab=bc[3]*bc[0];
  Tensor<1,3> gab;
  gab[0]=bc[3]-bc[0];
  gab[1]=-bc[0];
  if (ev[0]<ev[1]) {
    for (int i=1;i<order;i++) {
      Tensor<1,3> &gv=g[ijk++];
      gv[0]=gab[0]*hspv0[i-1]+ab*hsps0[i-1];
      gv[1]=gab[1]*hspv0[i-1];
      gv[2]=gv[1];
    }
  } else { // never called:
    for (int i=1;i<order;i++) {
      Tensor<1,3> &gv=g[ijk++];
      gv[0]=gab[0]*hspv3[i-1]-ab*hsps3[i-1];
      gv[1]=gab[1]*hspv3[i-1]-ab*hsps3[i-1];
      gv[2]=gv[1];
    }
  }
//line 1: v1 -> v2
  ab=bc[0]*bc[1];
  gab[0]=bc[1];
  gab[1]=bc[0];
  if (ev[1]<ev[2]) {
    for (int j=1;j<order;j++) {
      Tensor<1,3> &gv=g[ijk++];
      gv[0]=gab[0]*hspv1[j-1];
      gv[1]=gab[1]*hspv1[j-1]+ab*hsps1[j-1];
      gv[2]=0.;
    }
  } else { // never called:
    for (int j=1;j<order;j++) {
      Tensor<1,3> &gv=g[ijk++];
      gv[0]=gab[0]*hspv0[j-1]+ab*hsps0[j-1];
      gv[1]=gab[1]*hspv0[j-1];
      gv[2]=0.;
    }
  }
//line 2: v0 -> v2
  ab=bc[1]*bc[3];
  gab[0]=-bc[1];
  gab[1]=bc[3]-bc[1];
  if (ev[0]<ev[2]) {
    for (int j=1;j<order;j++) {
      Tensor<1,3> &gv=g[ijk++];
      gv[0]=gab[0]*hspv1[j-1];
      gv[1]=gab[1]*hspv1[j-1]+ab*hsps1[j-1];
      gv[2]=gv[0];
    }
  } else { // never called:
    for (int j=1;j<order;j++) {
      Tensor<1,3> &gv=g[ijk++];
      gv[0]=gab[0]*hspv3[j-1]+ab*hsps3[j-1];
      gv[1]=gab[1]*hspv3[j-1]-ab*hsps3[j-1];
      gv[2]=gv[0];
    }
  }
//line 3: v0 -> v3
  ab=bc[2]*bc[3];
  gab[0]=-bc[2];
  gab[2]=bc[3]-bc[2];
  if (ev[0]<ev[3]) {
    for (int k=1;k<order;k++) {
      Tensor<1,3> &gv=g[ijk++];
      gv[0]=gab[0]*hspv2[k-1];
      gv[1]=gv[0];
      gv[2]=gab[2]*hspv2[k-1]+ab*hsps2[k-1];
    }
  } else { // never called:
    for (int k=1;k<order;k++) {
      Tensor<1,3> &gv=g[ijk++];
      gv[0]=gab[0]*hspv3[k-1]+ab*hsps3[k-1];
      gv[1]=gv[0];
      gv[2]=gab[2]*hspv3[k-1]-ab*hsps3[k-1];
    }
  }
//line 4: v1 -> v3
  ab=bc[0]*bc[2];
  gab[0]=bc[2];
  gab[2]=bc[0];
  if (ev[1]<ev[3]) {
    for (int k=1;k<order;k++) {
      Tensor<1,3> &gv=g[ijk++];
      gv[0]=gab[0]*hspv2[k-1];
      gv[1]=0.;
      gv[2]=gab[2]*hspv2[k-1]+ab*hsps2[k-1];
    }
  } else { // never called:
    for (int k=1;k<order;k++) {
      Tensor<1,3> &gv=g[ijk++];
      gv[0]=gab[0]*hspv0[k-1]+ab*hsps0[k-1];
      gv[1]=0.;
      gv[2]=gab[2]*hspv0[k-1];
    }
  }
//line 5: v2 -> v3
  ab=bc[1]*bc[2];
  gab[1]=bc[2];
  gab[2]=bc[1];
  if (ev[2]<ev[3]) {
    for (int k=1;k<order;k++) {
      Tensor<1,3> &gv=g[ijk++];
      gv[0]=0.;
      gv[1]=gab[1]*hspv2[k-1];
      gv[2]=gab[2]*hspv2[k-1]+ab*hsps2[k-1];
    }
  } else {
    for (int k=1;k<order;k++) {
      Tensor<1,3> &gv=g[ijk++];
      gv[0]=0.;
      gv[1]=gab[1]*hspv1[k-1]+ab*hsps1[k-1];
      gv[2]=gab[2]*hspv1[k-1];
    }
  }

  if (order<3) return;
  LegendrePolynomial lp;
  NumPtr<REAL> lpv0(order-2);
  NumPtr<REAL> lpv1(order-2);
  NumPtr<REAL> lpv2(order-2);
  NumPtr<REAL> lpv3(order-2);
  lp.values(bc[0],lpv0);
  lp.values(bc[1],lpv1);
  lp.values(bc[2],lpv2);
  lp.values(bc[3],lpv3);
  NumPtr<REAL> lps0(order-2);
  NumPtr<REAL> lps1(order-2);
  NumPtr<REAL> lps2(order-2);
  NumPtr<REAL> lps3(order-2);
  lp.slopes(bc[0],lps0);
  lp.slopes(bc[1],lps1);
  lp.slopes(bc[2],lps2);
  lp.slopes(bc[3],lps3);
//face 0: v0, v2, v3
  REAL abc=bc[3]*bc[1]*bc[2];
  Tensor<1,3> gabc;
  gabc[0]=-bc[1]*bc[2];
  gabc[1]=bc[2]*(bc[3]-bc[1]);
  gabc[2]=bc[1]*(bc[3]-bc[2]);
  if (ev[0]<ev[2]) {
    if (ev[2]<ev[3]) {
      for (int m=0;m+3<=order;m++) {
        for (int j=0;j<=m;j++) {
          Tensor<1,3> &gv=g[ijk++];
          gv[0]=gabc[0]*lpv1[m-j]*lpv2[j];
          gv[1]=(gabc[1]*lpv1[m-j]+abc*lps1[m-j])*lpv2[j];
          gv[2]=(gabc[2]*lpv2[j]+abc*lps2[j])*lpv1[m-j];
        }
      }
    } else if (ev[0]<ev[3]) {
      for (int m=0;m+3<=order;m++) {
        for (int j=0;j<=m;j++) {
          Tensor<1,3> &gv=g[ijk++];
          gv[0]=gabc[0]*lpv2[m-j]*lpv1[j];
          gv[1]=(gabc[1]*lpv1[j]+abc*lps1[j])*lpv2[m-j];
          gv[2]=(gabc[2]*lpv2[m-j]+abc*lps2[m-j])*lpv1[j];
        }
      }
    } else { // never called:
      for (int m=0;m+3<=order;m++) {
        for (int j=0;j<=m;j++) {
          Tensor<1,3> &gv=g[ijk++];
          gv[0]=(gabc[0]*lpv3[m-j]-abc*lps3[m-j])*lpv1[j];
          gv[1]=gabc[1]*lpv3[m-j]*lpv1[j]
               +abc*(-lps3[m-j]*lpv1[j]+lpv3[m-j]*lps1[j]);
          gv[2]=(gabc[2]*lpv3[m-j]-abc*lps3[m-j])*lpv1[j];
        }
      }
    }
  } else if (ev[0]<ev[3]) { // never called:
    for (int m=0;m+3<=order;m++) {
      for (int j=0;j<=m;j++) {
        Tensor<1,3> &gv=g[ijk++];
        gv[0]=(gabc[0]*lpv3[m-j]-abc*lps3[m-j])*lpv2[j];
        gv[1]=(gabc[1]*lpv3[m-j]-abc*lps3[m-j])*lpv2[j];
        gv[2]=gabc[2]*lpv3[m-j]*lpv2[j]
             +abc*(-lps3[m-j]*lpv2[j]+lpv3[m-j]*lps2[j]);
      }
    }
  } else if (ev[2]<ev[3]) { // never called:
    for (int m=0;m+3<=order;m++) {
      for (int j=0;j<=m;j++) {
        Tensor<1,3> &gv=g[ijk++];
        gv[0]=(gabc[0]*lpv2[m-j]-abc*lpv2[m-j])*lps3[j];
        gv[1]=(gabc[1]*lpv2[m-j]-abc*lpv2[m-j])*lps3[j];
        gv[2]=gabc[2]*lpv2[m-j]*lpv3[j]
             +abc*(lps2[m-j]*lpv3[j]+lpv2[m-j]*lps3[j]);
      }
    }
  } else { // never called:
    for (int m=0;m+3<=order;m++) {
      for (int j=0;j<=m;j++) {
        Tensor<1,3> &gv=g[ijk++];
        gv[0]=(gabc[0]*lpv1[m-j]-abc*lpv1[m-j])*lps3[j];
        gv[1]=gabc[1]*lpv1[m-j]*lpv3[j]
             +abc*(lps1[m-j]*lpv3[j]-lpv1[m-j]*lps3[j]);
        gv[2]=(gabc[2]*lpv1[m-j]-abc*lpv1[m-j])*lps3[j];
      }
    }
  }
//face 1: v0, v1, v3
  abc=bc[3]*bc[0]*bc[2];
  gabc[0]=bc[2]*(bc[3]-bc[0]);
  gabc[1]=-bc[0]*bc[2];
  gabc[2]=bc[0]*(bc[3]-bc[2]);
  if (ev[0]<ev[1]) {
    if (ev[1]<ev[3]) {
      for (int m=0;m+3<=order;m++) {
        for (int j=0;j<=m;j++) {
          Tensor<1,3> &gv=g[ijk++];
          gv[0]=(gabc[0]*lpv0[m-j]+abc*lps0[m-j])*lpv2[j];
          gv[1]=gabc[1]*lpv0[m-j]*lpv2[j];
          gv[2]=(gabc[2]*lpv2[j]+abc*lps2[j])*lpv0[m-j];
        }
      }
    } else if (ev[0]<ev[3]) { // never called:
      for (int m=0;m+3<=order;m++) {
        for (int j=0;j<=m;j++) {
          Tensor<1,3> &gv=g[ijk++];
          gv[0]=(gabc[0]*lpv0[j]+abc*lps0[j])*lpv2[m-j];
          gv[1]=gabc[1]*lpv2[m-j]*lpv0[j];
          gv[0]=(gabc[0]*lpv2[m-j]+abc*lps2[m-j])*lpv0[j];
        }
      }
    } else { // never called:
      for (int m=0;m+3<=order;m++) {
        for (int j=0;j<=m;j++) {
          Tensor<1,3> &gv=g[ijk++];
          gv[0]=gabc[0]*lpv3[m-j]*lpv0[j]
               +abc*(-lps3[m-j]*lpv0[j]+lpv3[m-j]*lps0[j]);
          gv[1]=(gabc[1]*lpv3[m-j]-abc*lps3[m-j])*lpv0[j];
          gv[2]=(gabc[2]*lpv3[m-j]-abc*lps3[m-j])*lpv0[j];
        }
      }
    }
  } else if (ev[0]<ev[3]) { // never called:
    for (int m=0;m+3<=order;m++) {
      for (int j=0;j<=m;j++) {
        Tensor<1,3> &gv=g[ijk++];
        gv[0]=(gabc[0]*lpv3[m-j]-abc*lps3[m-j])*lpv2[j];
        gv[1]=(gabc[1]*lpv3[m-j]-abc*lps3[m-j])*lpv2[j];
        gv[2]=gabc[2]*lpv3[m-j]*lpv2[j]
             +abc*(-lps3[m-j]*lpv2[j]+lpv3[m-j]*lps2[j]);
      }
    }
  } else if (ev[1]<ev[3]) { // never called:
    for (int m=0;m+3<=order;m++) {
      for (int j=0;j<=m;j++) {
        Tensor<1,3> &gv=g[ijk++];
        gv[0]=(gabc[0]*lpv3[j]-abc*lps3[j])*lpv2[m-j];
        gv[1]=(gabc[1]*lpv3[j]-abc*lps3[j])*lpv2[m-j];
        gv[2]=gabc[2]*lpv2[m-j]*lpv3[j]
             +abc*(lps2[m-j]*lpv3[j]-lpv2[m-j]*lps3[j]);
      }
    }
  } else { // never called:
    for (int m=0;m+3<=order;m++) {
      for (int j=0;j<=m;j++) {
        Tensor<1,3> &gv=g[ijk++];
        gv[0]=gabc[0]*lpv0[m-j]*lpv3[j]
             +abc*(lps0[m-j]*lpv3[j]-lpv0[m-j]*lps3[j]);
        gv[1]=(gabc[1]*lpv3[j]-abc*lps3[j])*lpv0[m-j];
        gv[2]=(gabc[2]*lpv3[j]-abc*lps3[j])*lpv0[m-j];
      }
    }
  }
//face 2: v0, v1, v2
  abc=bc[3]*bc[0]*bc[1];
  gabc[0]=bc[1]*(bc[3]-bc[0]);
  gabc[1]=bc[0]*(bc[3]-bc[1]);
  gabc[2]=-bc[0]*bc[1];
  if (ev[0]<ev[1]) {
    if (ev[1]<ev[2]) {
      for (int m=0;m+3<=order;m++) {
        for (int j=0;j<=m;j++) {
          Tensor<1,3> &gv=g[ijk++];
          gv[0]=(gabc[0]*lpv0[m-j]+abc*lps0[m-j])*lpv1[j];
          gv[1]=(gabc[1]*lpv1[j]+abc*lps1[j])*lpv0[m-j];
          gv[2]=gabc[2]*lpv0[m-j]*lpv1[j];
        }
      }
    } else if (ev[0]<ev[2]) { // never called:
      for (int m=0;m+3<=order;m++) {
        for (int j=0;j<=m;j++) {
          Tensor<1,3> &gv=g[ijk++];
          gv[0]=(gabc[0]*lpv0[j]+abc*lps0[j])*lpv1[m-j];
          gv[1]=(gabc[1]*lpv1[m-j]+abc*lps1[m-j])*lpv0[j];
          gv[2]=gabc[2]*lpv0[m-j]*lpv1[j];
        }
      }
    } else { // never called:
      for (int m=0;m+3<=order;m++) {
        for (int j=0;j<=m;j++) {
          Tensor<1,3> &gv=g[ijk++];
          gv[0]=gabc[0]*lpv3[m-j]*lpv0[j]
               +abc*(-lps3[m-j]*lpv0[j]+lpv3[m-j]*lps0[j]);
          gv[1]=(gabc[1]*lpv3[m-j]-abc*lps3[m-j])*lpv0[j];
          gv[2]=(gabc[2]*lpv3[m-j]-abc*lps3[m-j])*lpv0[j];
        }
      }
    }
  } else if (ev[0]<ev[2]) { // never called:
    for (int m=0;m+3<=order;m++) {
      for (int j=0;j<=m;j++) {
        Tensor<1,3> &gv=g[ijk++];
        gv[0]=(gabc[0]*lpv3[m-j]-abc*lps3[m-j])*lpv1[j];
        gv[1]=gabc[1]*lpv3[m-j]*lpv1[j]
             +abc*(-lps3[m-j]*lpv1[j]+lpv3[m-j]*lps1[j]);
        gv[2]=(gabc[2]*lpv3[m-j]-abc*lps3[m-j])*lpv1[j];
      }
    }
  } else if (ev[1]<ev[2]) { // never called:
    for (int m=0;m+3<=order;m++) {
      for (int j=0;j<=m;j++) {
        Tensor<1,3> &gv=g[ijk++];
        gv[0]=(gabc[0]*lpv3[j]-abc*lps3[j])*lpv1[m-j];
        gv[1]=gabc[1]*lpv1[m-j]*lpv3[j]
             +abc*(lps1[m-j]*lpv3[j]-lpv1[m-j]*lps3[j]);
        gv[2]=(gabc[2]*lpv3[j]-abc*lps3[j])*lpv1[m-j];
      }
    }
  } else { // never called:
    for (int m=0;m+3<=order;m++) {
      for (int j=0;j<=m;j++) {
        Tensor<1,3> &gv=g[ijk++];
        gv[0]=gabc[0]*lpv0[m-j]*lpv3[j]
             +abc*(lps0[m-j]*lpv3[j]-lpv0[m-j]*lps3[j]);
        gv[1]=(gabc[1]*lpv3[j]-abc*lps3[j])*lpv0[m-j];
        gv[2]=(gabc[2]*lpv3[j]-abc*lps3[j])*lpv0[m-j];
      }
    }
  }
//face 3: v1, v2, v3
  abc=bc[0]*bc[1]*bc[2];
  gabc[0]=bc[1]*bc[2];
  gabc[1]=bc[0]*bc[2];
  gabc[2]=bc[0]*bc[1];
  if (ev[1]<ev[2]) {
    if (ev[2]<ev[3]) {
      for (int m=0;m+3<=order;m++) {
        for (int j=0;j<=m;j++) {
          Tensor<1,3> &gv=g[ijk++];
          gv[0]=gabc[0]*lpv1[m-j]*lpv2[j];
          gv[1]=(gabc[1]*lpv1[m-j]+abc*lps1[m-j])*lpv2[j];
          gv[2]=(gabc[2]*lpv2[j]+abc*lps2[j])*lpv1[m-j];
        }
      }
    } else if (ev[1]<ev[3]) {
      for (int m=0;m+3<=order;m++) {
        for (int j=0;j<=m;j++) {
          Tensor<1,3> &gv=g[ijk++];
          gv[0]=gabc[0]*lpv2[m-j]*lpv1[j];
          gv[1]=(gabc[1]*lpv1[j]+abc*lps1[j])*lpv2[m-j];
          gv[2]=(gabc[2]*lpv2[m-j]+abc*lps2[m-j])*lpv1[j];
        }
      }
    } else { // never called:
      for (int m=0;m+3<=order;m++) {
        for (int j=0;j<=m;j++) {
          Tensor<1,3> &gv=g[ijk++];
          gv[0]=(gabc[0]*lpv0[m-j]+abc*lps0[m-j])*lpv1[j];
          gv[1]=(gabc[1]*lpv1[j]+abc*lps1[j])*lpv0[m-j];
          gv[2]=gabc[2]*lpv0[m-j]*lpv1[j];
        }
      }
    }
  } else if (ev[1]<ev[3]) { // never called:
    for (int m=0;m+3<=order;m++) {
      for (int j=0;j<=m;j++) {
        Tensor<1,3> &gv=g[ijk++];
        gv[0]=(gabc[0]*lpv0[m-j]+abc*lps0[m-j])*lpv2[j];
        gv[1]=gabc[1]*lpv0[m-j]*lpv2[j];
        gv[2]=(gabc[2]*lpv2[j]+abc*lps2[j])*lpv0[m-j];
      }
    }
  } else if (ev[2]<ev[3]) { // never called:
    for (int m=0;m+3<=order;m++) {
      for (int j=0;j<=m;j++) {
        Tensor<1,3> &gv=g[ijk++];
        gv[0]=(gabc[0]*lpv0[j]+abc*lps0[j])*lpv2[m-j];
        gv[1]=gabc[1]*lpv2[m-j]*lpv0[j];
        gv[2]=(gabc[2]*lpv2[m-j]+abc*lps2[m-j])*lpv0[j];
      }
    }
  } else { // never called:
    for (int m=0;m+3<=order;m++) {
      for (int j=0;j<=m;j++) {
        Tensor<1,3> &gv=g[ijk++];
        gv[0]=(gabc[0]*lpv0[j]+abc*lps0[j])*lpv1[m-j];
        gv[1]=(gabc[1]*lpv1[m-j]+abc*lps1[m-j])*lpv0[j];
        gv[2]=gabc[2]*lpv1[m-j]*lpv0[j];
      }
    }
  }
  if (order<4) return;
  REAL abcd=bc[0]*bc[1]*bc[2]*bc[3];
  Tensor<1,3> gabcd;
  gabcd[0]=bc[1]*bc[2]*(bc[3]-bc[0]);
  gabcd[1]=bc[2]*bc[0]*(bc[3]-bc[1]);
  gabcd[2]=bc[0]*bc[1]*(bc[3]-bc[2]);
  for (int m=0;m+4<=order;m++) {
    for (int k=0;k<=m;k++) {
      for (int j=0;j+k<=m;j++) {
        Tensor<1,3> &gv=g[ijk++];
        gv[0]=(gabcd[0]*lpv0[m-j-k]+abcd*lps0[m-j-k])
             *lpv1[j]*lpv2[k];
        gv[1]=(gabcd[1]*lpv1[j]+abcd*lps1[j])*lpv0[m-j-k]*lpv2[k];
        gv[2]=(gabcd[2]*lpv2[k]+abcd*lps2[k])*lpv0[m-j-k]*lpv1[j];
      }
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalTetrahedronPolynomials::computeGradGradsFor(
const Shape *e,const Point<SPACEDIM> &p,
NumPtr<Tensor<2,SPACEDIM> > &gg) const {
  int n_polys=getNumber();
  gg.allocate(n_polys);
  Point<4> bc(p[0],p[1],p[2],1.-(p[0]+p[1]+p[2]));
  NumPtr<int> ev(e->numberVertices());
  for (int i=0;i<e->numberVertices();i++) ev[i]=e->vertexIndex(i);

  int ijk=0;
  for (int v=0;v<4;v++) {
    Tensor<2,3> &ggv=gg[ijk++];
    for (int j=0;j<3;j++) {
      for (int i=0;i<3;i++) {
        ggv[i][j]=0.;
      }
    }
  }
  if (order<2) return;
  HierarchicalSidePolynomial hsp;
  NumPtr<REAL> hspv0(order-1);
  NumPtr<REAL> hspv1(order-1);
  NumPtr<REAL> hspv2(order-1);
  NumPtr<REAL> hspv3(order-1);
  hsp.values(bc[0],hspv0);
  hsp.values(bc[1],hspv1);
  hsp.values(bc[2],hspv2);
  hsp.values(bc[3],hspv3);
  NumPtr<REAL> hsps0(order-1);
  NumPtr<REAL> hsps1(order-1);
  NumPtr<REAL> hsps2(order-1);
  NumPtr<REAL> hsps3(order-1);
  hsp.slopes(bc[0],hsps0);
  hsp.slopes(bc[1],hsps1);
  hsp.slopes(bc[2],hsps2);
  hsp.slopes(bc[3],hsps3);
  NumPtr<REAL> hspss0(order-1);
  NumPtr<REAL> hspss1(order-1);
  NumPtr<REAL> hspss2(order-1);
  NumPtr<REAL> hspss3(order-1);
  hsp.slope2s(bc[0],hspss0);
  hsp.slope2s(bc[1],hspss1);
  hsp.slope2s(bc[2],hspss2);
  hsp.slope2s(bc[3],hspss3);
//line 0: v0 -> v1
  REAL ab=bc[3]*bc[0];
  Tensor<1,3> gab;
  gab[0]=bc[3]-bc[0];
  gab[1]=-bc[0];
  gab[2]=gab[1];
  if (ev[0]<ev[1]) {
    for (int i=1;i<order;i++) {
      Tensor<2,3> &ggv=gg[ijk++];
      ggv[0][0]=-2.*hspv0[i-1]+2.*gab[0]*hsps0[i-1]+ab*hspss0[i-1];
      ggv[1][0]=-hspv0[i-1]+gab[1]*hsps0[i-1];
      ggv[2][0]=ggv[1][0];
      ggv[1][1]=0.;
      ggv[2][1]=0.;
      ggv[2][2]=0.;
      ggv[0][1]=ggv[1][0];
      ggv[0][2]=ggv[2][0];
      ggv[1][2]=ggv[2][1];
    }
  } else { // never called:
    for (int i=1;i<order;i++) {
      Tensor<2,3> &ggv=gg[ijk++];
      ggv[0][0]=-2.*hspv3[i-1]-2.*gab[0]*hsps3[i-1]+ab*hspss3[i-1];
      ggv[1][0]=-hspv3[i-1]-(gab[0]+gab[1])*hsps3[i-1]+ab*hspss3[i-1];
      ggv[2][0]=ggv[1][0];
      ggv[1][1]=-2.*gab[1]*hsps3[i-1]+ab*hspss3[i-1];
      ggv[2][1]=ggv[1][1];
      ggv[2][2]=-2.*gab[2]*hsps3[i-1]+ab*hspss3[i-1];
      ggv[0][1]=ggv[1][0];
      ggv[0][2]=ggv[2][0];
      ggv[1][2]=ggv[2][1];
    }
  }
//line 1: v1 -> v2
  ab=bc[0]*bc[1];
  gab[0]=bc[1];
  gab[1]=bc[0];
  if (ev[1]<ev[2]) {
    for (int j=1;j<order;j++) {
      Tensor<2,3> &ggv=gg[ijk++];
      ggv[0][0]=0.;
      ggv[1][0]=hspv1[j-1]+gab[0]*hsps1[j-1];
      ggv[2][0]=0.;
      ggv[1][1]=2.*gab[1]*hsps1[j-1]+ab*hspss1[j-1];
      ggv[2][1]=0.;
      ggv[2][2]=0.;
      ggv[0][1]=ggv[1][0];
      ggv[0][2]=ggv[2][0];
      ggv[1][2]=ggv[2][1];
    }
  } else { // never called:
    for (int j=1;j<order;j++) {
      Tensor<2,3> &ggv=gg[ijk++];
      ggv[0][0]=2.*gab[0]*hsps0[j-1]+ab*hspss0[j-1];
      ggv[1][0]=hspv0[j-1]+gab[1]*hsps0[j-1];
      ggv[2][0]=0.;
      ggv[1][1]=0.;
      ggv[2][1]=0.;
      ggv[2][2]=0.;
      ggv[0][1]=ggv[1][0];
      ggv[0][2]=ggv[2][0];
      ggv[1][2]=ggv[2][1];
    }
  }
//line 2: v0 -> v2
  ab=bc[1]*bc[3];
  gab[0]=-bc[1];
  gab[1]=bc[3]-bc[1];
  if (ev[0]<ev[2]) {
    for (int j=1;j<order;j++) {
      Tensor<2,3> &ggv=gg[ijk++];
      ggv[0][0]=0.;
      ggv[1][0]=-hspv1[j-1]+gab[0]*hsps1[j-1];
      ggv[2][0]=0.;
      ggv[1][1]=-2.*hspv1[j-1]+2.*gab[1]*hsps1[j-1]+ab*hspss1[j-1];
      ggv[2][1]=ggv[1][0];
      ggv[2][2]=0.;
      ggv[0][1]=ggv[1][0];
      ggv[0][2]=ggv[2][0];
      ggv[1][2]=ggv[2][1];
    }
  } else { // never called:
    for (int j=1;j<order;j++) {
      Tensor<2,3> &ggv=gg[ijk++];
      ggv[0][0]=-2.*gab[0]*hsps3[j-1]+ab*hspss3[j-1];
      ggv[1][0]=-hspv3[j-1]-(gab[1]+gab[0])*hsps3[j-1]+ab*hspss3[j-1];
      ggv[2][0]=ggv[0][0];
      ggv[1][1]=-2.*hspv3[j-1]-2.*gab[1]*hsps3[j-1]+ab*hspss3[j-1];
      ggv[2][1]=ggv[1][0];
      ggv[2][2]=ggv[0][0];
      ggv[0][1]=ggv[1][0];
      ggv[0][2]=ggv[2][0];
      ggv[1][2]=ggv[2][1];
    }
  }
//line 3: v0 -> v3
  ab=bc[2]*bc[3];
  gab[0]=-bc[2];
  gab[2]=bc[3]-bc[2];
  if (ev[0]<ev[3]) {
    for (int k=1;k<order;k++) {
      Tensor<2,3> &ggv=gg[ijk++];
      ggv[0][0]=0.;
      ggv[1][0]=0.;
      ggv[2][0]=-hspv2[k-1]+gab[0]*hsps2[k-1];
      ggv[1][1]=0.;
      ggv[2][1]=ggv[2][0];
      ggv[2][2]=-2.*hspv2[k-1]+2.*gab[2]*hsps2[k-1]+ab*hspss2[k-1];
      ggv[0][1]=ggv[1][0];
      ggv[0][2]=ggv[2][0];
      ggv[1][2]=ggv[2][1];
    }
  } else { // never called:
    for (int k=1;k<order;k++) {
      Tensor<2,3> &ggv=gg[ijk++];
      ggv[0][0]=2.*gab[0]*hsps3[k-1]+ab*hspss3[k-1];
      ggv[1][0]=0.;
      ggv[2][0]=-hspv3[k-1]-gab[2]*hsps3[k-1]+ab*hspss3[k-1];
      ggv[1][1]=ggv[0][0];
      ggv[2][1]=ggv[2][0];
      ggv[2][2]=-2.*hspv3[k-1]-2.*gab[2]*hsps3[k-1]+ab*hspss3[k-1];
      ggv[0][1]=ggv[1][0];
      ggv[0][2]=ggv[2][0];
      ggv[1][2]=ggv[2][1];
    }
  }
//line 4: v1 -> v3
  ab=bc[0]*bc[2];
  gab[0]=bc[2];
  gab[2]=bc[0];
  if (ev[1]<ev[3]) {
    for (int k=1;k<order;k++) {
      Tensor<2,3> &ggv=gg[ijk++];
      ggv[0][0]=0.;
      ggv[1][0]=0.;
      ggv[2][0]=hspv2[k-1]+gab[0]*hsps2[k-1];
      ggv[1][1]=0.;
      ggv[2][1]=0.;
      ggv[2][2]=2.*gab[2]*hsps2[k-1]+ab*hspss2[k-1];
      ggv[0][1]=ggv[1][0];
      ggv[0][2]=ggv[2][0];
      ggv[1][2]=ggv[2][1];
    }
  } else { // never called:
    for (int k=1;k<order;k++) {
      Tensor<2,3> &ggv=gg[ijk++];
      ggv[0][0]=2.*gab[0]*hsps0[k-1]+ab*hspss0[k-1];
      ggv[1][0]=0.;
      ggv[2][0]=hspv0[k-1]+gab[2]*hsps0[k-1];
      ggv[1][1]=0.;
      ggv[2][1]=0.;
      ggv[2][2]=0.;
      ggv[0][1]=ggv[1][0];
      ggv[0][2]=ggv[2][0];
      ggv[1][2]=ggv[2][1];
    }
  }
//line 5: v2 -> v3
  ab=bc[1]*bc[2];
  gab[1]=bc[2];
  gab[2]=bc[1];
  if (ev[2]<ev[3]) {
    for (int k=1;k<order;k++) {
      Tensor<2,3> &ggv=gg[ijk++];
      ggv[0][0]=0.;
      ggv[1][0]=0.;
      ggv[2][0]=0.;
      ggv[1][1]=0.;
      ggv[2][1]=hspv2[k-1]+gab[1]*hsps2[k-1];
      ggv[2][2]=2.*gab[2]*hsps2[k-1]+ab*hspss2[k-1];
      ggv[0][1]=ggv[1][0];
      ggv[0][2]=ggv[2][0];
      ggv[1][2]=ggv[2][1];
    }
  } else {
    for (int k=1;k<order;k++) {
      Tensor<2,3> &ggv=gg[ijk++];
      ggv[0][0]=0.;
      ggv[1][0]=0.;
      ggv[2][0]=0.;
      ggv[1][1]=2.*gab[1]*hsps1[k-1]+ab*hspss1[k-1];
      ggv[2][1]=hspv1[k-1]+gab[2]*hsps1[k-1];
      ggv[2][2]=0.;
      ggv[0][1]=ggv[1][0];
      ggv[0][2]=ggv[2][0];
      ggv[1][2]=ggv[2][1];
    }
  }

  if (order<3) return;
  LegendrePolynomial lp;
  NumPtr<REAL> lpv0(order-2);
  NumPtr<REAL> lpv1(order-2);
  NumPtr<REAL> lpv2(order-2);
  NumPtr<REAL> lpv3(order-2);
  lp.values(bc[0],lpv0);
  lp.values(bc[1],lpv1);
  lp.values(bc[2],lpv2);
  lp.values(bc[3],lpv3);
  NumPtr<REAL> lps0(order-2);
  NumPtr<REAL> lps1(order-2);
  NumPtr<REAL> lps2(order-2);
  NumPtr<REAL> lps3(order-2);
  lp.slopes(bc[0],lps0);
  lp.slopes(bc[1],lps1);
  lp.slopes(bc[2],lps2);
  lp.slopes(bc[3],lps3);
  NumPtr<REAL> lpss0(order-2);
  NumPtr<REAL> lpss1(order-2);
  NumPtr<REAL> lpss2(order-2);
  NumPtr<REAL> lpss3(order-2);
  lp.slope2s(bc[0],lpss0);
  lp.slope2s(bc[1],lpss1);
  lp.slope2s(bc[2],lpss2);
  lp.slope2s(bc[3],lpss3);
//face 0: v0, v2, v3
  REAL abc=bc[3]*bc[1]*bc[2];
  Tensor<1,3> gabc;
  gabc[0]=-bc[1]*bc[2];
  gabc[1]=bc[2]*(bc[3]-bc[1]);
  gabc[2]=bc[1]*(bc[3]-bc[2]);
  Tensor<2,3> ggabc;
  ggabc[0][0]=0.;
  ggabc[1][0]=-bc[2];
  ggabc[2][0]=-bc[1];
  ggabc[1][1]=-2.*bc[2];
  ggabc[2][1]=bc[3]-bc[1]-bc[2];
  ggabc[2][2]=-2.*bc[1];
  ggabc[0][1]=ggabc[1][0];
  ggabc[0][2]=ggabc[2][0];
  ggabc[1][2]=ggabc[2][1];
  if (ev[0]<ev[2]) {
    if (ev[2]<ev[3]) {
      for (int m=0;m+3<=order;m++) {
        for (int j=0;j<=m;j++) {
          Tensor<2,3> &ggv=gg[ijk++];
          ggv[0][0]=ggabc[0][0]*lpv1[m-j]*lpv2[j];
          ggv[1][0]=(ggabc[1][0]*lpv1[m-j]+gabc[0]*lps1[m-j])*lpv2[j];
          ggv[2][0]=(ggabc[2][0]*lpv2[j]+gabc[0]*lps2[j])*lpv1[m-j];
          ggv[1][1]=(ggabc[1][1]*lpv1[m-j]+2.*gabc[1]*lps1[m-j]
                    +abc*lpss1[m-j])*lpv2[j];
          ggv[2][1]=ggabc[2][1]*lpv1[m-j]*lpv2[j]+gabc[2]*lps1[m-j]*lpv2[j]
                   +gabc[1]*lpv1[m-j]*lps2[j]+abc*lps1[m-j]*lps2[j];
          ggv[2][2]=(ggabc[2][2]*lpv2[j]+2.*gabc[2]*lps2[j]
                    +abc*lpss2[j])*lpv1[m-j];
          ggv[0][1]=ggv[1][0];
          ggv[0][2]=ggv[2][0];
          ggv[1][2]=ggv[2][1];
        }
      }
    } else if (ev[0]<ev[3]) {
      for (int m=0;m+3<=order;m++) {
        for (int j=0;j<=m;j++) {
          Tensor<2,3> &ggv=gg[ijk++];
          ggv[0][0]=ggabc[0][0]*lpv2[m-j]*lpv1[j];
          ggv[1][0]=(ggabc[1][0]*lpv1[j]+gabc[0]*lps1[j])*lpv2[m-j];
          ggv[2][0]=(ggabc[2][0]*lpv2[m-j]+gabc[0]*lps2[m-j])*lpv1[j];
          ggv[1][1]=(ggabc[1][1]*lpv1[j]+2.*gabc[1]*lps1[j]
                    +abc*lpss1[j])*lpv2[m-j];
          ggv[2][1]=ggabc[2][1]*lpv2[m-j]*lpv1[j]
                   +gabc[2]*lpv2[m-j]*lps1[j]+gabc[1]*lps2[m-j]*lpv1[j]
                   +abc*lps2[m-j]*lps1[j];
          ggv[2][2]=(ggabc[2][2]*lpv2[m-j]+2.*gabc[2]*lps2[m-j]
                    +abc*lpss2[m-j])*lpv1[j];
          ggv[0][1]=ggv[1][0];
          ggv[0][2]=ggv[2][0];
          ggv[1][2]=ggv[2][1];
        }
      }
    } else { // never called:
      for (int m=0;m+3<=order;m++) {
        for (int j=0;j<=m;j++) {
          Tensor<2,3> &ggv=gg[ijk++];
          ggv[0][0]=(ggabc[0][0]*lpv3[m-j]-2.*gabc[0]*lps3[m-j]
                    +abc*lpss3[m-j])*lpv1[j];
          ggv[1][0]=ggabc[1][0]*lpv3[m-j]*lpv1[j]
                   -(gabc[1]+gabc[0])*lps3[m-j]*lpv1[j]
                   +abc*(lpss3[m-j]*lpv1[j]-lps3[m-j]*lps1[j]);
          ggv[2][0]=(ggabc[2][0]*lpv3[m-j]-(gabc[2]+gabc[0])*lps3[m-j]
                    +abc*lpss3[m-j])*lpv1[j];
          ggv[1][1]=ggabc[1][1]*lpv3[m-j]*lpv1[j]
                   +2.*gabc[1]*(-lps3[m-j]*lpv1[j]+lpv3[m-j]*lps1[j])
                   +abc*(lpss3[m-j]*lpv1[j]-2.*lps3[m-j]*lps1[j]
                        +lpv3[m-j]*lpss1[j]);
          ggv[2][1]=ggabc[2][1]*lpv3[m-j]*lpv1[j]
                   -gabc[2]*(-lps3[m-j]*lpv1[j]+lpv3[m-j]*lps1[j])
                   -gabc[1]*lpv3[m-j]*lps1[j]
                   +abc*(lpss3[m-j]*lpv1[j]-lps3[m-j]*lps1[j]);
          ggv[2][2]=(ggabc[2][2]*lpv3[m-j]-2.*gabc[2]*lps3[m-j]
                    +abc*lpss3[m-j])*lpv1[j];
          ggv[0][1]=ggv[1][0];
          ggv[0][2]=ggv[2][0];
          ggv[1][2]=ggv[2][1];
        }
      }
    }
  } else if (ev[0]<ev[3]) { // never called:
    OBSOLETE("not programmed");
  } else if (ev[2]<ev[3]) { // never called:
    OBSOLETE("not programmed");
  } else { // never called:
    OBSOLETE("not programmed");
  }
//face 1: v0, v1, v3
  abc=bc[3]*bc[0]*bc[2];
  gabc[0]=bc[2]*(bc[3]-bc[0]);
  gabc[1]=-bc[0]*bc[2];
  gabc[2]=bc[0]*(bc[3]-bc[2]);
  ggabc[0][0]=-2.*bc[2];
  ggabc[1][0]=-bc[2];
  ggabc[2][0]=bc[3]-bc[0]-bc[2];
  ggabc[1][1]=0.;
  ggabc[2][1]=-bc[0];
  ggabc[2][2]=-2.*bc[0];
  ggabc[0][1]=ggabc[1][0];
  ggabc[0][2]=ggabc[2][0];
  ggabc[1][2]=ggabc[2][1];
  if (ev[0]<ev[1]) {
    if (ev[1]<ev[3]) {
      for (int m=0;m+3<=order;m++) {
        for (int j=0;j<=m;j++) {
          Tensor<2,3> &ggv=gg[ijk++];
          ggv[0][0]=(ggabc[0][0]*lpv0[m-j]+2.*gabc[0]*lps0[m-j]
                    +abc*lpss0[m-j])*lpv2[j];
          ggv[1][0]=(ggabc[1][0]*lpv0[m-j]+gabc[1]*lps0[m-j])*lpv2[j];
          ggv[2][0]=ggabc[2][0]*lpv0[m-j]*lpv2[j]+gabc[0]*lpv0[m-j]*lps2[j]
                   +gabc[2]*lps0[m-j]*lpv2[j]+abc*lps0[m-j]*lps2[j];
          ggv[1][1]=ggabc[1][1]*lpv0[m-j]*lpv2[j];
          ggv[2][1]=(ggabc[2][1]*lpv2[j]+gabc[1]*lps2[j])*lpv0[m-j];
          ggv[2][2]=(ggabc[2][2]*lpv2[j]+2.*gabc[2]*lps2[j]
                    +abc*lpss2[j])*lpv0[m-j];
          ggv[0][1]=ggv[1][0];
          ggv[0][2]=ggv[2][0];
          ggv[1][2]=ggv[2][1];
        }
      }
    } else if (ev[0]<ev[3]) { // never called:
      OBSOLETE("not programmed");
    } else { // never called:
      OBSOLETE("not programmed");
    }
  } else if (ev[0]<ev[3]) { // never called:
    OBSOLETE("not programmed");
  } else if (ev[1]<ev[3]) { // never called:
    OBSOLETE("not programmed");
  } else { // never called:
    OBSOLETE("not programmed");
  }
//face 2: v0, v1, v2
  abc=bc[3]*bc[0]*bc[1];
  gabc[0]=bc[1]*(bc[3]-bc[0]);
  gabc[1]=bc[0]*(bc[3]-bc[1]);
  gabc[2]=-bc[0]*bc[1];
  ggabc[0][0]=-2.*bc[1];
  ggabc[1][0]=bc[3]-bc[0]-bc[1];
  ggabc[2][0]=-bc[1];
  ggabc[1][1]=-2.*bc[0];
  ggabc[2][1]=-bc[0];
  ggabc[2][2]=0.;
  ggabc[0][1]=ggabc[1][0];
  ggabc[0][2]=ggabc[2][0];
  ggabc[1][2]=ggabc[2][1];
  if (ev[0]<ev[1]) {
    if (ev[1]<ev[2]) {
      for (int m=0;m+3<=order;m++) {
        for (int j=0;j<=m;j++) {
          Tensor<2,3> &ggv=gg[ijk++];
          ggv[0][0]=(ggabc[0][0]*lpv0[m-j]+2.*gabc[0]*lps0[m-j]
                    +abc*lpss0[m-j])*lpv1[j];
          ggv[1][0]=ggabc[1][0]*lpv0[m-j]*lpv1[j]+gabc[0]*lpv0[m-j]*lps1[j]
                   +gabc[1]*lps0[m-j]*lpv1[j]+abc*lps0[m-j]*lps1[j];
          ggv[2][0]=(ggabc[2][0]*lpv0[m-j]+gabc[2]*lps0[m-j])*lpv1[j];
          ggv[1][1]=(ggabc[1][1]*lpv1[j]+2.*gabc[1]*lps1[j]
                    +abc*lpss1[j])*lpv0[m-j];
          ggv[2][1]=lpv0[m-j]*(ggabc[2][1]*lpv1[j]+gabc[2]*lps1[j]);
          ggv[2][2]=ggabc[2][2]*lpv0[m-j]*lpv1[j];
          ggv[0][1]=ggv[1][0];
          ggv[0][2]=ggv[2][0];
          ggv[1][2]=ggv[2][1];
        }
      }
    } else if (ev[0]<ev[2]) { // never called:
      OBSOLETE("not programmed");
    } else { // never called:
      OBSOLETE("not programmed");
    }
  } else if (ev[0]<ev[2]) { // never called:
    OBSOLETE("not programmed");
  } else if (ev[1]<ev[2]) { // never called:
    OBSOLETE("not programmed");
  } else { // never called:
    OBSOLETE("not programmed");
  }
//face 3: v1, v2, v3
  abc=bc[0]*bc[1]*bc[2];
  gabc[0]=bc[1]*bc[2];
  gabc[1]=bc[0]*bc[2];
  gabc[2]=bc[0]*bc[1];
  ggabc[0][0]=0.;
  ggabc[1][0]=bc[2];
  ggabc[2][0]=bc[1];
  ggabc[1][1]=0.;
  ggabc[2][1]=bc[0];
  ggabc[2][2]=0.;
  ggabc[0][1]=ggabc[1][0];
  ggabc[0][2]=ggabc[2][0];
  ggabc[1][2]=ggabc[2][1];
  if (ev[1]<ev[2]) {
    if (ev[2]<ev[3]) {
      for (int m=0;m+3<=order;m++) {
        for (int j=0;j<=m;j++) {
          Tensor<2,3> &ggv=gg[ijk++];
          ggv[0][0]=ggabc[0][0]*lpv1[m-j]*lpv2[j];
          ggv[1][0]=(ggabc[1][0]*lpv1[m-j]+gabc[0]*lps1[m-j])*lpv2[j];
          ggv[2][0]=(ggabc[2][0]*lpv2[j]+gabc[0]*lps2[j])*lpv1[m-j];
          ggv[1][1]=(ggabc[1][1]*lpv1[m-j]+2.*gabc[1]*lps1[m-j]
                    +abc*lpss1[m-j])*lpv2[j];
          ggv[2][1]=ggabc[2][1]*lpv1[m-j]*lpv2[j]
                   +gabc[1]*lpv1[m-j]*lps2[j]
                   +gabc[2]*lps1[m-j]*lpv2[j]
                   +abc*lps1[m-j]*lps2[j];
          ggv[2][2]=(ggabc[2][2]*lpv2[j]+2.*gabc[2]*lps2[j]
                    +abc*lpss2[j])*lpv1[m-j];
          ggv[0][1]=ggv[1][0];
          ggv[0][2]=ggv[2][0];
          ggv[1][2]=ggv[2][1];
        }
      }
    } else if (ev[1]<ev[3]) {
      for (int m=0;m+3<=order;m++) {
        for (int j=0;j<=m;j++) {
          Tensor<2,3> &ggv=gg[ijk++];
          ggv[0][0]=ggabc[0][0]*lpv2[m-j]*lpv1[j];
          ggv[1][0]=(ggabc[1][0]*lpv1[j]+gabc[0]*lps1[j])*lpv2[m-j];
          ggv[2][0]=(ggabc[2][0]*lpv2[m-j]+gabc[0]*lps2[m-j])*lpv1[j];
          ggv[1][1]=(ggabc[1][1]*lpv1[j]+2.*gabc[1]*lps1[j]
                    +abc*lpss1[j])*lpv2[m-j];
          ggv[2][1]=ggabc[2][1]*lpv2[m-j]*lpv1[j]
                   +gabc[2]*lpv2[m-j]*lps1[j]
                   +gabc[1]*lps2[m-j]*lpv1[j]
                   +abc*lps2[m-j]*lps1[j];
          ggv[2][2]=(ggabc[2][2]*lpv2[m-j]+2.*gabc[2]*lps2[m-j]
                    +abc*lpss2[m-j])*lpv1[j];
          ggv[0][1]=ggv[1][0];
          ggv[0][2]=ggv[2][0];
          ggv[1][2]=ggv[2][1];
        }
      }
    } else { // never called:
      OBSOLETE("not programmed");
    }
  } else if (ev[1]<ev[3]) { // never called:
    OBSOLETE("not programmed");
  } else if (ev[2]<ev[3]) { // never called:
    OBSOLETE("not programmed");
  } else { // never called:
    OBSOLETE("not programmed");
  }

  if (order<4) return;
  REAL abcd=bc[0]*bc[1]*bc[2]*bc[3];
  Tensor<1,3> gabcd;
  gabcd[0]=bc[1]*bc[2]*(bc[3]-bc[0]);
  gabcd[1]=bc[2]*bc[0]*(bc[3]-bc[1]);
  gabcd[2]=bc[0]*bc[1]*(bc[3]-bc[2]);
  Tensor<2,3> ggabcd;
  ggabcd[0][0]=-2.*bc[1]*bc[2];
  ggabcd[1][0]=bc[2]*(bc[3]-bc[0]-bc[1]);
  ggabcd[2][0]=bc[1]*(bc[3]-bc[0]-bc[2]);
  ggabcd[1][1]=-2.*bc[2]*bc[0];
  ggabcd[2][1]=bc[0]*(bc[3]-bc[1]-bc[2]);
  ggabcd[2][2]=-2.*bc[0]*bc[1];
  ggabcd[0][1]=ggabcd[1][0];
  ggabcd[0][2]=ggabcd[2][0];
  ggabcd[1][2]=ggabcd[2][1];
  for (int m=0;m+4<=order;m++) {
    for (int k=0;k<=m;k++) {
      for (int j=0;j+k<=m;j++) {
        Tensor<2,3> &ggv=gg[ijk++];
        ggv[0][0]=(ggabcd[0][0]*lpv0[m-j-k]+2.*gabcd[0]*lps0[m-j-k]
                  +abcd*lpss0[m-j-k])*lpv1[j]*lpv2[k];
        ggv[1][0]=(ggabcd[1][0]*lpv0[m-j-k]*lpv1[j]
                  +gabcd[0]*lpv0[m-j-k]*lps1[j]
                  +gabcd[1]*lps0[m-j-k]*lpv1[j]
                  +abcd*lps0[m-j-k]*lps1[j])*lpv2[k];
        ggv[2][0]=(ggabcd[2][0]*lpv0[m-j-k]*lpv2[k]
                  +gabcd[0]*lpv0[m-j-k]*lps2[k]
                  +gabcd[2]*lps0[m-j-k]*lpv2[k]
                  +abcd*lps0[m-j-k]*lps2[k])*lpv1[j];
        ggv[1][1]=(ggabcd[1][1]*lpv1[j]
                  +2.*gabcd[1]*lps1[j]
                  +abcd*lpss1[j])*lpv0[m-j-k]*lpv2[k];
        ggv[2][1]=(ggabcd[2][1]*lpv1[j]*lpv2[k]
                  +gabcd[1]*lpv1[j]*lps2[k]
                  +gabcd[2]*lps1[j]*lpv2[k]
                  +abcd*lps1[j]*lps2[k])*lpv0[m-j-k];
        ggv[2][2]=(ggabcd[2][2]*lpv2[k]+2.*gabcd[2]*lps2[k]
                  +abcd*lpss2[k])*lpv0[m-j-k]*lpv1[j];
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalTetrahedronPolynomials::shapeFunctionsOrder(
const Shape *element,NumPtr<int> &renumber) const {
//triangular faces have shape functions that are not symmetric in
//  barycentric coordinates, so cannot be reordered
//since requiresShapeToComputeValues=true, no reordering is needed
  int n=getNumber();
  renumber.allocate(n);
  for (int i=0;i<n;i++) renumber[i]=i;
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalTetrahedronPolynomials::nonzeroShapesOnFace(
const Shape1 *e,int f,NumPtr<int> &indices) const {
  const Shape3 *e3=dynamic_cast<const Shape3*>(e);
  int nv=e->numberVertices();
  int nvof=e3->numberVerticesOnFace(f);
  int lpf=linesPerFace(nv,f);
  int nieslpof=numberInteriorEquallySpacedLatticePointsOnFace(nv,f,order);
  indices.allocate(nvof+lpf*(order-1)+nieslpof);
  int iv=0;
  for (int v=0;v<nvof;v++) indices[iv++]=e->elementVertexFromFace(f,v);
  for (int l=0;l<lpf;l++) {
    int offset=nv+elementLineFromFace(nv,f,l)*(order-1)-1;
    for (int i=1;i<order;i++) indices[iv++]=offset+i;
  }
  int offset=nv+e3->numberCodimension2Shapes()*(order-1);
  for (int i=0;i<f;i++) {
    offset+=numberInteriorEquallySpacedLatticePointsOnFace(nv,i,order);
  }
  for (int i=0;i<nieslpof;i++) indices[iv++]=offset+i;
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalTetrahedronPolynomials::printOn(ostream &os) const {
  os << "HierarchicalTetrahedronPolynomials: order = " << order << endl;
  ShapeFunction<3>::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalPrismPolynomials::computeValues(const Point<3> &p,
NumPtr<REAL> &v) const {
  CHECK_SAME(order01,order2);
  CHECK_TEST(order01<=1);
  v.allocate(6);
  Point<3> bc(p[0],p[1],1.-(p[0]+p[1]));
  NumPtr<REAL> lpv(2);
  lp->values(p[2],lpv);
  int iv=0;
//vertices
  for (int k=0;k<2;k++) {
    v[iv++]=bc[2]*lpv[k];
    v[iv++]=bc[0]*lpv[k];
    v[iv++]=bc[1]*lpv[k];
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalPrismPolynomials::computeGrads(const Point<3> &p,
NumPtr<Tensor<1,3> > &g) const {
  CHECK_SAME(order01,order2);
  CHECK_TEST(order01<=1);
  g.allocate(6);
  NumPtr<REAL> lpv(2);
  NumPtr<REAL> lps(2);
  Point<3> bc(p[0],p[1],1.-(p[0]+p[1]));

  lp->values(p[2],lpv);
  lp->slopes(p[2],lps);
  int iv=0;
//vertices
  for (int k=0;k<2;k++) {
    REAL lpvk=lpv[k];
    REAL lpsk=lps[k];
    {
      Tensor<1,3> &gv=g[iv++];
      gv[0]=-lpvk;
      gv[1]=-lpvk;
      gv[2]=bc[2]*lpsk;
    }
    {
      Tensor<1,3> &gv=g[iv++];
      gv[0]=lpvk;
      gv[1]=0.;
      gv[2]=bc[0]*lpsk;
    }
    {
      Tensor<1,3> &gv=g[iv++];
      gv[0]=0.;
      gv[1]=lpvk;
      gv[2]=bc[1]*lpsk;
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalPrismPolynomials::computeGradGrads(const Point<3> &p,
NumPtr<Tensor<2,3> > &gg) const {
  CHECK_SAME(order01,order2);
  CHECK_TEST(order01<=1);
  gg.allocate(6);
//vertices
  for (int v=0;v<6;v++) gg[v].initialize(0.);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalPrismPolynomials::computeValuesFor(const Shape *e,
const Point<SPACEDIM> &p,NumPtr<REAL> &v) const {
  CHECK_SAME(order01,order2);
  int ntp=tp->getNumber();
  int nlp=order2+1;
  v.allocate(ntp*nlp);
  Point<3> bc(p[0],p[1],1.-(p[0]+p[1]));
  NumPtr<int> ev(e->numberVertices());
  for (int i=0;i<e->numberVertices();i++) ev[i]=e->vertexIndex(i);

  NumPtr<REAL> lpv(nlp);
  lp->values(p[2],lpv);
  int iv=0;
//vertices
  for (int k=0;k<2;k++) {
    v[iv++]=bc[2]*lpv[k];
    v[iv++]=bc[0]*lpv[k];
    v[iv++]=bc[1]*lpv[k];
  }

  if (order2<2) return;
  HierarchicalSidePolynomial hsp;
  NumPtr<REAL> hspv0(order2-1);
  NumPtr<REAL> hspv1(order2-1);
  NumPtr<REAL> hspv2(order2-1);
  hsp.values(bc[0],hspv0);
  hsp.values(bc[1],hspv1);
  hsp.values(bc[2],hspv2);
  {
    NumPtr<REAL> *hspv=(ev[0]<ev[1] ? &hspv0 : &hspv2);
    REAL ab=bc[2]*bc[0]*lpv[0];
    for (int i=1;i<order2;i++) {
      v[iv++]=ab*(*hspv)[i-1];
    }
  }
  {
    NumPtr<REAL> *hspv=(ev[1]<ev[2] ? &hspv1 : &hspv0);
    REAL ab=bc[0]*bc[1]*lpv[0];
    for (int j=1;j<order2;j++) {
      v[iv++]=ab*(*hspv)[j-1];
    }
  }
  {
    NumPtr<REAL> *hspv=(ev[0]<ev[2] ? &hspv1 : &hspv2);
    REAL ab=bc[2]*bc[1]*lpv[0];
    for (int j=1;j<order2;j++) {
      v[iv++]=ab*(*hspv)[j-1];
    }
  }
  {
    NumPtr<REAL> *hspv=(ev[3]<ev[4] ? &hspv0 : &hspv2);
    REAL ab=bc[2]*bc[0]*lpv[1];
    for (int i=1;i<order2;i++) {
      v[iv++]=ab*(*hspv)[i-1];
    }
  }
  {
    NumPtr<REAL> *hspv=(ev[4]<ev[5] ? &hspv1 : &hspv0);
    REAL ab=bc[0]*bc[1]*lpv[1];
    for (int j=1;j<order2;j++) {
      v[iv++]=ab*(*hspv)[j-1];
    }
  }
  {
    NumPtr<REAL> *hspv=(ev[3]<ev[5] ? &hspv1 : &hspv2);
    REAL ab=bc[2]*bc[1]*lpv[1];
    for (int j=1;j<order2;j++) {
      v[iv++]=ab*(*hspv)[j-1];
    }
  }
  {
    if (ev[0]<ev[3]) {
      for (int k=2;k<=order2;k++) {
        v[iv++]=bc[2]*lpv[k];
      }
    } else {
      REAL s=bc[2];
      for (int k=2;k<=order2;k++) {
        v[iv++]=lpv[k]*s;
        s=-s;
      }
    }
  }
  {
    if (ev[1]<ev[4]) {
      for (int k=2;k<=order2;k++) {
        v[iv++]=bc[0]*lpv[k];
      }
    } else {
      REAL s=bc[0];
      for (int k=2;k<=order2;k++) {
        v[iv++]=lpv[k]*s;
        s=-s;
      }
    }
  }
  {
    if (ev[2]<ev[5]) {
      for (int k=2;k<=order2;k++) {
        v[iv++]=bc[1]*lpv[k];
      }
    } else {
      REAL s=bc[1];
      for (int k=2;k<=order2;k++) {
        v[iv++]=lpv[k]*s;
        s=-s;
      }
    }
  }

  if (order2>=3) {
    LegendrePolynomial lp;
    NumPtr<REAL> lpv0(order2-2);
    NumPtr<REAL> lpv1(order2-2);
    NumPtr<REAL> lpv2(order2-2);
    lp.values(bc[0],lpv0);
    lp.values(bc[1],lpv1);
    lp.values(bc[2],lpv2);
    {
//    face 0: v0, v1, v2
      REAL abc=bc[2]*bc[0]*bc[1]*lpv[0];
      NumPtr<REAL> *lpva=0;
      NumPtr<REAL> *lpvb=0;
      if (ev[0]<ev[1]) {
        if (ev[1]<ev[2]) {
          lpva=&lpv0;
          lpvb=&lpv1;
        } else if (ev[0]<ev[2]) {
          lpva=&lpv1;
          lpvb=&lpv0;
        } else {
          lpva=&lpv2;
          lpvb=&lpv0;
        }
      } else if (ev[0]<ev[2]) {
        lpva=&lpv2;
        lpvb=&lpv1;
      } else if (ev[1]<ev[2]) {
        lpva=&lpv1;
        lpvb=&lpv2;
      } else {
        lpva=&lpv0;
        lpvb=&lpv2;
      }
      for (int m=0;m+3<=order2;m++) {
        for (int i=0;i<=m;i++) {
          v[iv++]=abc*(*lpva)[m-i]*(*lpvb)[i];
        }
      }
    }
    {
//    face 1: v3, v4, v5
      REAL abc=bc[2]*bc[0]*bc[1]*lpv[1];
      NumPtr<REAL> *lpva=0;
      NumPtr<REAL> *lpvb=0;
      if (ev[3]<ev[4]) {
        if (ev[4]<ev[5]) {
          lpva=&lpv0;
          lpvb=&lpv1;
        } else if (ev[3]<ev[5]) {
          lpva=&lpv1;
          lpvb=&lpv0;
        } else {
          lpva=&lpv2;
          lpvb=&lpv0;
        }
      } else if (ev[3]<ev[5]) {
        lpva=&lpv2;
        lpvb=&lpv1;
      } else if (ev[4]<ev[5]) {
        lpva=&lpv1;
        lpvb=&lpv2;
      } else {
        lpva=&lpv0;
        lpvb=&lpv2;
      }
      for (int m=0;m+3<=order2;m++) {
        for (int i=0;i<=m;i++) {
          v[iv++]=abc*(*lpva)[m-i]*(*lpvb)[i];
        }
      }
    }
  }
  {
//  face 2: v0, v2, v3, v5
    int vmin=min(min(ev[0],ev[2]),min(ev[3],ev[5]));
    REAL ab=bc[2]*bc[1];
    if (ev[0]==vmin) {
      if (ev[2]<ev[3]) {
        for (int k=2;k<=order2;k++) {
          REAL lpvk=ab*lpv[k];
          for (int j=1;j<order2;j++) {
            v[iv++]=hspv1[j-1]*lpvk;
          }
        }
      } else {
        for (int j=1;j<order2;j++) {
          REAL hspvj=ab*hspv1[j-1];
          for (int k=2;k<=order2;k++) {
            v[iv++]=hspvj*lpv[k];
          }
        }
      }
    } else if (ev[2]==vmin) {
      if (ev[0]<ev[5]) {
        for (int k=2;k<=order2;k++) {
          REAL lpvk=ab*lpv[k];
          REAL s=1.;
          for (int j=1;j<order2;j++) {
            v[iv++]=hspv1[j-1]*s*lpvk;
            s=-s;
          }
        }
      } else {
        REAL s=1.;
        for (int j=1;j<order2;j++) {
          REAL hspvj=ab*hspv1[j-1]*s;
          for (int k=2;k<=order2;k++) {
            v[iv++]=hspvj*lpv[k];
          }
          s=-s;
        }
      }
    } else if (ev[3]==vmin) {
      if (ev[0]<ev[5]) {
        for (int j=1;j<order2;j++) {
          REAL hspvj=ab*hspv1[j-1];
          REAL s=1.;
          for (int k=2;k<=order2;k++) {
            v[iv++]=hspvj*lpv[k]*s;
            s=-s;
          }
        }
      } else {
        REAL s=1.;
        for (int k=2;k<=order2;k++) {
          REAL lpvk=ab*lpv[k]*s;
          for (int j=1;j<order2;j++) {
            v[iv++]=hspv1[j-1]*lpvk;
          }
          s=-s;
        }
      }
    } else {
      if (ev[2]<ev[3]) {
        REAL sj=1.;
        for (int j=1;j<order2;j++) {
          REAL hspvj=ab*hspv1[j-1]*sj;
          REAL sk=1.;
          for (int k=2;k<=order2;k++) {
            v[iv++]=hspvj*lpv[k]*sk;
            sk=-sk;
          }
          sj=-sj;
        }
      } else {
        REAL sk=1.;
        for (int k=2;k<=order2;k++) {
          REAL lpvk=ab*lpv[k]*sk;
          REAL sj=1.;
          for (int j=1;j<order2;j++) {
            v[iv++]=hspv1[j-1]*sj*lpvk;
            sj=-sj;
          }
          sk=-sk;
        }
      }
    }
  }
  {
//  face 3: v0, v1, v3, v4
    int vmin=min(min(ev[0],ev[1]),min(ev[3],ev[4]));
    REAL ab=bc[2]*bc[0];
    if (ev[0]==vmin) {
      if (ev[1]<ev[3]) {
        for (int k=2;k<=order2;k++) {
          REAL lpvk=ab*lpv[k];
          for (int i=1;i<order2;i++) {
            v[iv++]=hspv0[i-1]*lpvk;
          }
        }
      } else {
        for (int i=1;i<order2;i++) {
          REAL hspvi=ab*hspv0[i-1];
          for (int k=2;k<=order2;k++) {
            v[iv++]=hspvi*lpv[k];
          }
        }
      }
    } else if (ev[1]==vmin) {
      if (ev[0]<ev[4]) {
        for (int k=2;k<=order2;k++) {
          REAL lpvk=ab*lpv[k];
          REAL s=1.;
          for (int i=1;i<order2;i++) {
            v[iv++]=hspv0[i-1]*s*lpvk;
            s=-s;
          }
        }
      } else {
        REAL s=1.;
        for (int i=1;i<order2;i++) {
          REAL hspvi=ab*hspv0[i-1]*s;
          for (int k=2;k<=order2;k++) {
            v[iv++]=hspvi*lpv[k];
          }
          s=-s;
        }
      }
    } else if (ev[3]==vmin) {
      if (ev[0]<ev[4]) {
        for (int i=1;i<order2;i++) {
          REAL hspvi=ab*hspv0[i-1];
          REAL s=1.;
          for (int k=2;k<=order2;k++) {
            v[iv++]=hspvi*lpv[k]*s;
            s=-s;
          }
        }
      } else {
        REAL s=1.;
        for (int k=2;k<=order2;k++) {
          REAL lpvk=ab*lpv[k]*s;
          for (int i=1;i<order2;i++) {
            v[iv++]=hspv0[i-1]*lpvk;
          }
          s=-s;
        }
      }
    } else {
      if (ev[1]<ev[3]) {
        REAL si=1.;
        for (int i=1;i<order2;i++) {
          REAL hspvi=ab*hspv0[i-1]*si;
          REAL sk=1.;
          for (int k=2;k<=order2;k++) {
            v[iv++]=hspvi*lpv[k]*sk;
            sk=-1;
          }
          si=-si;
        }
      } else {
        REAL sk=1.;
        for (int k=2;k<=order2;k++) {
          REAL lpvk=ab*lpv[k]*sk;
          REAL si=1.;
          for (int i=1;i<order2;i++) {
            v[iv++]=hspv0[i-1]*si*lpvk;
            si=-si;
          }
          sk=-sk;
        }
      }
    }
  }
  {
//  face 4: v1, v2, v4, v5
    int vmin=min(min(ev[1],ev[2]),min(ev[4],ev[5]));
    REAL ab=bc[0]*bc[1];
    if (ev[1]==vmin) {
      if (ev[2]<ev[4]) {
        for (int k=2;k<=order2;k++) {
          REAL lpvk=ab*lpv[k];
          for (int j=1;j<order2;j++) {
            v[iv++]=hspv1[j-1]*lpvk;
          }
        }
      } else {
        for (int j=1;j<order2;j++) {
          REAL hspvj=ab*hspv1[j-1];
          for (int k=2;k<=order2;k++) {
            v[iv++]=hspvj*lpv[k];
          }
        }
      }
    } else if (ev[2]==vmin) {
      if (ev[1]<ev[5]) {
        for (int k=2;k<=order2;k++) {
          REAL lpvk=ab*lpv[k];
          REAL s=1.;
          for (int j=1;j<order2;j++) {
            v[iv++]=hspv1[j-1]*s*lpvk;
            s=-s;
          }
        }
      } else {
        REAL s=1.;
        for (int j=1;j<order2;j++) {
          REAL hspvj=ab*hspv1[j-1]*s;
          for (int k=2;k<=order2;k++) {
            v[iv++]=hspvj*lpv[k];
          }
          s=-s;
        }
      }
    } else if (ev[4]==vmin) {
      if (ev[1]<ev[5]) {
        for (int j=1;j<order2;j++) {
          REAL hspvj=ab*hspv1[j-1];
          REAL s=1.;
          for (int k=2;k<=order2;k++) {
            v[iv++]=hspvj*lpv[k]*s;
            s=-s;
          }
        }
      } else {
        REAL s=1.;
        for (int k=2;k<=order2;k++) {
          REAL lpvk=ab*lpv[k]*s;
          for (int j=1;j<order2;j++) {
            v[iv++]=hspv1[j-1]*lpvk;
          }
          s=-s;
        }
      }
    } else {
      if (ev[2]<ev[4]) {
        REAL sj=1.;
        for (int j=1;j<order2;j++) {
          REAL hspvj=ab*hspv1[j-1]*sj;
          REAL sk=1.;
          for (int k=2;k<=order2;k++) {
            v[iv++]=hspvj*lpv[k]*sk;
            sk=-sk;
          }
          sj=-sj;
        }
      } else {
        REAL sk=1.;
        for (int k=2;k<=order2;k++) {
          REAL lpvk=ab*lpv[k]*sk;
          REAL sj=1.;
          for (int j=1;j<order2;j++) {
            v[iv++]=hspv1[j-1]*sj*lpvk;
            sj=-sj;
          }
          sk=-sk;
        }
      }
    }
  }

  if (order2<3) return;
  {
    REAL abc=bc[0]*bc[1]*bc[2];
    for (int k=2;k<=order2;k++) {
      REAL lpvk=abc*lpv[k];
      for (int m=0;m+3<=order2;m++) {
        for (int j=0;j<=m;j++) {
          v[iv++]=hspv0[m-j]*hspv1[j]*lpvk;
        }
      }
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalPrismPolynomials::computeGradsFor(const Shape *e,
const Point<SPACEDIM> &p,NumPtr<Tensor<1,SPACEDIM> > &g) const {
  CHECK_SAME(order01,order2);
  int ntp=tp->getNumber();
  int nlp=order2+1;
  g.allocate(ntp*nlp);
  NumPtr<REAL> lpv(nlp);
  NumPtr<REAL> lps(nlp);
  Point<3> bc(p[0],p[1],1.-(p[0]+p[1]));
  NumPtr<int> ev(e->numberVertices());
  for (int i=0;i<e->numberVertices();i++) ev[i]=e->vertexIndex(i);

  lp->values(p[2],lpv);
  lp->slopes(p[2],lps);
  int iv=0;
//vertices
  for (int k=0;k<2;k++) {
    REAL lpvk=lpv[k];
    REAL lpsk=lps[k];
    {
      Tensor<1,3> &gv=g[iv++];
      gv[0]=-lpvk;
      gv[1]=-lpvk;
      gv[2]=bc[2]*lpsk;
    }
    {
      Tensor<1,3> &gv=g[iv++];
      gv[0]=lpvk;
      gv[1]=0.;
      gv[2]=bc[0]*lpsk;
    }
    {
      Tensor<1,3> &gv=g[iv++];
      gv[0]=0.;
      gv[1]=lpvk;
      gv[2]=bc[1]*lpsk;
    }
  }
  if (order2<2) return;
  HierarchicalSidePolynomial hsp;
  NumPtr<REAL> hspv0(order2-1);
  NumPtr<REAL> hspv1(order2-1); 
  NumPtr<REAL> hspv2(order2-1);
  hsp.values(bc[0],hspv0);
  hsp.values(bc[1],hspv1);
  hsp.values(bc[2],hspv2);
  NumPtr<REAL> hsps0(order2-1);
  NumPtr<REAL> hsps1(order2-1); 
  NumPtr<REAL> hsps2(order2-1);
  hsp.slopes(bc[0],hsps0);
  hsp.slopes(bc[1],hsps1);
  hsp.slopes(bc[2],hsps2);
//line 0
  REAL lpvk=lpv[0];
  REAL lpsk=lps[0];
  REAL ab=bc[2]*bc[0];
  Tensor<1,2> gab;
  gab[0]=bc[2]-bc[0];
  gab[1]=-bc[0];
  if (ev[0]<ev[1]) {
    for (int i=1;i<order2;i++) {
      Tensor<1,3> &gv=g[iv++];
      gv[0]=(gab[0]*hspv0[i-1]+ab*hsps0[i-1])*lpvk;
      gv[1]=gab[1]*hspv0[i-1]*lpvk;
      gv[2]=ab*hspv0[i-1]*lpsk;
    }
  } else { // not used
    for (int i=1;i<order2;i++) {
      Tensor<1,3> &gv=g[iv++];
      gv[0]=(gab[0]*hspv2[i-1]-ab*hsps2[i-1])*lpvk;
      gv[1]=(gab[1]*hspv2[i-1]-ab*hsps2[i-1])*lpvk;
      gv[2]=ab*hspv2[i-1]*lpsk;
    }
  }
//line 1
  ab=bc[0]*bc[1];
  gab[0]=bc[1];
  gab[1]=bc[0];
  if (ev[1]<ev[2]) {
    for (int j=1;j<order2;j++) {
      Tensor<1,3> &gv=g[iv++];
      gv[0]=gab[0]*hspv1[j-1]*lpvk;
      gv[1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*lpvk;
      gv[2]=ab*hspv1[j-1]*lpsk;
    }
  } else {
    for (int j=1;j<order2;j++) {
      Tensor<1,3> &gv=g[iv++];
      gv[0]=(gab[0]*hspv0[j-1]+ab*hsps0[j-1])*lpvk;
      gv[1]=gab[1]*hspv0[j-1]*lpvk;
      gv[2]=ab*hspv0[j-1]*lpsk;
    }
  }
//line 2
  ab=bc[2]*bc[1];
  gab[0]=-bc[1];
  gab[1]=bc[2]-bc[1];
  if (ev[0]<ev[2]) {
    for (int j=1;j<order2;j++) {
      Tensor<1,3> &gv=g[iv++];
      gv[0]=gab[0]*hspv1[j-1]*lpvk;
      gv[1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*lpvk;
      gv[2]=ab*hspv1[j-1]*lpsk;
    }
  } else { // not used
    for (int j=1;j<order2;j++) {
      Tensor<1,3> &gv=g[iv++];
      gv[0]=(gab[0]*hspv2[j-1]-ab*hsps2[j-1])*lpvk;
      gv[1]=(gab[1]*hspv2[j-1]-ab*hsps2[j-1])*lpvk;
      gv[2]=ab*hspv2[j-1]*lpsk;
    }
  }
//line 3
  lpvk=lpv[1];
  lpsk=lps[1];
  ab=bc[2]*bc[0];
  gab[0]=bc[2]-bc[0];
  gab[1]=-bc[0];
  if (ev[3]<ev[4]) {
    for (int i=1;i<order2;i++) {
      Tensor<1,3> &gv=g[iv++];
      gv[0]=(gab[0]*hspv0[i-1]+ab*hsps0[i-1])*lpvk;
      gv[1]=gab[1]*hspv0[i-1]*lpvk;
      gv[2]=ab*hspv0[i-1]*lpsk;
    }
  } else {
    for (int i=1;i<order2;i++) {
      Tensor<1,3> &gv=g[iv++];
      gv[0]=(gab[0]*hspv2[i-1]-ab*hsps2[i-1])*lpvk;
      gv[1]=(gab[1]*hspv2[i-1]-ab*hsps2[i-1])*lpvk;
      gv[2]=ab*hspv2[i-1]*lpsk;
    }
  }
//line 4
  ab=bc[0]*bc[1];
  gab[0]=bc[1];
  gab[1]=bc[0];
  if (ev[4]<ev[5]) {
    for (int j=1;j<order2;j++) {
      Tensor<1,3> &gv=g[iv++];
      gv[0]=gab[0]*hspv1[j-1]*lpvk;
      gv[1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*lpvk;
      gv[2]=ab*hspv1[j-1]*lpsk;
    }
  } else {
    for (int j=1;j<order2;j++) {
      Tensor<1,3> &gv=g[iv++];
      gv[0]=(gab[0]*hspv0[j-1]+ab*hsps0[j-1])*lpvk;
      gv[1]=gab[1]*hspv0[j-1]*lpvk;
      gv[2]=ab*hspv0[j-1]*lpsk;
    }
  }
//line 5
  ab=bc[2]*bc[1];
  gab[0]=-bc[1];
  gab[1]=bc[2]-bc[1];
  if (ev[3]<ev[5]) {
    for (int j=1;j<order2;j++) {
      Tensor<1,3> &gv=g[iv++];
      gv[0]=gab[0]*hspv1[j-1]*lpvk;
      gv[1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*lpvk;
      gv[2]=ab*hspv1[j-1]*lpsk;
    }
  } else {
    for (int j=1;j<order2;j++) {
      Tensor<1,3> &gv=g[iv++];
      gv[0]=(gab[0]*hspv2[j-1]-ab*hsps2[j-1])*lpvk;
      gv[1]=(gab[1]*hspv2[j-1]-ab*hsps2[j-1])*lpvk;
      gv[2]=ab*hspv2[j-1]*lpsk;
    }
  }
//line 6
  if (ev[0]<ev[3]) {
    for (int k=2;k<=order2;k++) {
      Tensor<1,3> &gv=g[iv++];
      gv[0]=-lpv[k];
      gv[1]=gv[0];
      gv[2]=bc[2]*lps[k];
    }
  } else { // not used
    REAL s=1.;
    for (int k=2;k<=order2;k++) {
      Tensor<1,3> &gv=g[iv++];
      gv[0]=-lpv[k]*s;
      gv[1]=gv[0];
      gv[2]=bc[2]*lps[k]*s;
      s=-s;
    }
  }
//line 7
  if (ev[1]<ev[4]) {
    for (int k=2;k<=order2;k++) {
      Tensor<1,3> &gv=g[iv++];
      gv[0]=lpv[k];
      gv[1]=0.;
      gv[2]=bc[0]*lps[k];
    }
  } else {
    REAL s=1.;
    for (int k=2;k<=order2;k++) {
      Tensor<1,3> &gv=g[iv++];
      gv[0]=lpv[k]*s;
      gv[1]=0.;
      gv[2]=bc[0]*lps[k]*s;
      s=-s;
    }
  }
//line 8
  if (ev[2]<ev[5]) {
    for (int k=2;k<=order2;k++) {
      Tensor<1,3> &gv=g[iv++];
      gv[0]=0.;
      gv[1]=lpv[k];
      gv[2]=bc[1]*lps[k];
    }
  } else {
    REAL s=1.;
    for (int k=2;k<=order2;k++) {
      Tensor<1,3> &gv=g[iv++];
      gv[0]=0.;
      gv[1]=lpv[k]*s;
      gv[2]=bc[1]*lps[k]*s;
      s=-s;
    }
  }

  if (order2>=3) {
    LegendrePolynomial lp;
    NumPtr<REAL> lpv0(order2-2);
    NumPtr<REAL> lpv1(order2-2);
    NumPtr<REAL> lpv2(order2-2);
    lp.values(bc[0],lpv0);
    lp.values(bc[1],lpv1);
    lp.values(bc[2],lpv2);
    NumPtr<REAL> lps0(order2-2);
    NumPtr<REAL> lps1(order2-2);
    NumPtr<REAL> lps2(order2-2);
    lp.slopes(bc[0],lps0);
    lp.slopes(bc[1],lps1);
    lp.slopes(bc[2],lps2);
    REAL abc=bc[2]*bc[0]*bc[1];
    Tensor<1,2> gabc;
    gabc[0]=(bc[2]-bc[0])*bc[1];
    gabc[1]=(bc[2]-bc[1])*bc[0];
//  face 0
    lpvk=lpv[0];
    lpsk=lps[0];
    if (ev[0]<ev[1]) {
      if (ev[1]<ev[2]) {
        for (int m=0;m+3<=order2;m++) {
          for (int i=0;i<=m;i++) {
#ifdef DEBUG
//          cout << "\tiv = " << iv << endl;
#endif
            Tensor<1,3> &gv=g[iv++];
            gv[0]=(gabc[0]*lpv0[m-i]+abc*lps0[m-i])*lpv1[i]*lpvk;
            gv[1]=(gabc[1]*lpv1[i]+abc*lps1[i])*lpv0[m-i]*lpvk;
            gv[2]=abc*lpv0[m-i]*lpv1[i]*lpsk;
          }
        }
      } else if (ev[0]<ev[2]) {
        for (int m=0;m+3<=order2;m++) {
          for (int i=0;i<=m;i++) {
            Tensor<1,3> &gv=g[iv++];
            gv[0]=(gabc[0]*lpv0[i]+abc*lps0[i])*lpv1[m-i]*lpvk;
            gv[1]=(gabc[1]*lpv1[m-i]+abc*lps1[m-i])*lpv0[i]*lpvk;
            gv[2]=abc*lpv1[m-i]*lpv0[i]*lpsk;
          }
        }
      } else { // not used
        for (int m=0;m+3<=order2;m++) {
          for (int i=0;i<=m;i++) {
            Tensor<1,3> &gv=g[iv++];
            gv[0]=(gabc[0]*lpv2[m-i]*lpv0[i]
                  +abc*(-lps2[m-i]*lpv0[i]+lpv2[m-i]*lps0[i]))*lpvk;
            gv[1]=(gabc[1]*lpv2[m-i]-abc*lps2[m-i])*lpv0[i]*lpvk;
            gv[2]=abc*lpv2[m-i]*lpv0[i]*lpsk;
          }
        }
      }
    } else if (ev[0]<ev[2]) { // not used
      for (int m=0;m+3<=order2;m++) {
        for (int i=0;i<=m;i++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=(gabc[0]*lpv2[m-i]-abc*lps2[m-i])*lpv1[i]*lpvk;
          gv[1]=(gabc[1]*lpv2[m-i]*lpv1[i]
                +abc*(-lps2[m-i]*lpv1[i]+lpv2[m-i]*lps1[i]))*lpvk;
          gv[2]=abc*lpv2[m-i]*lpv1[i]*lpsk;
        }
      }
    } else if (ev[1]<ev[2]) { // not used
      for (int m=0;m+3<=order2;m++) {
        for (int i=0;i<=m;i++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=(gabc[0]*lpv2[i]-abc*lps2[i])*lpv1[m-i]*lpvk;
          gv[1]=(gabc[1]*lpv2[i]*lpv1[m-i]
                +abc*(lps1[m-i]*lpv2[i]-lpv1[m-i]*lps2[i]))*lpvk;
          gv[2]=abc*lpv1[m-i]*lpv2[i]*lpsk;
        }
      }
    } else { // not used
      for (int m=0;m+3<=order2;m++) {
        for (int i=0;i<=m;i++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=(gabc[0]*lpv2[i]*lpv0[m-i]
                +abc*(lps0[m-i]*lpv2[i]-lpv0[m-i]*lps2[i]))*lpvk;
          gv[1]=(gabc[1]*lpv2[i]-abc*lps2[i])*lpv0[m-i]*lpvk;
          gv[2]=abc*lpv0[m-i]*lpv2[i]*lpsk;
        }
      }
    }
//  face 1
    lpvk=lpv[1];
    lpsk=lps[1];
    if (ev[3]<ev[4]) {
      if (ev[4]<ev[5]) {
        for (int m=0;m+3<=order2;m++) {
          for (int i=0;i<=m;i++) {
            Tensor<1,3> &gv=g[iv++];
            gv[0]=(gabc[0]*lpv0[m-i]+abc*lps0[m-i])*lpv1[i]*lpvk;
            gv[1]=(gabc[1]*lpv1[i]+abc*lps1[i])*lpv0[m-i]*lpvk;
            gv[2]=abc*lpv0[m-i]*lpv1[i]*lpsk;
          }
        }
      } else if (ev[3]<ev[5]) {
        for (int m=0;m+3<=order2;m++) {
          for (int i=0;i<=m;i++) {
            Tensor<1,3> &gv=g[iv++];
            gv[0]=(gabc[0]*lpv0[i]+abc*lps0[i])*lpv1[m-i]*lpvk;
            gv[1]=(gabc[1]*lpv1[m-i]+abc*lps1[m-i])*lpv0[i]*lpvk;
            gv[2]=abc*lpv1[m-i]*lpv0[i]*lpsk;
          }
        }
      } else {
        for (int m=0;m+3<=order2;m++) {
          for (int i=0;i<=m;i++) {
            Tensor<1,3> &gv=g[iv++];
            gv[0]=(gabc[0]*lpv2[m-i]*lpv0[i]
                  +abc*(-lps2[m-i]*lpv0[i]+lpv2[m-i]*lps0[i]))*lpvk;
            gv[1]=(gabc[1]*lpv2[m-i]-abc*lps2[m-i])*lpv0[i]*lpvk;
            gv[2]=abc*lpv2[m-i]*lpv0[i]*lpsk;
          }
        }
      }
    } else if (ev[3]<ev[5]) {
      for (int m=0;m+3<=order2;m++) {
        for (int i=0;i<=m;i++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=(gabc[0]*lpv2[m-i]-abc*lps2[m-i])*lpv1[i]*lpvk;
          gv[1]=(gabc[1]*lpv2[m-i]*lpv1[i]
                +abc*(-lps2[m-i]*lpv1[i]+lpv2[m-i]*lps1[i]))*lpvk;
          gv[2]=abc*lpv2[m-i]*lpv1[i]*lpsk;
        }
      }
    } else if (ev[4]<ev[5]) {
      for (int m=0;m+3<=order2;m++) {
        for (int i=0;i<=m;i++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=(gabc[0]*lpv2[i]-abc*lps2[i])*lpv1[m-i]*lpvk;
          gv[1]=(gabc[1]*lpv2[i]*lpv1[m-i]
                +abc*(lps1[m-i]*lpv2[i]-lpv1[m-i]*lps2[i]))*lpvk;
          gv[2]=abc*lpv1[m-i]*lpv2[i]*lpsk;
        }
      }
    } else {
      for (int m=0;m+3<=order2;m++) {
        for (int i=0;i<=m;i++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=(gabc[0]*lpv2[i]*lpv0[m-i]
                +abc*(lps0[m-i]*lpv2[i]-lpv0[m-i]*lps2[i]))*lpvk;
          gv[1]=(gabc[1]*lpv2[i]-abc*lps2[i])*lpv0[m-i]*lpvk;
          gv[2]=abc*lpv0[m-i]*lpv2[i]*lpsk;
        }
      }
    }
  }
//face 2
  int vmin=min(min(ev[0],ev[2]),min(ev[3],ev[5]));
  ab=bc[2]*bc[1];
  gab[0]=-bc[1];
  gab[1]=bc[2]-bc[1];
  if (ev[0]==vmin) {
    if (ev[2]<ev[3]) {
      for (int k=2;k<=order2;k++) {
        for (int j=1;j<order2;j++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=gab[0]*hspv1[j-1]*lpv[k];
          gv[1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*lpv[k];
          gv[2]=ab*hspv1[j-1]*lps[k];
        }
      }
    } else {
      for (int j=1;j<order2;j++) {
        for (int k=2;k<=order2;k++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=gab[0]*hspv1[j-1]*lpv[k];
          gv[1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*lpv[k];
          gv[2]=ab*hspv1[j-1]*lps[k];
        }
      }
    }
  } else if (ev[2]==vmin) {
    if (ev[0]<ev[5]) { // not used
      for (int k=2;k<=order2;k++) {
        REAL s=1.;
        for (int j=1;j<order2;j++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=gab[0]*hspv1[j-1]*s*lpv[k];
          gv[1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*s*lpv[k];
          gv[2]=ab*hspv1[j-1]*s*lps[k];
          s=-s;
        }
      }
    } else { // not used
      REAL s=1.;
      for (int j=1;j<order2;j++) {
        for (int k=2;k<=order2;k++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=gab[0]*hspv1[j-1]*s*lpv[k];
          gv[1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*s*lpv[k];
          gv[2]=ab*hspv1[j-1]*s*lps[k];
        }
        s=-s;
      }
    }
  } else if (ev[3]==vmin) {
    if (ev[0]<ev[5]) { // not used
      for (int j=1;j<order2;j++) {
        REAL s=1.;
        for (int k=2;k<=order2;k++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=gab[0]*hspv1[j-1]*lpv[k]*s;
          gv[1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*lpv[k]*s;
          gv[2]=ab*hspv1[j-1]*lps[k]*s;
          s=-s;
        }
      }
    } else { // not used
      REAL s=1.;
      for (int k=2;k<=order2;k++) {
        for (int j=1;j<order2;j++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=gab[0]*hspv1[j-1]*lpv[k]*s;
          gv[1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*lpv[k]*s;
          gv[2]=ab*hspv1[j-1]*lps[k]*s;
        }
        s=-s;
      }
    }
  } else {
    if (ev[2]<ev[3]) { // not used
      REAL sj=1.;
      for (int j=1;j<order2;j++) {
        REAL sk=1.;
        for (int k=2;k<=order2;k++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=gab[0]*hspv1[j-1]*sj*lpv[k]*sk;
          gv[1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*sj*lpv[k]*sk;
          gv[2]=ab*hspv1[j-1]*sj*lps[k]*sk;
          sk=-sk;
        }
        sj=-sj;
      }
    } else { // not used
      REAL sk=1.;
      for (int k=2;k<=order2;k++) {
        REAL sj=1.;
        for (int j=1;j<order2;j++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=gab[0]*hspv1[j-1]*sj*lpv[k]*sk;
          gv[1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*sj*lpv[k]*sk;
          gv[2]=ab*hspv1[j-1]*sj*lps[k]*sk;
          sj=-sj;
        }
        sk=-sk;
      }
    }
  }
//face 3
  vmin=min(min(ev[0],ev[1]),min(ev[3],ev[4]));
  ab=bc[2]*bc[0];
  gab[0]=bc[2]-bc[0];
  gab[1]=-bc[0];
  if (ev[0]==vmin) {
    if (ev[1]<ev[3]) {
      for (int k=2;k<=order2;k++) {
        for (int i=1;i<order2;i++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=(gab[0]*hspv0[i-1]+ab*hsps0[i-1])*lpv[k];
          gv[1]=gab[1]*hspv0[i-1]*lpv[k];
          gv[2]=ab*hspv0[i-1]*lps[k];
        }
      }
    } else {
      for (int i=1;i<order2;i++) {
        for (int k=2;k<=order2;k++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=(gab[0]*hspv0[i-1]+ab*hsps0[i-1])*lpv[k];
          gv[1]=gab[1]*hspv0[i-1]*lpv[k];
          gv[2]=ab*hspv0[i-1]*lps[k];
        }
      }
    }
  } else if (ev[1]==vmin) { // not used
    if (ev[0]<ev[4]) {
      for (int k=2;k<=order2;k++) {
        REAL s=1.;
        for (int i=1;i<order2;i++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=(gab[0]*hspv0[i-1]+ab*hsps0[i-1])*s*lpv[k];
          gv[1]=gab[1]*hspv0[i-1]*s*lpv[k];
          gv[2]=ab*hspv0[i-1]*s*lps[k];
          s=-s;
        }
      }
    } else { // not used
      REAL s=1.;
      for (int i=1;i<order2;i++) {
        for (int k=2;k<=order2;k++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=(gab[0]*hspv0[i-1]+ab*hsps0[i-1])*s*lpv[k];
          gv[1]=gab[1]*hspv0[i-1]*s*lpv[k];
          gv[2]=ab*hspv0[i-1]*s*lps[k];
        }
        s=-s;
      }
    }
  } else if (ev[3]==vmin) { // not used
    if (ev[0]<ev[4]) {
      for (int i=1;i<order2;i++) {
        REAL s=1.;
        for (int k=2;k<=order2;k++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=(gab[0]*hspv0[i-1]+ab*hsps0[i-1])*lpv[k]*s;
          gv[1]=gab[1]*hspv0[i-1]*lpv[k]*s;
          gv[2]=ab*hspv0[i-1]*lps[k]*s;
          s=-s;
        }
      }
    } else { // not used
      REAL s=1.;
      for (int k=2;k<=order2;k++) {
        for (int i=1;i<order2;i++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=(gab[0]*hspv0[i-1]+ab*hsps0[i-1])*lpv[k]*s;
          gv[1]=gab[1]*hspv0[i-1]*lpv[k]*s;
          gv[2]=ab*hspv0[i-1]*lps[k]*s;
        }
        s=-s;
      }
    }
  } else {
    if (ev[1]<ev[3]) { // not used
      REAL si=1.;
      for (int i=1;i<order2;i++) {
        REAL sk=1.;
        for (int k=2;k<=order2;k++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=(gab[0]*hspv0[i-1]+ab*hsps0[i-1])*si*lpv[k]*sk;
          gv[1]=gab[1]*hspv0[i-1]*si*lpv[k]*sk;
          gv[2]=ab*hspv0[i-1]*si*lps[k]*sk;
          sk=-sk;
        }
        si=-si;
      }
    } else { // not used
      REAL sk=1.;
      for (int k=2;k<=order2;k++) {
        REAL si=1.;
        for (int i=1;i<order2;i++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=(gab[0]*hspv0[i-1]+ab*hsps0[i-1])*si*lpv[k]*sk;
          gv[1]=gab[1]*hspv0[i-1]*si*lpv[k]*sk;
          gv[2]=ab*hspv0[i-1]*si*lps[k]*sk;
          si=-si;
        }
        sk=-sk;
      }
    }
  }
//face 4
  vmin=min(min(ev[1],ev[2]),min(ev[4],ev[5]));
  ab=bc[0]*bc[1];
  gab[0]=bc[1];
  gab[1]=bc[0];
  if (ev[1]==vmin) {
    if (ev[2]<ev[4]) {
      for (int k=2;k<=order2;k++) {
        for (int j=1;j<order2;j++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=gab[0]*hspv1[j-1]*lpv[k];
          gv[1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*lpv[k];
          gv[2]=ab*hspv1[j-1]*lps[k];
        }
      }
    } else {
      for (int j=1;j<order2;j++) {
        for (int k=2;k<=order2;k++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=gab[0]*hspv1[j-1]*lpv[k];
          gv[1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*lpv[k];
          gv[2]=ab*hspv1[j-1]*lps[k];
        }
      }
    }
  } else if (ev[2]==vmin) {
    if (ev[1]<ev[5]) {
      for (int k=2;k<=order2;k++) {
        REAL s=1.;
        for (int j=1;j<order2;j++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=gab[0]*hspv1[j-1]*s*lpv[k];
          gv[1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*s*lpv[k];
          gv[2]=ab*hspv1[j-1]*s*lps[k];
          s=-s;
        }
      }
    } else {
      REAL s=1.;
      for (int j=1;j<order2;j++) {
        for (int k=2;k<=order2;k++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=gab[0]*hspv1[j-1]*s*lpv[k];
          gv[1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*s*lpv[k];
          gv[2]=ab*hspv1[j-1]*s*lps[k];
        }
        s=-s;
      }
    }
  } else if (ev[4]==vmin) {
    if (ev[1]<ev[5]) {
      for (int j=1;j<order2;j++) {
        REAL s=1.;
        for (int k=2;k<=order2;k++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=gab[0]*hspv1[j-1]*lpv[k]*s;
          gv[1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*lpv[k]*s;
          gv[2]=ab*hspv1[j-1]*lps[k]*s;
          s=-s;
        }
      }
    } else {
      REAL s=1.;
      for (int k=2;k<=order2;k++) {
        for (int j=1;j<order2;j++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=gab[0]*hspv1[j-1]*lpv[k]*s;
          gv[1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*lpv[k]*s;
          gv[2]=ab*hspv1[j-1]*lps[k]*s;
        }
        s=-s;
      }
    }
  } else {
    if (ev[2]<ev[4]) {
      REAL sj=1.;
      for (int j=1;j<order2;j++) {
        REAL sk=1.;
        for (int k=2;k<=order2;k++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=gab[0]*hspv1[j-1]*sj*lpv[k]*sk;
          gv[1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*sj*lpv[k]*sk;
          gv[2]=ab*hspv1[j-1]*sj*lps[k]*sk;
          sk=-sk;
        }
        sj=-sj;
      }
    } else {
      REAL sk=1.;
      for (int k=2;k<=order2;k++) {
        REAL sj=1.;
        for (int j=1;j<order2;j++) {
          Tensor<1,3> &gv=g[iv++];
          gv[0]=gab[0]*hspv1[j-1]*sj*lpv[k]*sk;
          gv[1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*sj*lpv[k]*sk;
          gv[2]=ab*hspv1[j-1]*sj*lps[k]*sk;
          sj=-sj;
        }
        sk=-sk;
      }
    }
  }

  if (order2<3) return;
  REAL abc=bc[0]*bc[1]*bc[2];
  Tensor<1,2> gabc;
  gabc[0]=(bc[2]-bc[0])*bc[1];
  gabc[1]=(bc[2]-bc[1])*bc[0];
//interior
  for (int k=2;k<=order2;k++) {
    for (int m=0;m+3<=order2;m++) {
      for (int j=0;j<=m;j++) {
        Tensor<1,3> &gv=g[iv++];
        gv[0]=(gabc[0]*hspv0[m-j]+abc*hsps0[m-j])*hspv1[j]*lpv[k];
        gv[1]=(gabc[1]*hspv1[j]+abc*hsps1[j])*hspv0[m-j]*lpv[k];
        gv[2]=abc*hspv0[m-j]*hspv1[j]*lps[k];
      }
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalPrismPolynomials::computeGradGradsFor(const Shape *e,
const Point<SPACEDIM> &p,NumPtr<Tensor<2,SPACEDIM> > &gg) const {
  CHECK_SAME(order01,order2);
  int ntp=tp->getNumber();
  int nlp=order2+1;
  gg.allocate(ntp*nlp);
  NumPtr<REAL> lpv(nlp);
  NumPtr<REAL> lps(nlp);
  NumPtr<REAL> lpss(nlp);
  Point<3> bc(p[0],p[1],1.-(p[0]+p[1]));
  NumPtr<int> ev(e->numberVertices());
  for (int i=0;i<e->numberVertices();i++) ev[i]=e->vertexIndex(i);

  lp->values(p[2],lpv);
  lp->slopes(p[2],lps);
  lp->slope2s(p[2],lpss);
  int iv=0;
//vertices
  for (int k=0;k<2;k++) {
    REAL lpvk=lpv[k];
    REAL lpsk=lps[k];
    REAL lpssk=lpss[k];
    {
      Tensor<2,3> &ggv=gg[iv++];
      ggv[0][0]=0.;
      ggv[1][0]=0.;
      ggv[2][0]=-lpsk;
      ggv[1][1]=0.;
      ggv[2][1]=-lpsk;
      ggv[2][2]=bc[2]*lpssk;
      ggv[0][1]=ggv[1][0];
      ggv[0][2]=ggv[2][0];
      ggv[1][2]=ggv[2][1];
    }
    {
      Tensor<2,3> &ggv=gg[iv++];
      ggv[0][0]=0.;
      ggv[1][0]=0.;
      ggv[2][0]=lpsk;
      ggv[1][1]=0.;
      ggv[2][1]=0.;
      ggv[2][2]=bc[0]*lpssk;
      ggv[0][1]=ggv[1][0];
      ggv[0][2]=ggv[2][0];
      ggv[1][2]=ggv[2][1];
    }
    {
      Tensor<2,3> &ggv=gg[iv++];
      ggv[0][0]=0.;
      ggv[1][0]=0.;
      ggv[2][0]=0.;
      ggv[1][1]=0.;
      ggv[2][1]=lpsk;
      ggv[2][2]=bc[1]*lpssk;
      ggv[0][1]=ggv[1][0];
      ggv[0][2]=ggv[2][0];
      ggv[1][2]=ggv[2][1];
    }
  }

  if (order2<2) return;
  HierarchicalSidePolynomial hsp;
  NumPtr<REAL> hspv0(order2-1);
  NumPtr<REAL> hspv1(order2-1); 
  NumPtr<REAL> hspv2(order2-1);
  hsp.values(bc[0],hspv0);
  hsp.values(bc[1],hspv1);
  hsp.values(bc[2],hspv2);
  NumPtr<REAL> hsps0(order2-1);
  NumPtr<REAL> hsps1(order2-1); 
  NumPtr<REAL> hsps2(order2-1);
  hsp.slopes(bc[0],hsps0);
  hsp.slopes(bc[1],hsps1);
  hsp.slopes(bc[2],hsps2);
  NumPtr<REAL> hspss0(order2-1);
  NumPtr<REAL> hspss1(order2-1); 
  NumPtr<REAL> hspss2(order2-1);
  hsp.slope2s(bc[0],hspss0);
  hsp.slope2s(bc[1],hspss1);
  hsp.slope2s(bc[2],hspss2);
  Tensor<1,2> gab;
  Tensor<2,2> ggab;
  REAL lpvk=lpv[0];
  REAL lpsk=lps[0];
  REAL lpssk=lpss[0];
  {
//  line 0
    REAL ab=bc[2]*bc[0];
    gab[0]=bc[2]-bc[0];
    gab[1]=-bc[0];
    ggab[0][0]=-2.;
    ggab[1][0]=-1.;
    ggab[0][1]=-1.;
    ggab[1][1]=0.;
    if (ev[0]<ev[1]) {
      for (int i=1;i<order2;i++) {
        Tensor<2,3> &ggv=gg[iv++];
        ggv[0][0]=(ggab[0][0]*hspv0[i-1]+2.*gab[0]*hsps0[i-1]
                  +ab*hspss0[i-1])*lpvk;
        ggv[1][0]=(ggab[1][0]*hspv0[i-1]+gab[1]*hsps0[i-1])*lpvk;
        ggv[2][0]=(gab[0]*hspv0[i-1]+ab*hsps0[i-1])*lpsk;
        ggv[1][1]=ggab[1][1]*hspv0[i-1]*lpvk;
        ggv[2][1]=gab[1]*hspv0[i-1]*lpsk;
        ggv[2][2]=ab*hspv0[i-1]*lpssk;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    } else {
      OBSOLETE("not programmed");
    }
  }
  {
//  line 1
    REAL ab=bc[0]*bc[1];
    gab[0]=bc[1];
    gab[1]=bc[0];
    ggab[0][0]=0.;
    ggab[1][0]=1.;
    ggab[0][1]=1.;
    ggab[1][1]=0.;
    if (ev[1]<ev[2]) {
      for (int j=1;j<order2;j++) {
        Tensor<2,3> &ggv=gg[iv++];
        ggv[0][0]=ggab[0][0]*hspv1[j-1]*lpvk;
        ggv[1][0]=(ggab[1][0]*hspv1[j-1]+gab[0]*hsps1[j-1])*lpvk;
        ggv[2][0]=gab[0]*hspv1[j-1]*lpsk;
        ggv[1][1]=(ggab[1][1]*hspv1[j-1]+2.*gab[1]*hsps1[j-1]
                  +ab*hspss1[j-1])*lpvk;
        ggv[2][1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*lpsk;
        ggv[2][2]=ab*hspv1[j-1]*lpssk;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    } else {
      for (int j=1;j<order2;j++) {
        Tensor<2,3> &ggv=gg[iv++];
        ggv[0][0]=(ggab[0][0]*hspv0[j-1]+2.*gab[0]*hsps0[j-1]
                  +ab*hspss0[j-1])*lpvk;
        ggv[1][0]=(ggab[1][0]*hspv0[j-1]+gab[1]*hsps0[j-1])*lpvk;
        ggv[2][0]=(gab[0]*hspv0[j-1]+ab*hsps0[j-1])*lpsk;
        ggv[1][1]=ggab[1][1]*hspv0[j-1]*lpvk;
        ggv[2][1]=gab[1]*hspv0[j-1]*lpsk;
        ggv[2][2]=ab*hspv0[j-1]*lpssk;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    }
  }
  {
//  line 2
    REAL ab=bc[2]*bc[1];
    gab[0]=-bc[1];
    gab[1]=bc[2]-bc[1];
    ggab[0][0]=0.;
    ggab[1][0]=-1.;
    ggab[0][1]=-1.;
    ggab[1][1]=-2.;
    if (ev[0]<ev[2]) {
      for (int j=1;j<order2;j++) {
        Tensor<2,3> &ggv=gg[iv++];
        ggv[0][0]=ggab[0][0]*hspv1[j-1]*lpvk;
        ggv[1][0]=(ggab[1][0]*hspv1[j-1]+gab[0]*hsps1[j-1])*lpvk;
        ggv[2][0]=gab[0]*hspv1[j-1]*lpsk;
        ggv[1][1]=(ggab[1][1]*hspv1[j-1]+2.*gab[1]*hsps1[j-1]
                  +ab*hspss1[j-1])*lpvk;
        ggv[2][1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*lpsk;
        ggv[2][2]=ab*hspv1[j-1]*lpssk;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    } else { // not used
      OBSOLETE("not programmed");
    }
  }
  {
//  line 3
    lpvk=lpv[1];
    lpsk=lps[1];
    lpssk=lpss[1];
    REAL ab=bc[2]*bc[0];
    gab[0]=bc[2]-bc[0];
    gab[1]=-bc[0];
    ggab[0][0]=-2.;
    ggab[1][0]=-1.;
    ggab[0][1]=-1.;
    ggab[1][1]=0.;
    if (ev[3]<ev[4]) {
      for (int i=1;i<order2;i++) {
        Tensor<2,3> &ggv=gg[iv++];
        ggv[0][0]=(ggab[0][0]*hspv0[i-1]+2.*gab[0]*hsps0[i-1]
                  +ab*hspss0[i-1])*lpvk;
        ggv[1][0]=(ggab[1][0]*hspv0[i-1]+gab[1]*hsps0[i-1])*lpvk;
        ggv[2][0]=(gab[0]*hspv0[i-1]+ab*hsps0[i-1])*lpsk;
        ggv[1][1]=ggab[1][1]*hspv0[i-1]*lpvk;
        ggv[2][1]=gab[1]*hspv0[i-1]*lpsk;
        ggv[2][2]=ab*hspv0[i-1]*lpssk;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    } else {
      for (int i=1;i<order2;i++) {
        Tensor<2,3> &ggv=gg[iv++];
        ggv[0][0]=(ggab[0][0]*hspv2[i-1]-2.*gab[0]*hsps2[i-1]
                  +ab*hspss2[i-1])*lpvk;
        ggv[1][0]=(ggab[1][0]*hspv2[i-1]-gab[1]*hsps2[i-1]
                  -gab[0]*hsps2[i-1]+ab*hspss2[i-1])*lpvk;
        ggv[2][0]=(gab[0]*hspv2[i-1]-ab*hsps2[i-1])*lpsk;
        ggv[1][1]=(ggab[1][1]*hspv2[i-1]-2.*gab[1]*hsps2[i-1]
                  +ab*hspss2[i-1])*lpvk;
        ggv[2][1]=(gab[1]*hspv2[i-1]-ab*hsps2[i-1])*lpsk;
        ggv[2][2]=ab*hspv2[i-1]*lpssk;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    }
  }
  {
//  line 4
    REAL ab=bc[0]*bc[1];
    gab[0]=bc[1];
    gab[1]=bc[0];
    ggab[0][0]=0.;
    ggab[1][0]=1.;
    ggab[0][1]=1.;
    ggab[1][1]=0.;
    if (ev[4]<ev[5]) {
      for (int j=1;j<order2;j++) {
        Tensor<2,3> &ggv=gg[iv++];
        ggv[0][0]=ggab[0][0]*hspv1[j-1]*lpvk;
        ggv[1][0]=(ggab[1][0]*hspv1[j-1]+gab[0]*hsps1[j-1])*lpvk;
        ggv[2][0]=gab[0]*hspv1[j-1]*lpsk;
        ggv[1][1]=(ggab[1][1]*hspv1[j-1]+2.*gab[1]*hsps1[j-1]
                  +ab*hspss1[j-1])*lpvk;
        ggv[2][1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*lpsk;
        ggv[2][2]=ab*hspv1[j-1]*lpssk;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    } else {
      for (int j=1;j<order2;j++) {
        Tensor<2,3> &ggv=gg[iv++];
        ggv[0][0]=(ggab[0][0]*hspv0[j-1]+2.*gab[0]*hsps0[j-1]
                  +ab*hspss0[j-1])*lpvk;
        ggv[1][0]=(ggab[1][0]*hspv0[j-1]+gab[1]*hsps0[j-1])*lpvk;
        ggv[2][0]=(gab[0]*hspv0[j-1]+ab*hsps0[j-1])*lpsk;
        ggv[1][1]=ggab[1][1]*hspv0[j-1]*lpvk;
        ggv[2][1]=gab[1]*hspv0[j-1]*lpsk;
        ggv[2][2]=ab*hspv0[j-1]*lpssk;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    }
  }
  {
//  line 5
    REAL ab=bc[2]*bc[1];
    gab[0]=-bc[1];
    gab[1]=bc[2]-bc[1];
    ggab[0][0]=0.;
    ggab[1][0]=-1.;
    ggab[0][1]=-1.;
    ggab[1][1]=-2.;
    if (ev[3]<ev[5]) {
      for (int j=1;j<order2;j++) {
        Tensor<2,3> &ggv=gg[iv++];
        ggv[0][0]=ggab[0][0]*hspv1[j-1]*lpvk;
        ggv[1][0]=(ggab[1][0]*hspv1[j-1]+gab[0]*hsps1[j-1])*lpvk;
        ggv[2][0]=gab[0]*hspv1[j-1]*lpsk;
        ggv[1][1]=(ggab[1][1]*hspv1[j-1]+2.*gab[1]*hsps1[j-1]
                  +ab*hspss1[j-1])*lpvk;
        ggv[2][1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*lpsk;
        ggv[2][2]=ab*hspv1[j-1]*lpssk;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    } else {
      for (int j=1;j<order2;j++) {
        Tensor<2,3> &ggv=gg[iv++];
        ggv[0][0]=(ggab[0][0]*hspv2[j-1]-2.*gab[0]*hsps2[j-1]
                  +ab*hspss2[j-1])*lpvk;
        ggv[1][0]=(ggab[1][0]*hspv2[j-1]-gab[1]*hsps2[j-1]
                  -gab[0]*hsps2[j-1]+ab*hspss2[j-1])*lpvk;
        ggv[2][0]=(gab[0]*hspv2[j-1]-ab*hsps2[j-1])*lpsk;
        ggv[1][1]=(ggab[1][1]*hspv2[j-1]-2.*gab[1]*hsps2[j-1]
                  +ab*hspss2[j-1])*lpvk;
        ggv[2][1]=(gab[1]*hspv2[j-1]-ab*hsps2[j-1])*lpsk;
        ggv[2][2]=ab*hspv2[j-1]*lpssk;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    }
  }
  {
//  line 6
    if (ev[0]<ev[3]) {
      for (int k=2;k<=order2;k++) {
        Tensor<2,3> &ggv=gg[iv++];
        ggv[0][0]=0.;
        ggv[1][0]=0.;
        ggv[2][0]=-lps[k];
        ggv[1][1]=0.;
        ggv[2][1]=-lps[k];
        ggv[2][2]=bc[2]*lpss[k];
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    } else { // not used
      OBSOLETE("not programmed");
    }
  }
  {
//  line 7
    if (ev[1]<ev[4]) {
      for (int k=2;k<=order2;k++) {
        Tensor<2,3> &ggv=gg[iv++];
        ggv[0][0]=0.;
        ggv[1][0]=0.;
        ggv[2][0]=lps[k];
        ggv[1][1]=0.;
        ggv[2][1]=0.;
        ggv[2][2]=bc[0]*lpss[k];
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    } else {
      REAL s=1.;
      for (int k=2;k<=order2;k++) {
        Tensor<2,3> &ggv=gg[iv++];
        ggv[0][0]=0.;
        ggv[1][0]=0.;
        ggv[2][0]=lps[k]*s;
        ggv[1][1]=0.;
        ggv[2][1]=0.;
        ggv[2][2]=bc[0]*lpss[k]*s;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
        s=-s;
      }
    }
  }
  {
//  line 8
    if (ev[2]<ev[5]) {
      for (int k=2;k<=order2;k++) {
        Tensor<2,3> &ggv=gg[iv++];
        ggv[0][0]=0.;
        ggv[1][0]=0.;
        ggv[2][0]=0.;
        ggv[1][1]=0.;
        ggv[2][1]=lps[k];
        ggv[2][2]=bc[1]*lpss[k];
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    } else {
      REAL s=1.;
      for (int k=2;k<=order2;k++) {
        Tensor<2,3> &ggv=gg[iv++];
        ggv[0][0]=0.;
        ggv[1][0]=0.;
        ggv[2][0]=0.;
        ggv[1][1]=0.;
        ggv[2][1]=lps[k]*s;
        ggv[2][2]=bc[1]*lpss[k]*s;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
        s=-s;
      }
    }
  }

  if (order2>=3) {
    LegendrePolynomial lp;
    NumPtr<REAL> lpv0(order2-2);
    NumPtr<REAL> lpv1(order2-2);
    NumPtr<REAL> lpv2(order2-2);
    lp.values(bc[0],lpv0);
    lp.values(bc[1],lpv1);
    lp.values(bc[2],lpv2);
    NumPtr<REAL> lps0(order2-2);
    NumPtr<REAL> lps1(order2-2);
    NumPtr<REAL> lps2(order2-2);
    lp.slopes(bc[0],lps0);
    lp.slopes(bc[1],lps1);
    lp.slopes(bc[2],lps2);
    NumPtr<REAL> lpss0(order2-2);
    NumPtr<REAL> lpss1(order2-2);
    NumPtr<REAL> lpss2(order2-2);
    lp.slope2s(bc[0],lpss0);
    lp.slope2s(bc[1],lpss1);
    lp.slope2s(bc[2],lpss2);
    REAL abc=bc[2]*bc[0]*bc[1];
    Tensor<1,2> gabc;
    gabc[0]=(bc[2]-bc[0])*bc[1];
    gabc[1]=(bc[2]-bc[1])*bc[0];
    Tensor<2,2> ggabc;
    ggabc[0][0]=-2.*bc[1];
    ggabc[1][0]=bc[2]-bc[0]-bc[1];
    ggabc[0][1]=ggabc[1][0];
    ggabc[1][1]=-2.*bc[0];
    {
//    face 0
      lpvk=lpv[0];
      lpsk=lps[0];
      lpssk=lpss[0];
      if (ev[0]<ev[1]) {
        if (ev[1]<ev[2]) {
          for (int m=0;m+3<=order2;m++) {
            for (int i=0;i<=m;i++) {
              Tensor<2,3> &ggv=gg[iv++];
              ggv[0][0]=(ggabc[0][0]*lpv0[m-i]+2.*gabc[0]*lps0[m-i]
                        +abc*lpss0[m-i])*lpv1[i]*lpvk;
              ggv[1][0]=((ggabc[1][0]*lpv1[i]+gabc[0]*lps1[i])*lpv0[m-i]
                        +(gabc[1]*lpv1[i]+abc*lps1[i])*lps0[m-i])*lpvk;
              ggv[2][0]=(gabc[0]*lpv0[m-i]+abc*lps0[m-i])*lpv1[i]*lpsk;
              ggv[1][1]=(ggabc[1][1]*lpv1[i]+2.*gabc[1]*lps1[i]
                        +abc*lpss1[i])*lpv0[m-i]*lpvk;
              ggv[2][1]=(gabc[1]*lpv1[i]+abc*lps1[i])*lpv0[m-i]*lpsk;
              ggv[2][2]=abc*lpv0[m-i]*lpv1[i]*lpssk;
              ggv[0][1]=ggv[1][0];
              ggv[0][2]=ggv[2][0];
              ggv[1][2]=ggv[2][1];
            }
          }
        } else if (ev[0]<ev[2]) {
          for (int m=0;m+3<=order2;m++) {
            for (int i=0;i<=m;i++) {
              Tensor<2,3> &ggv=gg[iv++];
              ggv[0][0]=(ggabc[0][0]*lpv0[i]+2.*gabc[0]*lps0[i]
                        +abc*lpss0[i])*lpv1[m-i]*lpvk;
              ggv[1][0]=((ggabc[1][0]*lpv1[m-i]+gabc[0]*lps1[m-i])*lpv0[i]
                        +(gabc[1]*lpv1[m-i]+abc*lps1[m-i])*lps0[i])*lpvk;
              ggv[2][0]=(gabc[0]*lpv0[i]+abc*lps0[i])*lpv1[m-i]*lpsk;
              ggv[1][1]=(ggabc[1][1]*lpv1[m-i]+2.*gabc[1]*lps1[m-i]
                        +abc*lpss1[m-i])*lpv0[i]*lpvk;
              ggv[2][1]=(gabc[1]*lpv1[m-i]+abc*lps1[m-i])*lpv0[i]*lpsk;
              ggv[2][2]=abc*lpv1[m-i]*lpv0[i]*lpssk;
              ggv[0][1]=ggv[1][0];
              ggv[0][2]=ggv[2][0];
              ggv[1][2]=ggv[2][1];
            }
          }
        } else { // not used
          OBSOLETE("not programmed");
        }
      } else if (ev[0]<ev[2]) { // not used
        OBSOLETE("not programmed");
      } else if (ev[1]<ev[2]) { // not used
        OBSOLETE("not programmed");
      } else { // not used
        OBSOLETE("not programmed");
      }
    }
    {
//    face 1
      lpvk=lpv[1];
      lpsk=lps[1];
      lpssk=lpss[1];
      if (ev[3]<ev[4]) {
        if (ev[4]<ev[5]) {
          for (int m=0;m+3<=order2;m++) {
            for (int i=0;i<=m;i++) {
              Tensor<2,3> &ggv=gg[iv++];
              ggv[0][0]=(ggabc[0][0]*lpv0[m-i]+2.*gabc[0]*lps0[m-i]
                        +abc*lpss0[m-i])*lpv1[i]*lpvk;
              ggv[1][0]=((ggabc[1][0]*lpv1[i]+gabc[0]*lps1[i])*lpv0[m-i]
                        +(gabc[1]*lpv1[i]+abc*lps1[i])*lps0[m-i])*lpvk;
              ggv[2][0]=(gabc[0]*lpv0[m-i]+abc*lps0[m-i])*lpv1[i]*lpsk;
              ggv[1][1]=(ggabc[1][1]*lpv1[i]+2.*gabc[1]*lps1[i]
                        +abc*lpss1[i])*lpv0[m-i]*lpvk;
              ggv[2][1]=(gabc[1]*lpv1[i]+abc*lps1[i])*lpv0[m-i]*lpsk;
              ggv[2][2]=abc*lpv0[m-i]*lpv1[i]*lpssk;
              ggv[0][1]=ggv[1][0];
              ggv[0][2]=ggv[2][0];
              ggv[1][2]=ggv[2][1];
            }
          }
        } else if (ev[3]<ev[5]) {
          for (int m=0;m+3<=order2;m++) {
            for (int i=0;i<=m;i++) {
              Tensor<2,3> &ggv=gg[iv++];
              ggv[0][0]=(ggabc[0][0]*lpv0[i]+2.*gabc[0]*lps0[i]
                        +abc*lpss0[i])*lpv1[m-i]*lpvk;
              ggv[1][0]=((ggabc[1][0]*lpv1[m-i]+gabc[0]*lps1[m-i])*lpv0[i]
                        +(gabc[1]*lpv1[m-i]+abc*lps1[m-i])*lps0[i])*lpvk;
              ggv[2][0]=(gabc[0]*lpv0[i]+abc*lps0[i])*lpv1[m-i]*lpsk;
              ggv[1][1]=(ggabc[1][1]*lpv1[m-i]+2.*gabc[1]*lps1[m-i]
                        +abc*lpss1[m-i])*lpv0[i]*lpvk;
              ggv[2][1]=(gabc[1]*lpv1[m-i]+abc*lps1[m-i])*lpv0[i]*lpsk;
              ggv[2][2]=abc*lpv1[m-i]*lpv0[i]*lpssk;
              ggv[0][1]=ggv[1][0];
              ggv[0][2]=ggv[2][0];
              ggv[1][2]=ggv[2][1];
            }
          }
        } else {
          for (int m=0;m+3<=order2;m++) {
            for (int i=0;i<=m;i++) {
              Tensor<2,3> &ggv=gg[iv++];
              ggv[0][0]=(ggabc[0][0]*lpv2[m-i]*lpv0[i]
                        -2.*gabc[0]*lps2[m-i]*lpv0[i]
                        +2.*gabc[0]*lpv2[m-i]*lps0[i]
                        +abc*(lpss2[m-i]*lpv0[i]-2.*lps2[m-i]*lps0[i]
                             +lpv2[m-i]*lpss0[i]))*lpvk;
              ggv[1][0]=((ggabc[1][0]*lpv2[m-i]-gabc[1]*lps2[m-i]
                         -gabc[0]*lps2[m-i]+abc*lpss2[m-i])*lpv0[i]
                        +(gabc[1]*lpv2[m-i]-abc*lps2[m-i])*lps0[i])*lpvk;
              ggv[2][0]=(gabc[0]*lpv2[m-i]*lpv0[i]
                        -abc*lps2[m-i]*lpv0[i]
                        +abc*lpv2[m-i]*lps0[i])*lpsk;
              ggv[1][1]=(ggabc[1][1]*lpv2[m-i]-2.*gabc[1]*lps2[m-i]
                        +abc*lpss2[m-i])*lpv0[i]*lpvk;
              ggv[2][1]=(gabc[1]*lpv2[m-i]-abc*lps2[m-i])*lpv0[i]*lpsk;
              ggv[2][2]=abc*lpv2[m-i]*lpv0[i]*lpssk;
              ggv[0][1]=ggv[1][0];
              ggv[0][2]=ggv[2][0];
              ggv[1][2]=ggv[2][1];
            }
          }
        }
      } else if (ev[3]<ev[5]) {
        for (int m=0;m+3<=order2;m++) {
          for (int i=0;i<=m;i++) {
            Tensor<2,3> &ggv=gg[iv++];
            ggv[0][0]=(ggabc[0][0]*lpv2[m-i]-2.*gabc[0]*lps2[m-i]
                      +abc*lpss2[m-i])*lpv1[i]*lpvk;
            ggv[1][0]=((ggabc[1][0]*lpv2[m-i]-gabc[1]*lps2[m-i])*lpv1[i]
                      +gabc[0]*(-lps2[m-i]*lpv1[i]+lpv2[m-i]*lps1[i])
                      +abc*(lpss2[m-i]*lpv1[i]-lps2[m-i]*lps1[i]))*lpvk;
            ggv[2][0]=(gabc[0]*lpv2[m-i]-abc*lps2[m-i])*lpv1[i]*lpsk;
            ggv[1][1]=(ggabc[1][1]*lpv2[m-i]*lpv1[i]
                      -2.*gabc[1]*lps2[m-i]*lpv1[i]
                      +2.*gabc[1]*lpv2[m-i]*lps1[i]
                      +abc*(lpss2[m-i]*lpv1[i]-2.*lps2[m-i]*lps1[i]
                           +lpv2[m-i]*lpss1[i]))*lpvk;
            ggv[2][1]=(gabc[1]*lpv2[m-i]*lpv1[i]
                      +abc*(-lps2[m-i]*lpv1[i]+lpv2[m-i]*lps1[i]))*lpsk;
            ggv[2][2]=abc*lpv2[m-i]*lpv1[i]*lpssk;
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
          }
        }
      } else if (ev[4]<ev[5]) {
        for (int m=0;m+3<=order2;m++) {
          for (int i=0;i<=m;i++) {
            Tensor<2,3> &ggv=gg[iv++];
            ggv[0][0]=(ggabc[0][0]*lpv2[i]-2.*gabc[0]*lps2[i]
                      +abc*lpss2[i])*lpv1[m-i]*lpvk;
            ggv[1][0]=((ggabc[1][0]*lpv2[i]-gabc[1]*lps2[i])*lpv1[m-i]
                      +gabc[0]*(lps1[m-i]*lpv2[i]-lpv1[m-i]*lps2[i])
                      +abc*(-lps1[m-i]*lps2[i]+lpv1[m-i]*lpss2[i]))*lpvk;
            ggv[2][0]=(gabc[0]*lpv2[i]-abc*lps2[i])*lpv1[m-i]*lpsk;
            ggv[1][1]=(ggabc[1][1]*lpv2[i]*lpv1[m-i]
                      -gabc[1]*lps2[i]*lpv1[m-i]
                      +gabc[1]*lpv2[i]*lps1[m-i]
                      +gabc[1]*(lps1[m-i]*lpv2[i]-lpv1[m-i]*lps2[i])
                      +abc*(lpss1[m-i]*lpv2[i]-2.*lps1[m-i]*lps2[i]
                           +lpv1[m-i]*lpss2[i]))*lpvk;
            ggv[2][1]=(gabc[1]*lpv1[m-i]*lpv2[i]
                      +abc*lps1[m-i]*lpv2[i]
                      -abc*lpv1[m-i]*lps2[i])*lpsk;
            ggv[2][2]=abc*lpv1[m-i]*lpv2[i]*lpssk;
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
          }
        }
      } else {
        for (int m=0;m+3<=order2;m++) {
          for (int i=0;i<=m;i++) {
            Tensor<2,3> &ggv=gg[iv++];
            ggv[0][0]=(ggabc[0][0]*lpv2[i]*lpv0[m-i]
                      -2.*gabc[0]*lps2[i]*lpv0[m-i]
                      +2.*gabc[0]*lpv2[i]*lps0[m-i]
                      +abc*(lpss0[m-i]*lpv2[i]-2.*lps0[m-i]*lps2[i]
                           +lpv0[m-i]*lpss2[i]))*lpvk;
            ggv[1][0]=((ggabc[1][0]*lpv2[i]-gabc[1]*lps2[i]
                       -gabc[0]*lps2[i]+abc*lpss2[i])*lpv0[m-i]
                      +(gabc[1]*lpv2[i]-abc*lps2[i])*lps0[m-i])*lpvk;
            ggv[2][0]=(gabc[0]*lpv0[m-i]*lpv2[i]
                      +abc*lps0[m-i]*lpv2[i]
                      -abc*lpv0[m-i]*lps2[i])*lpsk;
            ggv[1][1]=(ggabc[1][1]*lpv2[i]-2.*gabc[1]*lps2[i]
                      +abc*lpss2[i])*lpv0[m-i]*lpvk;
            ggv[2][1]=(gabc[1]*lpv2[i]-abc*lps2[i])*lpv0[m-i]*lpsk;
            ggv[2][2]=abc*lpv0[m-i]*lpv2[i]*lpssk;
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
          }
        }
      }
    }
  }
  {
//  face 2
    int vmin=min(min(ev[0],ev[2]),min(ev[3],ev[5]));
    REAL ab=bc[2]*bc[1];
    gab[0]=-bc[1];
    gab[1]=bc[2]-bc[1];
    ggab[0][0]=0.;
    ggab[1][0]=-1.;
    ggab[0][1]=-1.;
    ggab[1][1]=-2.;
    if (ev[0]==vmin) {
      if (ev[2]<ev[3]) {
        for (int k=2;k<=order2;k++) {
          for (int j=1;j<order2;j++) {
            Tensor<2,3> &ggv=gg[iv++];
            ggv[0][0]=ggab[0][0]*hspv1[j-1]*lpv[k];
            ggv[1][0]=(ggab[1][0]*hspv1[j-1]+gab[0]*hsps1[j-1])*lpv[k];
            ggv[2][0]=gab[0]*hspv1[j-1]*lps[k];
            ggv[1][1]=(ggab[1][1]*hspv1[j-1]+2.*gab[1]*hsps1[j-1]
                      +ab*hspss1[j-1])*lpv[k];
            ggv[2][1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*lps[k];
            ggv[2][2]=ab*hspv1[j-1]*lpss[k];
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
          }
        }
      } else {
        for (int j=1;j<order2;j++) {
          for (int k=2;k<=order2;k++) {
            Tensor<2,3> &ggv=gg[iv++];
            ggv[0][0]=ggab[0][0]*hspv1[j-1]*lpv[k];
            ggv[1][0]=(ggab[1][0]*hspv1[j-1]+gab[0]*hsps1[j-1])*lpv[k];
            ggv[2][0]=gab[0]*hspv1[j-1]*lps[k];
            ggv[1][1]=(ggab[1][1]*hspv1[j-1]+2.*gab[1]*hsps1[j-1]
                      +ab*hspss1[j-1])*lpv[k];
            ggv[2][1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*lps[k];
            ggv[2][2]=ab*hspv1[j-1]*lpss[k];
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
          }
        }
      }
    } else if (ev[2]==vmin) {
      if (ev[0]<ev[5]) { // not used
        OBSOLETE("not programmed");
      } else { // not used
        OBSOLETE("not programmed");
      }
    } else if (ev[3]==vmin) {
      if (ev[0]<ev[5]) { // not used
        OBSOLETE("not programmed");
      } else { // not used
        OBSOLETE("not programmed");
      }
    } else {
      if (ev[2]<ev[3]) { // not used
        OBSOLETE("not programmed");
      } else { // not used
        OBSOLETE("not programmed");
      }
    }
  }
  {
//  face 3
    int vmin=min(min(ev[0],ev[1]),min(ev[3],ev[4]));
    REAL ab=bc[2]*bc[0];
    gab[0]=bc[2]-bc[0];
    gab[1]=-bc[0];
    ggab[0][0]=-2.;
    ggab[1][0]=-1.;
    ggab[0][1]=-1.;
    ggab[1][1]=0.;
    if (ev[0]==vmin) {
      if (ev[1]<ev[3]) {
        for (int k=2;k<=order2;k++) {
          for (int i=1;i<order2;i++) {
            Tensor<2,3> &ggv=gg[iv++];
            ggv[0][0]=(ggab[0][0]*hspv0[i-1]+2.*gab[0]*hsps0[i-1]
                      +ab*hspss0[i-1])*lpv[k];
            ggv[1][0]=(ggab[1][0]*hspv0[i-1]+gab[1]*hsps0[i-1])*lpv[k];
            ggv[2][0]=(gab[0]*hspv0[i-1]+ab*hsps0[i-1])*lps[k];
            ggv[1][1]=ggab[1][1]*hspv0[i-1]*lpv[k];
            ggv[2][1]=gab[1]*hspv0[i-1]*lps[k];
            ggv[2][2]=ab*hspv0[i-1]*lpss[k];
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
          }
        }
      } else {
        for (int i=1;i<order2;i++) {
          for (int k=2;k<=order2;k++) {
            Tensor<2,3> &ggv=gg[iv++];
            ggv[0][0]=(ggab[0][0]*hspv0[i-1]+2.*gab[0]*hsps0[i-1]
                      +ab*hspss0[i-1])*lpv[k];
            ggv[1][0]=(ggab[1][0]*hspv0[i-1]+gab[1]*hsps0[i-1])*lpv[k];
            ggv[2][0]=(gab[0]*hspv0[i-1]+ab*hsps0[i-1])*lps[k];
            ggv[1][1]=ggab[1][1]*hspv0[i-1]*lpv[k];
            ggv[2][1]=gab[1]*hspv0[i-1]*lps[k];
            ggv[2][2]=ab*hspv0[i-1]*lpss[k];
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
          }
        }
      }
    } else if (ev[1]==vmin) {
      if (ev[0]<ev[4]) { // not used
        OBSOLETE("not programmed");
      } else { // not used
        OBSOLETE("not programmed");
      }
    } else if (ev[3]==vmin) {
      if (ev[0]<ev[4]) { // not used
        OBSOLETE("not programmed");
      } else { // not used
        OBSOLETE("not programmed");
      }
    } else {
      if (ev[1]<ev[3]) { // not used
        OBSOLETE("not programmed");
      } else { // not used
        OBSOLETE("not programmed");
      }
    }
  }
  {
//  face 4
    int vmin=min(min(ev[1],ev[2]),min(ev[4],ev[5]));
    REAL ab=bc[0]*bc[1];
    gab[0]=bc[1];
    gab[1]=bc[0];
    ggab[0][0]=0.;
    ggab[1][0]=1.;
    ggab[0][1]=1.;
    ggab[1][1]=0.;
    if (ev[1]==vmin) {
      if (ev[2]<ev[4]) {
        for (int k=2;k<=order2;k++) {
          for (int j=1;j<order2;j++) {
            Tensor<2,3> &ggv=gg[iv++];
            ggv[0][0]=ggab[0][0]*hspv1[j-1]*lpv[k];
            ggv[1][0]=(ggab[1][0]*hspv1[j-1]+gab[0]*hsps1[j-1])*lpv[k];
            ggv[2][0]=gab[0]*hspv1[j-1]*lps[k];
            ggv[1][1]=(ggab[1][1]*hspv1[j-1]+2.*gab[1]*hsps1[j-1]
                      +ab*hspss1[j-1])*lpv[k];
            ggv[2][1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*lps[k];
            ggv[2][2]=ab*hspv1[j-1]*lpss[k];
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
          }
        }
      } else {
        for (int j=1;j<order2;j++) {
          for (int k=2;k<=order2;k++) {
            Tensor<2,3> &ggv=gg[iv++];
            ggv[0][0]=ggab[0][0]*hspv1[j-1]*lpv[k];
            ggv[1][0]=(ggab[1][0]*hspv1[j-1]+gab[0]*hsps1[j-1])*lpv[k];
            ggv[2][0]=gab[0]*hspv1[j-1]*lps[k];
            ggv[1][1]=(ggab[1][1]*hspv1[j-1]+2.*gab[1]*hsps1[j-1]
                      +ab*hspss1[j-1])*lpv[k];
            ggv[2][1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*lps[k];
            ggv[2][2]=ab*hspv1[j-1]*lpss[k];
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
          }
        }
      }
    } else if (ev[2]==vmin) {
      if (ev[1]<ev[5]) {
        for (int k=2;k<=order2;k++) {
          REAL s=1.;
          for (int j=1;j<order2;j++) {
            Tensor<2,3> &ggv=gg[iv++];
            ggv[0][0]=ggab[0][0]*hspv1[j-1]*s*lpv[k];
            ggv[1][0]=(ggab[1][0]*hspv1[j-1]+gab[0]*hsps1[j-1])*s*lpv[k];
            ggv[2][0]=gab[0]*hspv1[j-1]*s*lps[k];
            ggv[1][1]=(ggab[1][1]*hspv1[j-1]+2.*gab[1]*hsps1[j-1]
                      +ab*hspss1[j-1])*s*lpv[k];
            ggv[2][1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*s*lps[k];
            ggv[2][2]=ab*hspv1[j-1]*s*lpss[k];
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
            s=-s;
          }
        }
      } else {
        REAL s=1.;
        for (int j=1;j<order2;j++) {
          for (int k=2;k<=order2;k++) {
            Tensor<2,3> &ggv=gg[iv++];
            ggv[0][0]=ggab[0][0]*hspv1[j-1]*s*lpv[k];
            ggv[1][0]=(ggab[1][0]*hspv1[j-1]+gab[0]*hsps1[j-1])*s*lpv[k];
            ggv[2][0]=gab[0]*hspv1[j-1]*s*lps[k];
            ggv[1][1]=(ggab[1][1]*hspv1[j-1]+2.*gab[1]*hsps1[j-1]
                      +ab*hspss1[j-1])*s*lpv[k];
            ggv[2][1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*s*lps[k];
            ggv[2][2]=ab*hspv1[j-1]*s*lpss[k];
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
          }
          s=-s;
        }
      }
    } else if (ev[4]==vmin) {
      if (ev[1]<ev[5]) {
        for (int j=1;j<order2;j++) {
          REAL s=1.;
          for (int k=2;k<=order2;k++) {
            Tensor<2,3> &ggv=gg[iv++];
            ggv[0][0]=ggab[0][0]*hspv1[j-1]*lpv[k]*s;
            ggv[1][0]=(ggab[1][0]*hspv1[j-1]+gab[0]*hsps1[j-1])*lpv[k]*s;
            ggv[2][0]=gab[0]*hspv1[j-1]*lps[k]*s;
            ggv[1][1]=(ggab[1][1]*hspv1[j-1]+2.*gab[1]*hsps1[j-1]
                      +ab*hspss1[j-1])*lpv[k]*s;
            ggv[2][1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*lps[k]*s;
            ggv[2][2]=ab*hspv1[j-1]*lpss[k]*s;
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
            s=-s;
          }
        }
      } else {
        REAL s=1.;
        for (int k=2;k<=order2;k++) {
          for (int j=1;j<order2;j++) {
            Tensor<2,3> &ggv=gg[iv++];
            ggv[0][0]=ggab[0][0]*hspv1[j-1]*lpv[k]*s;
            ggv[1][0]=(ggab[1][0]*hspv1[j-1]+gab[0]*hsps1[j-1])*lpv[k]*s;
            ggv[2][0]=gab[0]*hspv1[j-1]*lps[k]*s;
            ggv[1][1]=(ggab[1][1]*hspv1[j-1]+2.*gab[1]*hsps1[j-1]
                      +ab*hspss1[j-1])*lpv[k]*s;
            ggv[2][1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*lps[k]*s;
            ggv[2][2]=ab*hspv1[j-1]*lpss[k]*s;
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
          }
          s=-s;
        }
      }
    } else {
      if (ev[2]<ev[4]) {
        REAL sj=1.;
        for (int j=1;j<order2;j++) {
          REAL sk=1.;
          for (int k=2;k<=order2;k++) {
            Tensor<2,3> &ggv=gg[iv++];
            ggv[0][0]=ggab[0][0]*hspv1[j-1]*sj*lpv[k]*sk;
            ggv[1][0]=(ggab[1][0]*hspv1[j-1]
                      +gab[0]*hsps1[j-1])*sj*lpv[k]*sk;
            ggv[2][0]=gab[0]*hspv1[j-1]*sj*lps[k]*sk;
            ggv[1][1]=(ggab[1][1]*hspv1[j-1]+2.*gab[1]*hsps1[j-1]
                      +ab*hspss1[j-1])*sj*lpv[k]*sk;
            ggv[2][1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*sj*lps[k]*sk;
            ggv[2][2]=ab*hspv1[j-1]*sj*lpss[k]*sk;
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
            sk=-sk;
          }
          sj=-sj;
        }
      } else {
        REAL sk=1.;
        for (int k=2;k<=order2;k++) {
          REAL sj=1.;
          for (int j=1;j<order2;j++) {
            Tensor<2,3> &ggv=gg[iv++];
            ggv[0][0]=ggab[0][0]*hspv1[j-1]*sj*lpv[k]*sk;
            ggv[1][0]=(ggab[1][0]*hspv1[j-1]
                      +gab[0]*hsps1[j-1])*sj*lpv[k]*sk;
            ggv[2][0]=gab[0]*hspv1[j-1]*sj*lps[k]*sk;
            ggv[1][1]=(ggab[1][1]*hspv1[j-1]+2.*gab[1]*hsps1[j-1]
                      +ab*hspss1[j-1])*sj*lpv[k]*sk;
            ggv[2][1]=(gab[1]*hspv1[j-1]+ab*hsps1[j-1])*sj*lps[k]*sk;
            ggv[2][2]=ab*hspv1[j-1]*sj*lpss[k]*sk;
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
            sj=-sj;
          }
          sk=-sk;
        }
      }
    }
  }
  if (order2<3) return;
  {
//  interior
    REAL abc=bc[0]*bc[1]*bc[2];
    Tensor<1,2> gabc;
    gabc[0]=(bc[2]-bc[0])*bc[1];
    gabc[1]=(bc[2]-bc[1])*bc[0];
    Tensor<2,2> ggabc;
    ggabc[0][0]=-2.*bc[1];
    ggabc[1][0]=bc[2]-bc[1]-bc[0];
    ggabc[0][1]=ggabc[1][0];
    ggabc[1][1]=-2.*bc[0];
    for (int k=2;k<=order2;k++) {
      for (int m=0;m+3<=order2;m++) {
        for (int j=0;j<=m;j++) {
          Tensor<2,3> &ggv=gg[iv++];
          ggv[0][0]=(ggabc[0][0]*hspv0[m-j]+2.*gabc[0]*hsps0[m-j]
                    +abc*hspss0[m-j])*hspv1[j]*lpv[k];
          ggv[1][0]=((ggabc[1][0]*hspv1[j]+gabc[0]*hsps1[j])*hspv0[m-j]
                    +(gabc[1]*hspv1[j]+abc*hsps1[j])*hsps0[m-j])*lpv[k];
          ggv[2][0]=(gabc[0]*hspv0[m-j]+abc*hsps0[m-j])*hspv1[j]*lps[k];
          ggv[1][1]=(ggabc[1][1]*hspv1[j]+2.*gabc[1]*hsps1[j]
                    +abc*hspss1[j])*hspv0[m-j]*lpv[k];
          ggv[2][1]=(gabc[1]*hspv1[j]+abc*hsps1[j])*hspv0[m-j]*lps[k];
          ggv[2][2]=abc*hspv0[m-j]*hspv1[j]*lpss[k];
          ggv[0][1]=ggv[1][0];
          ggv[0][2]=ggv[2][0];
          ggv[1][2]=ggv[2][1];
        }
      }
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalPrismPolynomials::shapeFunctionsOrder(
const Shape *element,NumPtr<int> &renumber) const {
//triangular faces have shape functions that are not symmetric in
//  barycentric coordinates, so cannot be reordered
//since requiresShapeToComputeValues=true, no reordering is needed
  int n=getNumber();
  renumber.allocate(n);
  for (int i=0;i<n;i++) renumber[i]=i;
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalPrismPolynomials::nonzeroShapesOnFace(
const Shape1 *e,int f,NumPtr<int> &indices) const {
  CHECK_SAME(order01,order2);
  const Shape3 *e3=dynamic_cast<const Shape3*>(e);
  int nv=e->numberVertices();
  int nvof=e3->numberVerticesOnFace(f);
  int lpf=linesPerFace(nv,f);
  int nieslpof=numberInteriorEquallySpacedLatticePointsOnFace(nv,f,order2);
  indices.allocate(nvof+lpf*(order2-1)+nieslpof);
  int iv=0;
  for (int v=0;v<nvof;v++) indices[iv++]=e->elementVertexFromFace(f,v);
  for (int l=0;l<lpf;l++) {
    int offset=nv+elementLineFromFace(nv,f,l)*(order2-1)-1;
    for (int i=1;i<order2;i++) indices[iv++]=offset+i;
  }
  int offset=nv+e3->numberCodimension2Shapes()*(order2-1);
  for (int i=0;i<f;i++) {
    offset+=numberInteriorEquallySpacedLatticePointsOnFace(nv,i,order2);
  }
  for (int i=0;i<nieslpof;i++) indices[iv++]=offset+i;
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalPrismPolynomials::printOn(ostream &os) const {
  os << "HierarchicalPrismPolynomials:" << endl;
  PrismPolynomials::printOn(os);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalHexahedronPolynomials::computeValues(
const Point<SPACEDIM> &p,NumPtr<REAL> &v) const {
  CHECK_TEST(n_pols_1d<=1);
  v.allocate(n_pols);
  NumPtr<NumPtr<REAL> > hpv(3);
  HierarchicalPolynomial hp;
  for (int d=0;d<3;d++) {
    hpv[d].allocate(n_pols_1d);
    hp.values(p[d],hpv[d]);
  }
  NumPtr<REAL> &hpv0=hpv[0];
  NumPtr<REAL> &hpv1=hpv[1];
  NumPtr<REAL> &hpv2=hpv[2];
  int iv=0;
//vertices
  for (int k=0;k<2;k++) {
    REAL hpvk=hpv2[k];
    for (int j=0;j<2;j++) {
      REAL hpvjk=hpv1[j]*hpvk;
      for (int i=0;i<2;i++) {
        v[iv++]=hpv0[i]*hpvjk;
      }
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalHexahedronPolynomials::computeGrads(
const Point<SPACEDIM> &p,NumPtr<Tensor<1,SPACEDIM> > &g) const {
  CHECK_TEST(n_pols_1d<=1);
  g.allocate(n_pols);
#ifdef INDEF
  Tensor<1,3> bogus_tensor;
  bogus_tensor[0]=HUGE_VAL;
  bogus_tensor[1]=HUGE_VAL;
  bogus_tensor[2]=HUGE_VAL;
  g.initialize(bogus_tensor);
#endif
  NumPtr<NumPtr<REAL> > hpv(3);
  NumPtr<NumPtr<REAL> > hps(3);
  HierarchicalPolynomial hp;
  for (int d=0;d<3;d++) {
    hpv[d].allocate(n_pols_1d);
    hps[d].allocate(n_pols_1d);
    hp.values(p[d],hpv[d]);
    hp.slopes(p[d],hps[d]);
  }
  NumPtr<REAL> &hpv0=hpv[0];
  NumPtr<REAL> &hpv1=hpv[1];
  NumPtr<REAL> &hpv2=hpv[2];
  NumPtr<REAL> &hps0=hps[0];
  NumPtr<REAL> &hps1=hps[1];
  NumPtr<REAL> &hps2=hps[2];
  int iv=0;
//vertices
  for (int k=0;k<2;k++) {
    const REAL &hpvk=hpv2[k];
    const REAL &hpsk=hps2[k];
    for (int j=0;j<2;j++) {
      const REAL &hpvj=hpv1[j];
      REAL hpvv=hpvj*hpvk;
      REAL hpsv=hps1[j]*hpvk;
      REAL hpvs=hpvj*hpsk;
      for (int i=0;i<2;i++) {
        Tensor<1,SPACEDIM> &gv=g[iv++];
        const REAL &hpvi=hpv0[i];
        gv[0]=hps0[i]*hpvv;
        gv[1]=hpvi*hpsv;
        gv[2]=hpvi*hpvs;
      }
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalHexahedronPolynomials::computeGradGrads(
const Point<SPACEDIM> &p,NumPtr<Tensor<2,SPACEDIM> > &gg) const {
  CHECK_TEST(n_pols_1d<=1);
  gg.allocate(n_pols);
  NumPtr<NumPtr<REAL> > hpv(3);
  NumPtr<NumPtr<REAL> > hps(3);
  NumPtr<NumPtr<REAL> > hpss(3);
  HierarchicalPolynomial hp;
  for (int d=0;d<3;d++) {
    hpv[d].allocate(n_pols_1d);
    hps[d].allocate(n_pols_1d);
    hpss[d].allocate(n_pols_1d);
    hp.values(p[d],hpv[d]);
    hp.slopes(p[d],hps[d]);
    hp.slope2s(p[d],hpss[d]);
  }
  NumPtr<REAL> &hpv0=hpv[0];
  NumPtr<REAL> &hpv1=hpv[1];
  NumPtr<REAL> &hpv2=hpv[2];
  NumPtr<REAL> &hps0=hps[0];
  NumPtr<REAL> &hps1=hps[1];
  NumPtr<REAL> &hps2=hps[2];
  NumPtr<REAL> &hpss0=hpss[0];
  NumPtr<REAL> &hpss1=hpss[1];
  NumPtr<REAL> &hpss2=hpss[2];
  int iv=0;
//vertices
  for (int k=0;k<2;k++) {
    const REAL &hpvk=hpv2[k];
    const REAL &hpsk=hps2[k];
    const REAL &hpssk=hpss2[k];
    for (int j=0;j<2;j++) {
      const REAL &hpvj=hpv1[j];
      const REAL &hpsj=hps1[j];
      REAL hpvjvk=hpvj*hpvk;
      REAL hpvjsk=hpvj*hpsk;
      REAL hpvjssk=hpvj*hpssk;
      REAL hpsjvk=hpsj*hpvk;
      REAL hpsjsk=hpsj*hpsk;
      REAL hpssjvk=hpss1[j]*hpvk;
      for (int i=0;i<2;i++) {
        Tensor<2,SPACEDIM> &ggv=gg[iv++];
        const REAL &hpvi=hpv0[i];
        const REAL &hpsi=hps0[i];
        ggv[0][0]=hpss0[i]*hpvjvk;
        ggv[1][0]=hpsi*hpsjvk;
        ggv[2][0]=hpsi*hpvjsk;
        ggv[1][1]=hpvi*hpssjvk;
        ggv[2][1]=hpvi*hpsjsk;
        ggv[2][2]=hpvi*hpvjssk;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalHexahedronPolynomials::computeValuesFor(const Shape *e,
const Point<SPACEDIM> &p,NumPtr<REAL> &v) const {
  NumPtr<int> ev(e->numberVertices());
  for (int i=0;i<e->numberVertices();i++) ev[i]=e->vertexIndex(i);

  v.allocate(n_pols);
  NumPtr<NumPtr<REAL> > hpv(3);
  HierarchicalPolynomial hp;
  for (int d=0;d<3;d++) {
    hpv[d].allocate(n_pols_1d);
    hp.values(p[d],hpv[d]);
  }
  NumPtr<REAL> &hpv0=hpv[0];
  NumPtr<REAL> &hpv1=hpv[1];
  NumPtr<REAL> &hpv2=hpv[2];
  int iv=0;
  {
    for (int k=0;k<2;k++) {
      REAL hpvk=hpv2[k];
      for (int j=0;j<2;j++) {
        REAL hpvjk=hpv1[j]*hpvk;
        for (int i=0;i<2;i++) {
          v[iv++]=hpv0[i]*hpvjk;
        }
      }
    }
  }
  if (n_pols_1d<2) return;
  {
    CHECK_TEST(ev[0]<ev[2]);
    REAL hpvik=hpv0[0]*hpv2[0];
    for (int j=2;j<n_pols_1d;j++) {
      v[iv++]=hpv1[j]*hpvik;
    }
  }
  {
    REAL hpvik=hpv0[1]*hpv2[0];
    if (ev[1]<ev[3]) {
      for (int j=2;j<n_pols_1d;j++) {
        v[iv++]=hpv1[j]*hpvik;
      }
    } else {
      for (int j=2;j<n_pols_1d;j++) {
        v[iv++]=hpv1[j]*hpvik;
        hpvik=-hpvik;
      }
    }
  }
  {
    CHECK_TEST(ev[0]<ev[1]);
    REAL hpvjk=hpv1[0]*hpv2[0];
    for (int i=2;i<n_pols_1d;i++) {
      v[iv++]=hpv0[i]*hpvjk;
    }
  }
  {
    REAL hpvjk=hpv1[1]*hpv2[0];
    if (ev[2]<ev[3]) {
      for (int i=2;i<n_pols_1d;i++) {
        v[iv++]=hpv0[i]*hpvjk;
      }
    } else {
      for (int i=2;i<n_pols_1d;i++) {
        v[iv++]=hpv0[i]*hpvjk;
        hpvjk=-hpvjk;
      }
    }
  }
  {
    REAL hpvik=hpv0[0]*hpv2[1];
    if (ev[4]<ev[6]) {
      for (int j=2;j<n_pols_1d;j++) {
        v[iv++]=hpv1[j]*hpvik;
      }
    } else {
      for (int j=2;j<n_pols_1d;j++) {
        v[iv++]=hpv1[j]*hpvik;
        hpvik=-hpvik;
      }
    }
  }
  {
    REAL hpvik=hpv0[1]*hpv2[1];
    if (ev[5]<ev[7]) {
      for (int j=2;j<n_pols_1d;j++) {
        v[iv++]=hpv1[j]*hpvik;
      }
    } else {
      for (int j=2;j<n_pols_1d;j++) {
        v[iv++]=hpv1[j]*hpvik;
        hpvik=-hpvik;
      }
    }
  }
  {
    REAL hpvjk=hpv1[0]*hpv2[1];
    if (ev[4]<ev[5]) {
      for (int i=2;i<n_pols_1d;i++) {
        v[iv++]=hpv0[i]*hpvjk;
      }
    } else {
      for (int i=2;i<n_pols_1d;i++) {
        v[iv++]=hpv0[i]*hpvjk;
        hpvjk=-hpvjk;
      }
    }
  }
  {
    REAL hpvjk=hpv1[1]*hpv2[1];
    if (ev[6]<ev[7]) {
      for (int i=2;i<n_pols_1d;i++) {
        v[iv++]=hpv0[i]*hpvjk;
      }
    } else {
      for (int i=2;i<n_pols_1d;i++) {
        v[iv++]=hpv0[i]*hpvjk;
        hpvjk=-hpvjk;
      }
    }
  }
  {
    CHECK_TEST(ev[0]<ev[4]);
    REAL hpvij=hpv0[0]*hpv1[0];
    for (int k=2;k<n_pols_1d;k++) {
      v[iv++]=hpvij*hpv2[k];
    }
  }
  {
    REAL hpvij=hpv0[1]*hpv1[0];
    if (ev[1]<ev[5]) {
      for (int k=2;k<n_pols_1d;k++) {
        v[iv++]=hpvij*hpv2[k];
      }
    } else {
      for (int k=2;k<n_pols_1d;k++) {
        v[iv++]=hpvij*hpv2[k];
        hpvij=-hpvij;
      }
    }
  }
  {
    REAL hpvij=hpv0[0]*hpv1[1];
    if (ev[2]<ev[6]) {
      for (int k=2;k<n_pols_1d;k++) {
        v[iv++]=hpvij*hpv2[k];
      }
    } else {
      for (int k=2;k<n_pols_1d;k++) {
        v[iv++]=hpvij*hpv2[k];
        hpvij=-hpvij;
      }
    }
  }
  {
    REAL hpvij=hpv0[1]*hpv1[1];
    if (ev[3]<ev[7]) {
      for (int k=2;k<n_pols_1d;k++) {
        v[iv++]=hpvij*hpv2[k];
      }
    } else {
      for (int k=2;k<n_pols_1d;k++) {
        v[iv++]=hpvij*hpv2[k];
        hpvij=-hpvij;
      }
    }
  }
  { // face 0
    CHECK_TEST(ev[0]<min(ev[2],min(ev[4],ev[6])));
    const REAL &hpvi=hpv0[0];
    if (ev[2]<ev[4]) {
      for (int k=2;k<n_pols_1d;k++) {
        REAL hpvik=hpvi*hpv2[k];
        for (int j=2;j<n_pols_1d;j++) {
          v[iv++]=hpvik*hpv1[j];
        }
      }
    } else {
      for (int j=2;j<n_pols_1d;j++) {
        REAL hpvij=hpvi*hpv1[j];
        for (int k=2;k<n_pols_1d;k++) {
          v[iv++]=hpvij*hpv2[k];
        }
      }
    }
  }
  { // face 1
    REAL hpvi=hpv0[1];
    int vmin=min(min(ev[1],ev[3]),min(ev[5],ev[7]));
    if (ev[1]==vmin) {
      if (ev[3]<ev[5]) {
        for (int k=2;k<n_pols_1d;k++) {
          REAL hpvik=hpvi*hpv2[k];
          for (int j=2;j<n_pols_1d;j++) {
            v[iv++]=hpvik*hpv1[j];
          }
        }
      } else {
        for (int j=2;j<n_pols_1d;j++) {
          REAL hpvij=hpvi*hpv1[j];
          for (int k=2;k<n_pols_1d;k++) {
            v[iv++]=hpvij*hpv2[k];
          }
        }
      }
    } else if (ev[3]==vmin) {
      if (ev[1]<ev[7]) {
        for (int k=2;k<n_pols_1d;k++) {
          REAL hpvik=hpvi*hpv2[k];
          for (int j=2;j<n_pols_1d;j++) {
            v[iv++]=hpvik*hpv1[j];
            hpvik=-hpvik;
          }
        }
      } else {
        for (int j=2;j<n_pols_1d;j++) {
          REAL hpvij=hpv1[j]*hpvi;
          for (int k=2;k<n_pols_1d;k++) {
            v[iv++]=hpvij*hpv2[k];
          }
          hpvi=-hpvi;
        }
      }
    } else if (ev[5]==vmin) {
      if (ev[1]<ev[7]) {
        for (int j=2;j<n_pols_1d;j++) {
          REAL hpvij=hpv1[j]*hpvi;
          for (int k=2;k<n_pols_1d;k++) {
            v[iv++]=hpvij*hpv2[k];
            hpvij=-hpvij;
          }
        }
      } else {
        for (int k=2;k<n_pols_1d;k++) {
          REAL hpvik=hpvi*hpv2[k];
          for (int j=2;j<n_pols_1d;j++) {
            v[iv++]=hpvik*hpv1[j];
          }
          hpvi=-hpvi;
        }
      }
    } else {
      if (ev[3]<ev[5]) {
        for (int j=2;j<n_pols_1d;j++) {
          REAL hpvij=hpv1[j]*hpvi;
          for (int k=2;k<n_pols_1d;k++) {
            v[iv++]=hpvij*hpv2[k];
            hpvij=-hpvij;
          }
          hpvi=-hpvi;
        }
      } else {
        for (int k=2;k<n_pols_1d;k++) {
          REAL hpvik=hpv2[k]*hpvi;
          for (int j=2;j<n_pols_1d;j++) {
            v[iv++]=hpvik*hpv1[j];
            hpvik=-hpvik;
          }
          hpvi=-hpvi;
        }
      }
    }
  }
  { // face 2
    CHECK_TEST(ev[0]<min(ev[1],min(ev[4],ev[5])) && ev[1]<ev[4]);
    const REAL &hpvj=hpv1[0];
    for (int k=2;k<n_pols_1d;k++) {
      REAL hpvjk=hpvj*hpv2[k];
      for (int i=2;i<n_pols_1d;i++) {
        v[iv++]=hpv0[i]*hpvjk;
      }
    }
  }
  { // face 3
    REAL hpvj=hpv1[1];
    int vmin=min(min(ev[2],ev[3]),min(ev[6],ev[7]));
    if (ev[2]==vmin) {
      if (ev[3]<ev[6]) {
        for (int k=2;k<n_pols_1d;k++) {
          REAL hpvjk=hpvj*hpv2[k];
          for (int i=2;i<n_pols_1d;i++) {
            v[iv++]=hpv0[i]*hpvjk;
          }
        }
      } else {
        for (int i=2;i<n_pols_1d;i++) {
          REAL hpvij=hpv0[i]*hpvj;
          for (int k=2;k<n_pols_1d;k++) {
            v[iv++]=hpv2[k]*hpvij;
          }
        }
      }
    } else if (ev[3]==vmin) {
      if (ev[2]<ev[7]) {
        for (int k=2;k<n_pols_1d;k++) {
          REAL hpvjk=hpvj*hpv2[k];
          for (int i=2;i<n_pols_1d;i++) {
            v[iv++]=hpv0[i]*hpvjk;
            hpvjk=-hpvjk;
          }
        }
      } else {
        for (int i=2;i<n_pols_1d;i++) {
          REAL hpvij=hpv0[i]*hpvj;
          for (int k=2;k<n_pols_1d;k++) {
            v[iv++]=hpv2[k]*hpvij;
          }
          hpvj=-hpvj;
        }
      }
    } else if (ev[6]==vmin) {
      if (ev[2]<ev[7]) {
        for (int i=2;i<n_pols_1d;i++) {
          REAL hpvij=hpv0[i]*hpvj;
          for (int k=2;k<n_pols_1d;k++) {
            v[iv++]=hpv2[k]*hpvij;
            hpvij=-hpvij;
          }
        }
      } else {
        for (int k=2;k<n_pols_1d;k++) {
          REAL hpvjk=hpvj*hpv2[k];
          for (int i=2;i<n_pols_1d;i++) {
            v[iv++]=hpv0[i]*hpvjk;
          }
          hpvj=-hpvj;
        }
      }
    } else {
      if (ev[3]<ev[6]) {
        for (int i=2;i<n_pols_1d;i++) {
          REAL hpvij=hpv0[i]*hpvj;
          for (int k=2;k<n_pols_1d;k++) {
            v[iv++]=hpv2[k]*hpvij;
            hpvij=-hpvij;
          }
          hpvj=-hpvj;
        }
      } else {
        for (int k=2;k<n_pols_1d;k++) {
          REAL hpvjk=hpvj*hpv2[k];
          for (int i=2;i<n_pols_1d;i++) {
            v[iv++]=hpv0[i]*hpvjk;
            hpvjk=-hpvjk;
          }
          hpvj=-hpvj;
        }
      }
    }
  }
  { // face 4
    CHECK_TEST(ev[0]<min(ev[1],min(ev[2],ev[3])) && ev[1]<ev[2]);
    REAL hpvk=hpv2[0];
    for (int j=2;j<n_pols_1d;j++) {
      REAL hpvjk=hpv1[j]*hpvk;
      for (int i=2;i<n_pols_1d;i++) {
        v[iv++]=hpv0[i]*hpvjk;
      }
    }
  }
  { // face 5
    REAL hpvk=hpv2[1];
    int vmin=min(min(ev[4],ev[5]),min(ev[6],ev[7]));
    if (ev[4]==vmin) {
      if (ev[5]<ev[6]) {
        for (int j=2;j<n_pols_1d;j++) {
          REAL hpvjk=hpv1[j]*hpvk;
          for (int i=2;i<n_pols_1d;i++) {
            v[iv++]=hpv0[i]*hpvjk;
          }
        }
      } else {
        for (int i=2;i<n_pols_1d;i++) {
          REAL hpvik=hpv0[i]*hpvk;
          for (int j=2;j<n_pols_1d;j++) {
            v[iv++]=hpv1[j]*hpvik;
          }
        }
      }
    } else if (ev[5]==vmin) {
      if (ev[4]<ev[7]) {
        for (int j=2;j<n_pols_1d;j++) {
          REAL hpvjk=hpv1[j]*hpvk;
          for (int i=2;i<n_pols_1d;i++) {
            v[iv++]=hpv0[i]*hpvjk;
            hpvjk=-hpvjk;
          }
        }
      } else {
        for (int i=2;i<n_pols_1d;i++) {
          REAL hpvik=hpv0[i]*hpvk;
          for (int j=2;j<n_pols_1d;j++) {
            v[iv++]=hpv1[j]*hpvik;
          }
          hpvk=-hpvk;
        }
      }
    } else if (ev[6]==vmin) {
      if (ev[4]<ev[7]) {
        for (int i=2;i<n_pols_1d;i++) {
          REAL hpvik=hpv0[i]*hpvk;
          for (int j=2;j<n_pols_1d;j++) {
            v[iv++]=hpv1[j]*hpvik;
            hpvik=-hpvik;
          }
        }
      } else {
        for (int j=2;j<n_pols_1d;j++) {
          REAL hpvjk=hpv1[j]*hpvk;
          for (int i=2;i<n_pols_1d;i++) {
            v[iv++]=hpv0[i]*hpvjk;
          }
          hpvk=-hpvk;
        }
      }
    } else {
      if (ev[5]<ev[6]) {
        for (int j=2;j<n_pols_1d;j++) {
          REAL hpvjk=hpv0[j]*hpvk;
          for (int i=2;i<n_pols_1d;i++) {
            v[iv++]=hpv1[i]*hpvjk;
            hpvjk=-hpvjk;
          }
          hpvk=-hpvk;
        }
      } else {
        for (int i=2;i<n_pols_1d;i++) {
          REAL hpvik=hpv1[i]*hpvk;
          for (int j=2;j<n_pols_1d;j++) {
            v[iv++]=hpv0[j]*hpvik;
            hpvik=-hpvik;
          }
          hpvk=-hpvk;
        }
      }
    }
  }
  {
    for (int k=2;k<n_pols_1d;k++) {
      const REAL &hpvk=hpv2[k];
      for (int j=2;j<n_pols_1d;j++) {
        REAL hpvjk=hpv1[j]*hpvk;
        for (int i=2;i<n_pols_1d;i++) {
          v[iv++]=hpv0[i]*hpvjk;
        }
      }
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalHexahedronPolynomials::computeGradsFor(const Shape *e,
const Point<SPACEDIM> &p,NumPtr<Tensor<1,SPACEDIM> > &g) const {
  NumPtr<int> ev(e->numberVertices());
  for (int i=0;i<e->numberVertices();i++) ev[i]=e->vertexIndex(i);

  g.allocate(n_pols);
#ifdef INDEF
  Tensor<1,3> bogus_tensor;
  bogus_tensor[0]=HUGE_VAL;
  bogus_tensor[1]=HUGE_VAL;
  bogus_tensor[2]=HUGE_VAL;
  g.initialize(bogus_tensor);
#endif
  NumPtr<NumPtr<REAL> > hpv(3);
  NumPtr<NumPtr<REAL> > hps(3);
  HierarchicalPolynomial hp;
  for (int d=0;d<3;d++) {
    hpv[d].allocate(n_pols_1d);
    hps[d].allocate(n_pols_1d);
    hp.values(p[d],hpv[d]);
    hp.slopes(p[d],hps[d]);
  }
  NumPtr<REAL> &hpv0=hpv[0];
  NumPtr<REAL> &hpv1=hpv[1];
  NumPtr<REAL> &hpv2=hpv[2];
  NumPtr<REAL> &hps0=hps[0];
  NumPtr<REAL> &hps1=hps[1];
  NumPtr<REAL> &hps2=hps[2];
  int iv=0;
  {
    for (int k=0;k<2;k++) {
      REAL hpvk=hpv2[k];
      REAL hpsk=hps2[k];
      for (int j=0;j<2;j++) {
        REAL hpvv=hpv1[j]*hpvk;
        REAL hpsv=hps1[j]*hpvk;
        REAL hpvs=hpv1[j]*hpsk;
        for (int i=0;i<2;i++) {
          Tensor<1,SPACEDIM> &gv=g[iv++];
          REAL hpvi=hpv0[i];
          gv[0]=hps0[i]*hpvv;
          gv[1]=hpvi*hpsv;
          gv[2]=hpvi*hpvs;
        }
      }
    }
  }
  if (n_pols_1d<2) return;
  {
    CHECK_TEST(ev[0]<ev[2]);
    REAL hpvv=hpv0[0]*hpv2[0];
    REAL hpsv=hps0[0]*hpv2[0];
    REAL hpvs=hpv0[0]*hps2[0];
    for (int j=2;j<n_pols_1d;j++) {
      Tensor<1,SPACEDIM> &gv=g[iv++];
      REAL hpvj=hpv1[j];
      gv[0]=hpsv*hpvj;
      gv[1]=hpvv*hps1[j];
      gv[2]=hpvs*hpvj;
    }
  }
  {
    REAL hpvv=hpv0[1]*hpv2[0];
    REAL hpsv=hps0[1]*hpv2[0];
    REAL hpvs=hpv0[1]*hps2[0];
    if (ev[1]<ev[3]) {
      for (int j=2;j<n_pols_1d;j++) {
        Tensor<1,SPACEDIM> &gv=g[iv++];
        REAL hpvj=hpv1[j];
        gv[0]=hpsv*hpvj;
        gv[1]=hpvv*hps1[j];
        gv[2]=hpvs*hpvj;
      }
    } else {
      for (int j=2;j<n_pols_1d;j++) {
        Tensor<1,SPACEDIM> &gv=g[iv++];
        REAL hpvj=hpv1[j];
        gv[0]=hpsv*hpvj;
        gv[1]=hpvv*hps1[j];
        gv[2]=hpvs*hpvj;
        hpvv=-hpvv;
        hpsv=-hpsv;
        hpvs=-hpvs;
      }
    }
  }
  {
    CHECK_TEST(ev[0]<ev[1]);
    REAL hpvv=hpv1[0]*hpv2[0];
    REAL hpsv=hps1[0]*hpv2[0];
    REAL hpvs=hpv1[0]*hps2[0];
    for (int i=2;i<n_pols_1d;i++) {
      Tensor<1,SPACEDIM> &gv=g[iv++];
      REAL hpvi=hpv0[i];
      gv[0]=hps0[i]*hpvv;
      gv[1]=hpvi*hpsv;
      gv[2]=hpvi*hpvs;
    }
  }
  {
    REAL hpvv=hpv1[1]*hpv2[0];
    REAL hpsv=hps1[1]*hpv2[0];
    REAL hpvs=hpv1[1]*hps2[0];
    if (ev[2]<ev[3]) {
      for (int i=2;i<n_pols_1d;i++) {
        Tensor<1,SPACEDIM> &gv=g[iv++];
        REAL hpvi=hpv0[i];
        gv[0]=hps0[i]*hpvv;
        gv[1]=hpvi*hpsv;
        gv[2]=hpvi*hpvs;
      }
    } else {
      for (int i=2;i<n_pols_1d;i++) {
        Tensor<1,SPACEDIM> &gv=g[iv++];
        REAL hpvi=hpv0[i];
        gv[0]=hps0[i]*hpvv;
        gv[1]=hpvi*hpsv;
        gv[2]=hpvi*hpvs;
        hpvv=-hpvv;
        hpsv=-hpsv;
        hpvs=-hpvs;
      }
    }
  }
  {
    REAL hpvv=hpv0[0]*hpv2[1];
    REAL hpsv=hps0[0]*hpv2[1];
    REAL hpvs=hpv0[0]*hps2[1];
    if (ev[4]<ev[6]) {
      for (int j=2;j<n_pols_1d;j++) {
        Tensor<1,SPACEDIM> &gv=g[iv++];
        REAL hpvj=hpv1[j];
        gv[0]=hpvj*hpsv;
        gv[1]=hps1[j]*hpvv;
        gv[2]=hpvj*hpvs;
      }
    } else {
      for (int j=2;j<n_pols_1d;j++) {
        Tensor<1,SPACEDIM> &gv=g[iv++];
        REAL hpvj=hpv1[j];
        gv[0]=hpvj*hpsv;
        gv[1]=hps1[j]*hpvv;
        gv[2]=hpvj*hpvs;
        hpvv=-hpvv;
        hpsv=-hpsv;
        hpvs=-hpvs;
      }
    }
  }
  {
    REAL hpvv=hpv0[1]*hpv2[1];
    REAL hpsv=hps0[1]*hpv2[1];
    REAL hpvs=hpv0[1]*hps2[1];
    if (ev[5]<ev[7]) {
      for (int j=2;j<n_pols_1d;j++) {
        Tensor<1,SPACEDIM> &gv=g[iv++];
        REAL hpvj=hpv1[j];
        gv[0]=hpvj*hpsv;
        gv[1]=hps1[j]*hpvv;
        gv[2]=hpvj*hpvs;
      }
    } else {
      for (int j=2;j<n_pols_1d;j++) {
        Tensor<1,SPACEDIM> &gv=g[iv++];
        REAL hpvj=hpv1[j];
        gv[0]=hpvj*hpsv;
        gv[1]=hps1[j]*hpvv;
        gv[2]=hpvj*hpvs;
        hpvv=-hpvv;
        hpsv=-hpsv;
        hpvs=-hpvs;
      }
    }
  }
  {
    REAL hpvv=hpv1[0]*hpv2[1];
    REAL hpsv=hps1[0]*hpv2[1];
    REAL hpvs=hpv1[0]*hps2[1];
    if (ev[4]<ev[5]) {
      for (int i=2;i<n_pols_1d;i++) {
        Tensor<1,SPACEDIM> &gv=g[iv++];
        REAL hpvi=hpv0[i];
        gv[0]=hps0[i]*hpvv;
        gv[1]=hpvi*hpsv;
        gv[2]=hpvi*hpvs;
      }
    } else {
      for (int i=2;i<n_pols_1d;i++) {
        Tensor<1,SPACEDIM> &gv=g[iv++];
        REAL hpvi=hpv0[i];
        gv[0]=hps0[i]*hpvv;
        gv[1]=hpvi*hpsv;
        gv[2]=hpvi*hpvs;
        hpvv=-hpvv;
        hpsv=-hpsv;
        hpvs=-hpvs;
      }
    }
  }
  {
    REAL hpvv=hpv1[1]*hpv2[1];
    REAL hpsv=hps1[1]*hpv2[1];
    REAL hpvs=hpv1[1]*hps2[1];
    if (ev[6]<ev[7]) {
      for (int i=2;i<n_pols_1d;i++) {
        Tensor<1,SPACEDIM> &gv=g[iv++];
        REAL hpvi=hpv0[i];
        gv[0]=hps0[i]*hpvv;
        gv[1]=hpvi*hpsv;
        gv[2]=hpvi*hpvs;
      }
    } else {
      for (int i=2;i<n_pols_1d;i++) {
        Tensor<1,SPACEDIM> &gv=g[iv++];
        REAL hpvi=hpv0[i];
        gv[0]=hps0[i]*hpvv;
        gv[1]=hpvi*hpsv;
        gv[2]=hpvi*hpvs;
        hpvv=-hpvv;
        hpsv=-hpsv;
        hpvs=-hpvs;
      }
    }
  }
  {
    CHECK_TEST(ev[0]<ev[4]);
    REAL hpvv=hpv0[0]*hpv1[0];
    REAL hpsv=hps0[0]*hpv1[0];
    REAL hpvs=hpv0[0]*hps1[0];
    for (int k=2;k<n_pols_1d;k++) {
      Tensor<1,SPACEDIM> &gv=g[iv++];
      REAL hpvk=hpv2[k];
      gv[0]=hpsv*hpvk;
      gv[1]=hpvs*hpvk;
      gv[2]=hpvv*hps2[k];
    }
  }
  {
    REAL hpvv=hpv0[1]*hpv1[0];
    REAL hpsv=hps0[1]*hpv1[0];
    REAL hpvs=hpv0[1]*hps1[0];
    if (ev[1]<ev[5]) {
      for (int k=2;k<n_pols_1d;k++) {
        Tensor<1,SPACEDIM> &gv=g[iv++];
        REAL hpvk=hpv2[k];
        gv[0]=hpsv*hpvk;
        gv[1]=hpvs*hpvk;
        gv[2]=hpvv*hps2[k];
      }
    } else {
      for (int k=2;k<n_pols_1d;k++) {
        Tensor<1,SPACEDIM> &gv=g[iv++];
        REAL hpvk=hpv2[k];
        gv[0]=hpsv*hpvk;
        gv[1]=hpvs*hpvk;
        gv[2]=hpvv*hps2[k];
        hpvv=-hpvv;
        hpsv=-hpsv;
        hpvs=-hpvs;
      }
    }
  }
  {
    REAL hpvv=hpv0[0]*hpv1[1];
    REAL hpsv=hps0[0]*hpv1[1];
    REAL hpvs=hpv0[0]*hps1[1];
    if (ev[2]<ev[6]) {
      for (int k=2;k<n_pols_1d;k++) {
        Tensor<1,SPACEDIM> &gv=g[iv++];
        REAL hpvk=hpv2[k];
        gv[0]=hpsv*hpvk;
        gv[1]=hpvs*hpvk;
        gv[2]=hpvv*hps2[k];
      }
    } else {
      for (int k=2;k<n_pols_1d;k++) {
        Tensor<1,SPACEDIM> &gv=g[iv++];
        REAL hpvk=hpv2[k];
        gv[0]=hpsv*hpvk;
        gv[1]=hpvs*hpvk;
        gv[2]=hpvv*hps2[k];
        hpvv=-hpvv;
        hpsv=-hpsv;
        hpvs=-hpvs;
      }
    }
  }
  {
    REAL hpvv=hpv0[1]*hpv1[1];
    REAL hpsv=hps0[1]*hpv1[1];
    REAL hpvs=hpv0[1]*hps1[1];
    if (ev[3]<ev[7]) {
      for (int k=2;k<n_pols_1d;k++) {
        Tensor<1,SPACEDIM> &gv=g[iv++];
        REAL hpvk=hpv2[k];
        gv[0]=hpsv*hpvk;
        gv[1]=hpvs*hpvk;
        gv[2]=hpvv*hps2[k];
      }
    } else {
      for (int k=2;k<n_pols_1d;k++) {
        Tensor<1,SPACEDIM> &gv=g[iv++];
        REAL hpvk=hpv2[k];
        gv[0]=hpsv*hpvk;
        gv[1]=hpvs*hpvk;
        gv[2]=hpvv*hps2[k];
        hpvv=-hpvv;
        hpsv=-hpsv;
        hpvs=-hpvs;
      }
    }
  }
  {
    CHECK_TEST(ev[0]<min(ev[2],min(ev[4],ev[6])));
    if (ev[2]<ev[4]) {
      for (int k=2;k<n_pols_1d;k++) {
        REAL hpvv=hpv0[0]*hpv2[k];
        REAL hpsv=hps0[0]*hpv2[k];
        REAL hpvs=hpv0[0]*hps2[k];
        for (int j=2;j<n_pols_1d;j++) {
          Tensor<1,SPACEDIM> &gv=g[iv++];
          REAL hpvj=hpv1[j];
          gv[0]=hpsv*hpvj;
          gv[1]=hpvv*hps1[j];
          gv[2]=hpvs*hpvj;
        }
      }
    } else {
      for (int j=2;j<n_pols_1d;j++) {
        REAL hpvv=hpv0[0]*hpv1[j];
        REAL hpsv=hps0[0]*hpv1[j];
        REAL hpvs=hpv0[0]*hps1[j];
        for (int k=2;k<n_pols_1d;k++) {
          Tensor<1,SPACEDIM> &gv=g[iv++];
          REAL hpvk=hpv2[k];
          gv[0]=hpsv*hpvk;
          gv[1]=hpvs*hpvk;
          gv[2]=hpvv*hps2[k];
        }
      }
    }
  }
  {
    int vmin=min(min(ev[1],ev[3]),min(ev[5],ev[7]));
    if (ev[1]==vmin) {
      if (ev[3]<ev[5]) {
        for (int k=2;k<n_pols_1d;k++) {
          REAL hpvv=hpv0[1]*hpv2[k];
          REAL hpsv=hps0[1]*hpv2[k];
          REAL hpvs=hpv0[1]*hps2[k];
          for (int j=2;j<n_pols_1d;j++) {
            Tensor<1,SPACEDIM> &gv=g[iv++];
            REAL hpvj=hpv1[j];
            gv[0]=hpsv*hpvj;
            gv[1]=hpvv*hps1[j];
            gv[2]=hpvs*hpvj;
          }
        }
      } else {
        for (int j=2;j<n_pols_1d;j++) {
          REAL hpvv=hpv0[1]*hpv1[j];
          REAL hpsv=hps0[1]*hpv1[j];
          REAL hpvs=hpv0[1]*hps1[j];
          for (int k=2;k<n_pols_1d;k++) {
            Tensor<1,SPACEDIM> &gv=g[iv++];
            REAL hpvk=hpv2[k];
            gv[0]=hpsv*hpvk;
            gv[1]=hpvs*hpvk;
            gv[2]=hpvv*hps2[k];
          }
        }
      }
    } else if (ev[3]==vmin) {
      if (ev[1]<ev[7]) {
        for (int k=2;k<n_pols_1d;k++) {
          REAL hpvv=hpv0[1]*hpv2[k];
          REAL hpsv=hps0[1]*hpv2[k];
          REAL hpvs=hpv0[1]*hps2[k];
          for (int j=2;j<n_pols_1d;j++) {
            Tensor<1,SPACEDIM> &gv=g[iv++];
            REAL hpvj=hpv1[j];
            gv[0]=hpsv*hpvj;
            gv[1]=hpvv*hps1[j];
            gv[2]=hpvs*hpvj;
            hpvv=-hpvv;
            hpsv=-hpsv;
            hpvs=-hpvs;
          }
        }
      } else {
        REAL hpvi=hpv0[1];
        REAL hpsi=hps0[1];
        for (int j=2;j<n_pols_1d;j++) {
          REAL hpvv=hpvi*hpv1[j];
          REAL hpsv=hpsi*hpv1[j];
          REAL hpvs=hpvi*hps1[j];
          for (int k=2;k<n_pols_1d;k++) {
            Tensor<1,SPACEDIM> &gv=g[iv++];
            REAL hpvk=hpv2[k];
            gv[0]=hpsv*hpvk;
            gv[1]=hpvs*hpvk;
            gv[2]=hpvv*hps2[k];
          }
          hpvi=-hpvi;
          hpsi=-hpsi;
        }
      }
    } else if (ev[5]==vmin) {
      if (ev[1]<ev[7]) {
        REAL hpvi=hpv0[1];
        REAL hpsi=hps0[1];
        for (int j=2;j<n_pols_1d;j++) {
          REAL hpvv=hpvi*hpv1[j];
          REAL hpsv=hpsi*hpv1[j];
          REAL hpvs=hpvi*hps1[j];
          for (int k=2;k<n_pols_1d;k++) {
            Tensor<1,SPACEDIM> &gv=g[iv++];
            REAL hpvk=hpv2[k];
            gv[0]=hpsv*hpvk;
            gv[1]=hpvs*hpvk;
            gv[2]=hpvv*hps2[k];
            hpvv=-hpvv;
            hpsv=-hpsv;
            hpvs=-hpvs;
          }
        }
      } else {
        REAL hpvi=hpv0[1];
        REAL hpsi=hps0[1];
        for (int k=2;k<n_pols_1d;k++) {
          REAL hpvv=hpvi*hpv2[k];
          REAL hpsv=hpsi*hpv2[k];
          REAL hpvs=hpvi*hps2[k];
          for (int j=2;j<n_pols_1d;j++) {
            Tensor<1,SPACEDIM> &gv=g[iv++];
            REAL hpvj=hpv1[j];
            gv[0]=hpsv*hpvj;
            gv[1]=hpvv*hps1[j];
            gv[2]=hpvs*hpvj;
          }
          hpvi=-hpvi;
          hpsi=-hpsi;
        }
      }
    } else {
      if (ev[3]<ev[5]) {
        REAL hpvi=hpv0[1];
        REAL hpsi=hps0[1];
        for (int j=2;j<n_pols_1d;j++) {
          REAL hpvv=hpvi*hpv1[j];
          REAL hpsv=hpsi*hpv1[j];
          REAL hpvs=hpvi*hps1[j];
          for (int k=2;k<n_pols_1d;k++) {
            Tensor<1,SPACEDIM> &gv=g[iv++];
            REAL hpvk=hpv2[k];
            gv[0]=hpsv*hpvk;
            gv[1]=hpvs*hpvk;
            gv[2]=hpvv*hps2[k];
            hpvv=-hpvv;
            hpsv=-hpsv;
            hpvs=-hpvs;
          }
          hpvi=-hpvi;
          hpsi=-hpsi;
        }
      } else {
        REAL hpvi=hpv0[1];
        REAL hpsi=hps0[1];
        for (int k=2;k<n_pols_1d;k++) {
          REAL hpvv=hpvi*hpv2[k];
          REAL hpsv=hpsi*hpv2[k];
          REAL hpvs=hpvi*hps2[k];
          for (int j=2;j<n_pols_1d;j++) {
            Tensor<1,SPACEDIM> &gv=g[iv++];
            REAL hpvj=hpv1[j];
            gv[0]=hpsv*hpvj;
            gv[1]=hpvv*hps1[j];
            gv[2]=hpvs*hpvj;
            hpvv=-hpvv;
            hpsv=-hpsv;
            hpvs=-hpvs;
          }
          hpvi=-hpvi;
          hpsi=-hpsi;
        }
      }
    }
  }
  {
    CHECK_TEST(ev[0]<min(ev[1],min(ev[4],ev[5])) && ev[1]<ev[4]);
    const REAL &hpvj=hpv1[0];
    for (int k=2;k<n_pols_1d;k++) {
      REAL hpvv=hpv1[0]*hpv2[k];
      REAL hpsv=hps1[0]*hpv2[k];
      REAL hpvs=hpv1[0]*hps2[k];
      for (int i=2;i<n_pols_1d;i++) {
        Tensor<1,SPACEDIM> &gv=g[iv++];
        REAL hpvi=hpv0[i];
        gv[0]=hps0[i]*hpvv;
        gv[1]=hpvi*hpsv;
        gv[2]=hpvi*hpvs;
      }
    }
  }
  {
    REAL hpvj=hpv1[1];
    int vmin=min(min(ev[2],ev[3]),min(ev[6],ev[7]));
    if (ev[2]==vmin) {
      if (ev[3]<ev[6]) {
        for (int k=2;k<n_pols_1d;k++) {
          REAL hpvv=hpv1[1]*hpv2[k];
          REAL hpsv=hps1[1]*hpv2[k];
          REAL hpvs=hpv1[1]*hps2[k];
          for (int i=2;i<n_pols_1d;i++) {
            Tensor<1,SPACEDIM> &gv=g[iv++];
            REAL hpvi=hpv0[i];
            gv[0]=hps0[i]*hpvv;
            gv[1]=hpvi*hpsv;
            gv[2]=hpvi*hpvs;
          }
        }
      } else {
        for (int i=2;i<n_pols_1d;i++) {
          REAL hpvv=hpv0[i]*hpv1[1];
          REAL hpsv=hps0[i]*hpv1[1];
          REAL hpvs=hpv0[i]*hps1[1];
          for (int k=2;k<n_pols_1d;k++) {
            Tensor<1,SPACEDIM> &gv=g[iv++];
            REAL hpvk=hpv2[k];
            gv[0]=hpsv*hpvk;
            gv[1]=hpvs*hpvk;
            gv[2]=hpvv*hps2[k];
          }
        }
      }
    } else if (ev[3]==vmin) {
      if (ev[2]<ev[7]) {
        for (int k=2;k<n_pols_1d;k++) {
          REAL hpvv=hpv1[1]*hpv2[k];
          REAL hpsv=hps1[1]*hpv2[k];
          REAL hpvs=hpv1[1]*hps2[k];
          for (int i=2;i<n_pols_1d;i++) {
            Tensor<1,SPACEDIM> &gv=g[iv++];
            REAL hpvi=hpv0[i];
            gv[0]=hps0[i]*hpvv;
            gv[1]=hpvi*hpsv;
            gv[2]=hpvi*hpvs;
            hpvv=-hpvv;
            hpsv=-hpsv;
            hpvs=-hpvs;
          }
        }
      } else {
        REAL hpvj=hpv1[1];
        REAL hpsj=hps1[1];
        for (int i=2;i<n_pols_1d;i++) {
          REAL hpvv=hpv0[i]*hpvj;
          REAL hpsv=hps0[i]*hpvj;
          REAL hpvs=hpv0[i]*hpsj;
          for (int k=2;k<n_pols_1d;k++) {
            Tensor<1,SPACEDIM> &gv=g[iv++];
            REAL hpvk=hpv2[k];
            gv[0]=hpsv*hpvk;
            gv[1]=hpvs*hpvk;
            gv[2]=hpvv*hps2[k];
          }
          hpvj=-hpvj;
          hpsj=-hpsj;
        }
      }
    } else if (ev[6]==vmin) {
      if (ev[2]<ev[7]) {
        for (int i=2;i<n_pols_1d;i++) {
          REAL hpvv=hpv0[i]*hpv1[1];
          REAL hpsv=hps0[i]*hpv1[1];
          REAL hpvs=hpv0[i]*hps1[1];
          for (int k=2;k<n_pols_1d;k++) {
            Tensor<1,SPACEDIM> &gv=g[iv++];
            REAL hpvk=hpv2[k];
            gv[0]=hpsv*hpvk;
            gv[1]=hpvs*hpvk;
            gv[2]=hpvv*hps2[k];
            hpvv=-hpvv;
            hpsv=-hpsv;
            hpvs=-hpvs;
          }
        }
      } else {
        REAL hpvj=hpv1[1];
        REAL hpsj=hps1[1];
        for (int k=2;k<n_pols_1d;k++) {
          REAL hpvv=hpvj*hpv2[k];
          REAL hpsv=hpsj*hpv2[k];
          REAL hpvs=hpvj*hps2[k];
          for (int i=2;i<n_pols_1d;i++) {
            Tensor<1,SPACEDIM> &gv=g[iv++];
            REAL hpvi=hpv0[i];
            gv[0]=hps0[i]*hpvv;
            gv[1]=hpvi*hpsv;
            gv[2]=hpvi*hpvs;
          }
          hpvj=-hpvj;
          hpsj=-hpsj;
        }
      }
    } else {
      if (ev[3]<ev[6]) {
        REAL hpvj=hpv1[1];
        REAL hpsj=hps1[1];
        for (int i=2;i<n_pols_1d;i++) {
          REAL hpvv=hpv0[i]*hpvj;
          REAL hpsv=hps0[i]*hpvj;
          REAL hpvs=hpv0[i]*hpsj;
          for (int k=2;k<n_pols_1d;k++) {
            Tensor<1,SPACEDIM> &gv=g[iv++];
            REAL hpvk=hpv2[k];
            gv[0]=hpsv*hpvk;
            gv[1]=hpvs*hpvk;
            gv[2]=hpvv*hps2[k];
            hpvv=-hpvv;
            hpsv=-hpsv;
            hpvs=-hpvs;
          }
          hpvj=-hpvj;
          hpsj=-hpsj;
        }
      } else {
        REAL hpvj=hpv1[1];
        REAL hpsj=hps1[1];
        for (int k=2;k<n_pols_1d;k++) {
          REAL hpvv=hpvj*hpv2[k];
          REAL hpsv=hpsj*hpv2[k];
          REAL hpvs=hpvj*hps2[k];
          for (int i=2;i<n_pols_1d;i++) {
            Tensor<1,SPACEDIM> &gv=g[iv++];
            REAL hpvi=hpv0[i];
            gv[0]=hps0[i]*hpvv;
            gv[1]=hpvi*hpsv;
            gv[2]=hpvi*hpvs;
            hpvv=-hpvv;
            hpsv=-hpsv;
            hpvs=-hpvs;
          }
          hpvj=-hpvj;
          hpsj=-hpsj;
        }
      }
    }
  }
  {
    CHECK_TEST(ev[0]<min(ev[1],min(ev[2],ev[3])) && ev[1]<ev[2]);
    for (int j=2;j<n_pols_1d;j++) {
      REAL hpvv=hpv1[j]*hpv2[0];
      REAL hpsv=hps1[j]*hpv2[0];
      REAL hpvs=hpv1[j]*hps2[0];
      for (int i=2;i<n_pols_1d;i++) {
        Tensor<1,SPACEDIM> &gv=g[iv++];
        REAL hpvi=hpv0[i];
        gv[0]=hps0[i]*hpvv;
        gv[1]=hpvi*hpsv;
        gv[2]=hpvi*hpvs;
      }
    }
  }
  {
    REAL hpvk=hpv2[1];
    int vmin=min(min(ev[4],ev[5]),min(ev[6],ev[7]));
    if (ev[4]==vmin) {
      if (ev[5]<ev[6]) {
        for (int j=2;j<n_pols_1d;j++) {
          REAL hpvv=hpv1[j]*hpv2[1];
          REAL hpsv=hps1[j]*hpv2[1];
          REAL hpvs=hpv1[j]*hps2[1];
          for (int i=2;i<n_pols_1d;i++) {
            Tensor<1,SPACEDIM> &gv=g[iv++];
            REAL hpvi=hpv0[i];
            gv[0]=hps0[i]*hpvv;
            gv[1]=hpvi*hpsv;
            gv[2]=hpvi*hpvs;
          }
        }
      } else {
        for (int i=2;i<n_pols_1d;i++) {
          REAL hpvv=hpv0[i]*hpv2[1];
          REAL hpsv=hps0[i]*hpv2[1];
          REAL hpvs=hpv0[i]*hps2[1];
          for (int j=2;j<n_pols_1d;j++) {
            Tensor<1,SPACEDIM> &gv=g[iv++];
            REAL hpvj=hpv1[j];
            gv[0]=hpvj*hpsv;
            gv[1]=hps1[j]*hpvv;
            gv[2]=hpvj*hpvs;
          }
        }
      }
    } else if (ev[5]==vmin) {
      if (ev[4]<ev[7]) {
        for (int j=2;j<n_pols_1d;j++) {
          REAL hpvv=hpv1[j]*hpv2[1];
          REAL hpsv=hps1[j]*hpv2[1];
          REAL hpvs=hpv1[j]*hps2[1];
          for (int i=2;i<n_pols_1d;i++) {
            Tensor<1,SPACEDIM> &gv=g[iv++];
            REAL hpvi=hpv0[i];
            gv[0]=hps0[i]*hpvv;
            gv[1]=hpvi*hpsv;
            gv[2]=hpvi*hpvs;
            hpvv=-hpvv;
            hpsv=-hpsv;
            hpvs=-hpvs;
          }
        }
      } else {
        REAL hpvk=hpv2[1];
        REAL hpsk=hps2[1];
        for (int i=2;i<n_pols_1d;i++) {
          REAL hpvv=hpv0[i]*hpvk;
          REAL hpsv=hps0[i]*hpvk;
          REAL hpvs=hpv0[i]*hpsk;
          for (int j=2;j<n_pols_1d;j++) {
            Tensor<1,SPACEDIM> &gv=g[iv++];
            REAL hpvj=hpv1[j];
            gv[0]=hpvj*hpsv;
            gv[1]=hps1[j]*hpvv;
            gv[2]=hpvj*hpvs;
          }
          hpvk=-hpvk;
          hpsk=-hpsk;
        }
      }
    } else if (ev[6]==vmin) {
      if (ev[4]<ev[7]) {
        for (int i=2;i<n_pols_1d;i++) {
          REAL hpvv=hpv0[i]*hpv2[1];
          REAL hpsv=hps0[i]*hpv2[1];
          REAL hpvs=hpv0[i]*hps2[1];
          for (int j=2;j<n_pols_1d;j++) {
            Tensor<1,SPACEDIM> &gv=g[iv++];
            REAL hpvj=hpv1[j];
            gv[0]=hpvj*hpsv;
            gv[1]=hps1[j]*hpvv;
            gv[2]=hpvj*hpvs;
            hpvv=-hpvv;
            hpsv=-hpsv;
            hpvs=-hpvs;
          }
        }
      } else {
        REAL hpvk=hpv2[1];
        REAL hpsk=hps2[1];
        for (int j=2;j<n_pols_1d;j++) {
          REAL hpvv=hpv1[j]*hpvk;
          REAL hpsv=hps1[j]*hpvk;
          REAL hpvs=hpv1[j]*hpsk;
          for (int i=2;i<n_pols_1d;i++) {
            Tensor<1,SPACEDIM> &gv=g[iv++];
            REAL hpvi=hpv0[i];
            gv[0]=hps0[i]*hpvv;
            gv[1]=hpvi*hpsv;
            gv[2]=hpvi*hpvs;
          }
          hpvk=-hpvk;
          hpsk=-hpsk;
        }
      }
    } else {
      if (ev[5]<ev[6]) {
        REAL hpvk=hpv2[1];
        REAL hpsk=hps2[1];
        for (int j=2;j<n_pols_1d;j++) {
          REAL hpvv=hpv0[j]*hpvk;
          REAL hpsv=hps0[j]*hpvk;
          REAL hpvs=hpv0[j]*hpsk;
          for (int i=2;i<n_pols_1d;i++) {
            Tensor<1,SPACEDIM> &gv=g[iv++];
            REAL hpvi=hpv1[i];
            gv[1]=hps1[i]*hpvv;
            gv[0]=hpvi*hpsv;
            gv[2]=hpvi*hpvs;
            hpvv=-hpvv;
            hpsv=-hpsv;
            hpvs=-hpvs;
          }
          hpvk=-hpvk;
          hpsk=-hpsk;
        }
      } else {
        REAL hpvk=hpv2[1];
        REAL hpsk=hps2[1];
        for (int i=2;i<n_pols_1d;i++) {
          REAL hpvv=hpv1[i]*hpvk;
          REAL hpsv=hps1[i]*hpvk;
          REAL hpvs=hpv1[i]*hpsk;
          for (int j=2;j<n_pols_1d;j++) {
            Tensor<1,SPACEDIM> &gv=g[iv++];
            REAL hpvj=hpv0[j];
            gv[1]=hpvj*hpsv;
            gv[0]=hps0[j]*hpvv;
            gv[2]=hpvj*hpvs;
            hpvv=-hpvv;
            hpsv=-hpsv;
            hpvs=-hpvs;
          }
          hpvk=-hpvk;
          hpsk=-hpsk;
        }
      }
    }
  }
  {
    for (int k=2;k<n_pols_1d;k++) {
      const REAL &hpvk=hpv2[k];
      const REAL &hpsk=hps2[k];
      for (int j=2;j<n_pols_1d;j++) {
        REAL hpvv=hpv1[j]*hpvk;
        REAL hpsv=hps1[j]*hpvk;
        REAL hpvs=hpv1[j]*hpsk;
        for (int i=2;i<n_pols_1d;i++) {
          Tensor<1,SPACEDIM> &gv=g[iv++];
          REAL hpvi=hpv0[i];
          gv[0]=hps0[i]*hpvv;
          gv[1]=hpvi*hpsv;
          gv[2]=hpvi*hpvs;
        }
      }
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalHexahedronPolynomials::computeGradGradsFor(
const Shape *e,const Point<SPACEDIM> &p,
NumPtr<Tensor<2,SPACEDIM> > &gg) const {
  NumPtr<int> ev(e->numberVertices());
  for (int i=0;i<e->numberVertices();i++) ev[i]=e->vertexIndex(i);

  gg.allocate(n_pols);
#ifdef INDEF
  Tensor<2,SPACEDIM> bogus_tensor;
  bogus_tensor.initialize(HUGE_VAL);
  gg.initialize(bogus_tensor);
#endif
  NumPtr<NumPtr<REAL> > hpv(3);
  NumPtr<NumPtr<REAL> > hps(3);
  NumPtr<NumPtr<REAL> > hpss(3);
  HierarchicalPolynomial hp;
  for (int d=0;d<3;d++) {
    hpv[d].allocate(n_pols_1d);
    hps[d].allocate(n_pols_1d);
    hpss[d].allocate(n_pols_1d);
    hp.values(p[d],hpv[d]);
    hp.slopes(p[d],hps[d]);
    hp.slope2s(p[d],hpss[d]);
  }
  NumPtr<REAL> &hpv0=hpv[0];
  NumPtr<REAL> &hpv1=hpv[1];
  NumPtr<REAL> &hpv2=hpv[2];
  NumPtr<REAL> &hps0=hps[0];
  NumPtr<REAL> &hps1=hps[1];
  NumPtr<REAL> &hps2=hps[2];
  NumPtr<REAL> &hpss0=hpss[0];
  NumPtr<REAL> &hpss1=hpss[1];
  NumPtr<REAL> &hpss2=hpss[2];
  int iv=0;
  {
    for (int k=0;k<2;k++) {
      REAL hpvk=hpv2[k];
      REAL hpsk=hps2[k];
      REAL hpssk=hpss2[k];
      for (int j=0;j<2;j++) {
        REAL hpvv=hpv1[j]*hpvk;
        REAL hpsv=hps1[j]*hpvk;
        REAL hpvs=hpv1[j]*hpsk;
        REAL hpss=hps1[j]*hpsk;
        REAL hpssv=hpss1[j]*hpvk;
        REAL hpvss=hpv1[j]*hpssk;
        for (int i=0;i<2;i++) {
          Tensor<2,3> &ggv=gg[iv++];
          REAL hpvi=hpv0[i];
          REAL hpsi=hps0[i];
          ggv[0][0]=hpss0[i]*hpvv;
          ggv[1][0]=hpsi*hpsv;
          ggv[2][0]=hpsi*hpvs;
          ggv[1][1]=hpvi*hpssv;
          ggv[2][1]=hpvi*hpss;
          ggv[2][2]=hpvi*hpvss;
          ggv[0][1]=ggv[1][0];
          ggv[0][2]=ggv[2][0];
          ggv[1][2]=ggv[2][1];
        }
      }
    }
  }
  if (n_pols_1d<2) return;
  {
    CHECK_TEST(ev[0]<ev[2]);
    REAL hpvv=hpv0[0]*hpv2[0];
    REAL hpsv=hps0[0]*hpv2[0];
    REAL hpvs=hpv0[0]*hps2[0];
    REAL hpss=hps0[0]*hps2[0];
    REAL hpssv=hpss0[0]*hpv2[0];
    REAL hpvss=hpv0[0]*hpss2[0];
    for (int j=2;j<n_pols_1d;j++) {
      Tensor<2,3> &ggv=gg[iv++];
      REAL hpvj=hpv1[j];
      REAL hpsj=hps1[j];
      ggv[0][0]=hpssv*hpvj;
      ggv[1][0]=hpsv*hpsj;
      ggv[2][0]=hpss*hpvj;
      ggv[1][1]=hpvv*hpss1[j];
      ggv[2][1]=hpvs*hpsj;
      ggv[2][2]=hpvss*hpvj;
      ggv[0][1]=ggv[1][0];
      ggv[0][2]=ggv[2][0];
      ggv[1][2]=ggv[2][1];
    }
  }
  {
    REAL hpvv=hpv0[1]*hpv2[0];
    REAL hpsv=hps0[1]*hpv2[0];
    REAL hpvs=hpv0[1]*hps2[0];
    REAL hpss=hps0[1]*hps2[0];
    REAL hpssv=hpss0[1]*hpv2[0];
    REAL hpvss=hpv0[1]*hpss2[0];
    if (ev[1]<ev[3]) {
      for (int j=2;j<n_pols_1d;j++) {
        Tensor<2,3> &ggv=gg[iv++];
        REAL hpvj=hpv1[j];
        REAL hpsj=hps1[j];
        ggv[0][0]=hpssv*hpvj;
        ggv[1][0]=hpsv*hpsj;
        ggv[2][0]=hpss*hpvj;
        ggv[1][1]=hpvv*hpss1[j];
        ggv[2][1]=hpvs*hpsj;
        ggv[2][2]=hpvss*hpvj;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    } else {
      for (int j=2;j<n_pols_1d;j++) {
        Tensor<2,3> &ggv=gg[iv++];
        REAL hpvj=hpv1[j];
        REAL hpsj=hps1[j];
        ggv[0][0]=hpssv*hpvj;
        ggv[1][0]=hpsv*hpsj;
        ggv[2][0]=hpss*hpvj;
        ggv[1][1]=hpvv*hpss1[j];
        ggv[2][1]=hpvs*hpsj;
        ggv[2][2]=hpvss*hpvj;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
        hpvv=-hpvv;
        hpsv=-hpsv;
        hpvs=-hpvs;
        hpss=-hpss;
        hpssv=-hpssv;
        hpvss=-hpvss;
      }
    }
  }
  {
    CHECK_TEST(ev[0]<ev[1]);
    REAL hpvv=hpv1[0]*hpv2[0];
    REAL hpsv=hps1[0]*hpv2[0];
    REAL hpvs=hpv1[0]*hps2[0];
    REAL hpss=hps1[0]*hps2[0];
    REAL hpssv=hpss1[0]*hpv2[0];
    REAL hpvss=hpv1[0]*hpss2[0];
    for (int i=2;i<n_pols_1d;i++) {
      Tensor<2,3> &ggv=gg[iv++];
      REAL hpvi=hpv0[i];
      REAL hpsi=hps0[i];
      ggv[0][0]=hpss0[i]*hpvv;
      ggv[1][0]=hpsi*hpsv;
      ggv[2][0]=hpsi*hpvs;
      ggv[1][1]=hpvi*hpssv;
      ggv[2][1]=hpvi*hpss;
      ggv[2][2]=hpvi*hpvss;
      ggv[0][1]=ggv[1][0];
      ggv[0][2]=ggv[2][0];
      ggv[1][2]=ggv[2][1];
    }
  }
  {
    REAL hpvv=hpv1[1]*hpv2[0];
    REAL hpsv=hps1[1]*hpv2[0];
    REAL hpvs=hpv1[1]*hps2[0];
    REAL hpss=hps1[1]*hps2[0];
    REAL hpssv=hpss1[1]*hpv2[0];
    REAL hpvss=hpv1[1]*hpss2[0];
    if (ev[2]<ev[3]) {
      for (int i=2;i<n_pols_1d;i++) {
        Tensor<2,3> &ggv=gg[iv++];
        REAL hpvi=hpv0[i];
        REAL hpsi=hps0[i];
        ggv[0][0]=hpss0[i]*hpvv;
        ggv[1][0]=hpsi*hpsv;
        ggv[2][0]=hpsi*hpvs;
        ggv[1][1]=hpvi*hpssv;
        ggv[2][1]=hpvi*hpss;
        ggv[2][2]=hpvi*hpvss;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    } else {
      for (int i=2;i<n_pols_1d;i++) {
        Tensor<2,3> &ggv=gg[iv++];
        REAL hpvi=hpv0[i];
        REAL hpsi=hps0[i];
        ggv[0][0]=hpss0[i]*hpvv;
        ggv[1][0]=hpsi*hpsv;
        ggv[2][0]=hpsi*hpvs;
        ggv[1][1]=hpvi*hpssv;
        ggv[2][1]=hpvi*hpss;
        ggv[2][2]=hpvi*hpvss;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
        hpvv=-hpvv;
        hpsv=-hpsv;
        hpvs=-hpvs;
        hpss=-hpss;
        hpssv=-hpssv;
        hpvss=-hpvss;
      }
    }
  }
  {
    REAL hpvv=hpv0[0]*hpv2[1];
    REAL hpsv=hps0[0]*hpv2[1];
    REAL hpvs=hpv0[0]*hps2[1];
    REAL hpss=hps0[0]*hps2[1];
    REAL hpssv=hpss0[0]*hpv2[1];
    REAL hpvss=hpv0[0]*hpss2[1];
    if (ev[4]<ev[6]) {
      for (int j=2;j<n_pols_1d;j++) {
        Tensor<2,3> &ggv=gg[iv++];
        REAL hpvj=hpv1[j];
        REAL hpsj=hps1[j];
        ggv[0][0]=hpssv*hpvj;
        ggv[1][0]=hpsv*hpsj;
        ggv[2][0]=hpss*hpvj;
        ggv[1][1]=hpvv*hpss1[j];
        ggv[2][1]=hpvs*hpsj;
        ggv[2][2]=hpvv*hpvss;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    } else {
      for (int j=2;j<n_pols_1d;j++) {
        Tensor<2,3> &ggv=gg[iv++];
        REAL hpvj=hpv1[j];
        REAL hpsj=hps1[j];
        ggv[0][0]=hpssv*hpvj;
        ggv[1][0]=hpsv*hpsj;
        ggv[2][0]=hpss*hpvj;
        ggv[1][1]=hpvv*hpss1[j];
        ggv[2][1]=hpvs*hpsj;
        ggv[2][2]=hpvv*hpvss;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
        hpvv=-hpvv;
        hpsv=-hpsv;
        hpvs=-hpvs;
        hpss=-hpss;
        hpssv=-hpssv;
        hpvss=-hpvss;
      }
    }
  }
  {
    REAL hpvv=hpv0[1]*hpv2[1];
    REAL hpsv=hps0[1]*hpv2[1];
    REAL hpvs=hpv0[1]*hps2[1];
    REAL hpss=hps0[1]*hps2[1];
    REAL hpssv=hpss0[1]*hpv2[1];
    REAL hpvss=hpv0[1]*hpss2[1];
    if (ev[5]<ev[7]) {
      for (int j=2;j<n_pols_1d;j++) {
        Tensor<2,3> &ggv=gg[iv++];
        REAL hpvj=hpv1[j];
        REAL hpsj=hps1[j];
        ggv[0][0]=hpssv*hpvj;
        ggv[1][0]=hpsv*hpsj;
        ggv[2][0]=hpss*hpvj;
        ggv[1][1]=hpvv*hpss1[j];
        ggv[2][1]=hpvs*hpsj;
        ggv[2][2]=hpvv*hpvss;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    } else {
      for (int j=2;j<n_pols_1d;j++) {
        Tensor<2,3> &ggv=gg[iv++];
        REAL hpvj=hpv1[j];
        REAL hpsj=hps1[j];
        ggv[0][0]=hpssv*hpvj;
        ggv[1][0]=hpsv*hpsj;
        ggv[2][0]=hpss*hpvj;
        ggv[1][1]=hpvv*hpss1[j];
        ggv[2][1]=hpvs*hpsj;
        ggv[2][2]=hpvv*hpvss;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
        hpvv=-hpvv;
        hpsv=-hpsv;
        hpvs=-hpvs;
        hpss=-hpss;
        hpssv=-hpssv;
        hpvss=-hpvss;
      }
    }
  }
  {
    REAL hpvv=hpv1[0]*hpv2[1];
    REAL hpsv=hps1[0]*hpv2[1];
    REAL hpvs=hpv1[0]*hps2[1];
    REAL hpss=hps1[0]*hps2[1];
    REAL hpssv=hpss1[0]*hpv2[1];
    REAL hpvss=hpv1[0]*hpss2[1];
    if (ev[4]<ev[5]) {
      for (int i=2;i<n_pols_1d;i++) {
        Tensor<2,3> &ggv=gg[iv++];
        REAL hpvi=hpv0[i];
        REAL hpsi=hps0[i];
        ggv[0][0]=hpss0[i]*hpvv;
        ggv[1][0]=hpsi*hpsv;
        ggv[2][0]=hpsi*hpvs;
        ggv[1][1]=hpvi*hpssv;
        ggv[2][1]=hpvi*hpss;
        ggv[2][2]=hpvi*hpvss;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    } else {
      for (int i=2;i<n_pols_1d;i++) {
        Tensor<2,3> &ggv=gg[iv++];
        REAL hpvi=hpv0[i];
        REAL hpsi=hps0[i];
        ggv[0][0]=hpss0[i]*hpvv;
        ggv[1][0]=hpsi*hpsv;
        ggv[2][0]=hpsi*hpvs;
        ggv[1][1]=hpvi*hpssv;
        ggv[2][1]=hpvi*hpss;
        ggv[2][2]=hpvi*hpvss;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
        hpvv=-hpvv;
        hpsv=-hpsv;
        hpvs=-hpvs;
        hpss=-hpss;
        hpssv=-hpssv;
        hpvss=-hpvss;
      }
    }
  }
  {
    REAL hpvv=hpv1[1]*hpv2[1];
    REAL hpsv=hps1[1]*hpv2[1];
    REAL hpvs=hpv1[1]*hps2[1];
    REAL hpss=hps1[1]*hps2[1];
    REAL hpssv=hpss1[1]*hpv2[1];
    REAL hpvss=hpv1[1]*hpss2[1];
    if (ev[6]<ev[7]) {
      for (int i=2;i<n_pols_1d;i++) {
        Tensor<2,3> &ggv=gg[iv++];
        REAL hpvi=hpv0[i];
        REAL hpsi=hps0[i];
        ggv[0][0]=hpss0[i]*hpvv;
        ggv[1][0]=hpsi*hpsv;
        ggv[2][0]=hpsi*hpvs;
        ggv[1][1]=hpvi*hpssv;
        ggv[2][1]=hpvi*hpss;
        ggv[2][2]=hpvi*hpvss;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    } else {
      for (int i=2;i<n_pols_1d;i++) {
        Tensor<2,3> &ggv=gg[iv++];
        REAL hpvi=hpv0[i];
        REAL hpsi=hps0[i];
        ggv[0][0]=hpss0[i]*hpvv;
        ggv[1][0]=hpsi*hpsv;
        ggv[2][0]=hpsi*hpvs;
        ggv[1][1]=hpvi*hpssv;
        ggv[2][1]=hpvi*hpss;
        ggv[2][2]=hpvi*hpvss;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
        hpvv=-hpvv;
        hpsv=-hpsv;
        hpvs=-hpvs;
        hpss=-hpss;
        hpssv=-hpssv;
        hpvss=-hpvss;
      }
    }
  }
  {
    CHECK_TEST(ev[0]<ev[4]);
    REAL hpvv=hpv0[0]*hpv1[0];
    REAL hpsv=hps0[0]*hpv1[0];
    REAL hpvs=hpv0[0]*hps1[0];
    REAL hpss=hps0[0]*hps1[0];
    REAL hpssv=hpss0[0]*hpv1[0];
    REAL hpvss=hpv0[0]*hpss1[0];
    for (int k=2;k<n_pols_1d;k++) {
      Tensor<2,3> &ggv=gg[iv++];
      REAL hpvk=hpv2[k];
      REAL hpsk=hps2[k];
      ggv[0][0]=hpssv*hpvk;
      ggv[1][0]=hpss*hpvk;
      ggv[2][0]=hpsv*hpsk;
      ggv[1][1]=hpvss*hpvk;
      ggv[2][1]=hpvs*hpsk;
      ggv[2][2]=hpvv*hpss2[k];
      ggv[0][1]=ggv[1][0];
      ggv[0][2]=ggv[2][0];
      ggv[1][2]=ggv[2][1];
    }
  }
  {
    REAL hpvv=hpv0[1]*hpv1[0];
    REAL hpsv=hps0[1]*hpv1[0];
    REAL hpvs=hpv0[1]*hps1[0];
    REAL hpss=hps0[1]*hps1[0];
    REAL hpssv=hpss0[1]*hpv1[0];
    REAL hpvss=hpv0[1]*hpss1[0];
    if (ev[1]<ev[5]) {
      for (int k=2;k<n_pols_1d;k++) {
        Tensor<2,3> &ggv=gg[iv++];
        REAL hpvk=hpv2[k];
        REAL hpsk=hps2[k];
        ggv[0][0]=hpssv*hpvk;
        ggv[1][0]=hpss*hpvk;
        ggv[2][0]=hpsv*hpsk;
        ggv[1][1]=hpvss*hpvk;
        ggv[2][1]=hpvs*hpsk;
        ggv[2][2]=hpvv*hpss2[k];
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    } else {
      for (int k=2;k<n_pols_1d;k++) {
        Tensor<2,3> &ggv=gg[iv++];
        REAL hpvk=hpv2[k];
        REAL hpsk=hps2[k];
        ggv[0][0]=hpssv*hpvk;
        ggv[1][0]=hpss*hpvk;
        ggv[2][0]=hpsv*hpsk;
        ggv[1][1]=hpvss*hpvk;
        ggv[2][1]=hpvs*hpsk;
        ggv[2][2]=hpvv*hpss2[k];
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
        hpvv=-hpvv;
        hpsv=-hpsv;
        hpvs=-hpvs;
        hpss=-hpss;
        hpssv=-hpssv;
        hpvss=-hpvss;
      }
    }
  }
  {
    REAL hpvv=hpv0[0]*hpv1[1];
    REAL hpsv=hps0[0]*hpv1[1];
    REAL hpvs=hpv0[0]*hps1[1];
    REAL hpss=hps0[0]*hps1[1];
    REAL hpssv=hpss0[0]*hpv1[1];
    REAL hpvss=hpv0[0]*hpss1[1];
    if (ev[2]<ev[6]) {
      for (int k=2;k<n_pols_1d;k++) {
        Tensor<2,3> &ggv=gg[iv++];
        REAL hpvk=hpv2[k];
        REAL hpsk=hps2[k];
        ggv[0][0]=hpssv*hpvk;
        ggv[1][0]=hpss*hpvk;
        ggv[2][0]=hpsv*hpsk;
        ggv[1][1]=hpvss*hpvk;
        ggv[2][1]=hpvs*hpsk;
        ggv[2][2]=hpvv*hpss2[k];
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    } else {
      for (int k=2;k<n_pols_1d;k++) {
        Tensor<2,3> &ggv=gg[iv++];
        REAL hpvk=hpv2[k];
        REAL hpsk=hps2[k];
        ggv[0][0]=hpssv*hpvk;
        ggv[1][0]=hpss*hpvk;
        ggv[2][0]=hpsv*hpsk;
        ggv[1][1]=hpvss*hpvk;
        ggv[2][1]=hpvs*hpsk;
        ggv[2][2]=hpvv*hpss2[k];
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
        hpvv=-hpvv;
        hpsv=-hpsv;
        hpvs=-hpvs;
        hpss=-hpss;
        hpssv=-hpssv;
        hpvss=-hpvss;
      }
    }
  }
  {
    REAL hpvv=hpv0[1]*hpv1[1];
    REAL hpsv=hps0[1]*hpv1[1];
    REAL hpvs=hpv0[1]*hps1[1];
    REAL hpss=hps0[1]*hps1[1];
    REAL hpssv=hpss0[1]*hpv1[1];
    REAL hpvss=hpv0[1]*hpss1[1];
    if (ev[3]<ev[7]) {
      for (int k=2;k<n_pols_1d;k++) {
        Tensor<2,3> &ggv=gg[iv++];
        REAL hpvk=hpv2[k];
        REAL hpsk=hps2[k];
        ggv[0][0]=hpssv*hpvk;
        ggv[1][0]=hpss*hpvk;
        ggv[2][0]=hpsv*hpsk;
        ggv[1][1]=hpvss*hpvk;
        ggv[2][1]=hpvs*hpsk;
        ggv[2][2]=hpvv*hpss2[k];
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    } else {
      for (int k=2;k<n_pols_1d;k++) {
        Tensor<2,3> &ggv=gg[iv++];
        REAL hpvk=hpv2[k];
        REAL hpsk=hps2[k];
        ggv[0][0]=hpssv*hpvk;
        ggv[1][0]=hpss*hpvk;
        ggv[2][0]=hpsv*hpsk;
        ggv[1][1]=hpvss*hpvk;
        ggv[2][1]=hpvs*hpsk;
        ggv[2][2]=hpvv*hpss2[k];
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
        hpvv=-hpvv;
        hpsv=-hpsv;
        hpvs=-hpvs;
        hpss=-hpss;
        hpssv=-hpssv;
        hpvss=-hpvss;
      }
    }
  }
  {
    CHECK_TEST(ev[0]<min(ev[2],min(ev[4],ev[6])));
    if (ev[2]<ev[4]) {
      for (int k=2;k<n_pols_1d;k++) {
        REAL hpvv=hpv0[0]*hpv2[k];
        REAL hpsv=hps0[0]*hpv2[k];
        REAL hpvs=hpv0[0]*hps2[k];
        REAL hpss=hps0[0]*hps2[k];
        REAL hpssv=hpss0[0]*hpv2[k];
        REAL hpvss=hpv0[0]*hpss2[k];
        for (int j=2;j<n_pols_1d;j++) {
          Tensor<2,3> &ggv=gg[iv++];
          REAL hpvj=hpv1[j];
          REAL hpsj=hps1[j];
          ggv[0][0]=hpssv*hpvj;
          ggv[1][0]=hpsv*hpsj;
          ggv[2][0]=hpss*hpvj;
          ggv[1][1]=hpvv*hpss1[j];
          ggv[2][1]=hpvs*hpsj;
          ggv[2][2]=hpvss*hpvj;
          ggv[0][1]=ggv[1][0];
          ggv[0][2]=ggv[2][0];
          ggv[1][2]=ggv[2][1];
        }
      }
    } else {
      for (int j=2;j<n_pols_1d;j++) {
        REAL hpvv=hpv0[0]*hpv1[j];
        REAL hpsv=hps0[0]*hpv1[j];
        REAL hpvs=hpv0[0]*hps1[j];
        REAL hpss=hps0[0]*hps1[j];
        REAL hpssv=hpss0[0]*hpv1[j];
        REAL hpvss=hpv0[0]*hpss1[j];
        for (int k=2;k<n_pols_1d;k++) {
          Tensor<2,3> &ggv=gg[iv++];
          REAL hpvk=hpv2[k];
          REAL hpsk=hps2[k];
          ggv[0][0]=hpssv*hpvk;
          ggv[1][0]=hpss*hpvk;
          ggv[2][0]=hpsv*hpsk;
          ggv[1][1]=hpvss*hpvk;
          ggv[2][1]=hpvs*hpsk;
          ggv[2][2]=hpvv*hpss2[k];
          ggv[0][1]=ggv[1][0];
          ggv[0][2]=ggv[2][0];
          ggv[1][2]=ggv[2][1];
        }
      }
    }
  }
  {
    int vmin=min(min(ev[1],ev[3]),min(ev[5],ev[7]));
    if (ev[1]==vmin) {
      if (ev[3]<ev[5]) {
        for (int k=2;k<n_pols_1d;k++) {
          REAL hpvv=hpv0[1]*hpv2[k];
          REAL hpsv=hps0[1]*hpv2[k];
          REAL hpvs=hpv0[1]*hps2[k];
          REAL hpss=hps0[1]*hps2[k];
          REAL hpssv=hpss0[1]*hpv2[k];
          REAL hpvss=hpv0[1]*hpss2[k];
          for (int j=2;j<n_pols_1d;j++) {
            Tensor<2,3> &ggv=gg[iv++];
            REAL hpvj=hpv1[j];
            REAL hpsj=hps1[j];
            ggv[0][0]=hpssv*hpvj;
            ggv[1][0]=hpsv*hpsj;
            ggv[2][0]=hpss*hpvj;
            ggv[1][1]=hpvv*hpss1[j];
            ggv[2][1]=hpvs*hpsj;
            ggv[2][2]=hpvss*hpvj;
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
          }
        }
      } else {
        for (int j=2;j<n_pols_1d;j++) {
          REAL hpvv=hpv0[1]*hpv1[j];
          REAL hpsv=hps0[1]*hpv1[j];
          REAL hpvs=hpv0[1]*hps1[j];
          REAL hpss=hps0[1]*hps1[j];
          REAL hpssv=hpss0[1]*hpv1[j];
          REAL hpvss=hpv0[1]*hpss1[j];
          for (int k=2;k<n_pols_1d;k++) {
            Tensor<2,3> &ggv=gg[iv++];
            REAL hpvk=hpv2[k];
            REAL hpsk=hps2[k];
            ggv[0][0]=hpssv*hpvk;
            ggv[1][0]=hpss*hpvk;
            ggv[2][0]=hpsv*hpsk;
            ggv[1][1]=hpvss*hpvk;
            ggv[2][1]=hpvs*hpsk;
            ggv[2][2]=hpvv*hpss2[k];
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
          }
        }
      }
    } else if (ev[3]==vmin) {
      if (ev[1]<ev[7]) {
        for (int k=2;k<n_pols_1d;k++) {
          REAL hpvv=hpv0[1]*hpv2[k];
          REAL hpsv=hps0[1]*hpv2[k];
          REAL hpvs=hpv0[1]*hps2[k];
          REAL hpss=hps0[1]*hps2[k];
          REAL hpssv=hpss0[1]*hpv2[k];
          REAL hpvss=hpv0[1]*hpss2[k];
          for (int j=2;j<n_pols_1d;j++) {
            Tensor<2,3> &ggv=gg[iv++];
            REAL hpvj=hpv1[j];
            REAL hpsj=hps1[j];
            ggv[0][0]=hpssv*hpvj;
            ggv[1][0]=hpsv*hpsj;
            ggv[2][0]=hpss*hpvj;
            ggv[1][1]=hpvv*hpss1[j];
            ggv[2][1]=hpvs*hpsj;
            ggv[2][2]=hpvss*hpvj;
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
            hpvv=-hpvv;
            hpsv=-hpsv;
            hpvs=-hpvs;
            hpss=-hpss;
            hpssv=-hpssv;
            hpvss=-hpvss;
          }
        }
      } else {
        REAL hpvi=hpv0[1];
        REAL hpsi=hps0[1];
        REAL hpssi=hpss0[1];
        for (int j=2;j<n_pols_1d;j++) {
          REAL hpvv=hpvi*hpv1[j];
          REAL hpsv=hpsi*hpv1[j];
          REAL hpvs=hpvi*hps1[j];
          REAL hpss=hpsi*hps1[j];
          REAL hpssv=hpssi*hpv1[j];
          REAL hpvss=hpvi*hpss1[j];
          for (int k=2;k<n_pols_1d;k++) {
            Tensor<2,3> &ggv=gg[iv++];
            REAL hpvk=hpv2[k];
            REAL hpsk=hps2[k];
            ggv[0][0]=hpssv*hpvk;
            ggv[1][0]=hpss*hpvk;
            ggv[2][0]=hpsv*hpsk;
            ggv[1][1]=hpvss*hpvk;
            ggv[2][1]=hpvs*hpsk;
            ggv[2][2]=hpvv*hpss2[k];
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
          }
          hpvi=-hpvi;
          hpsi=-hpsi;
          hpssi=-hpssi;
        }
      }
    } else if (ev[5]==vmin) {
      if (ev[1]<ev[7]) {
        REAL hpvi=hpv0[1];
        REAL hpsi=hps0[1];
        REAL hpssi=hpss0[1];
        for (int j=2;j<n_pols_1d;j++) {
          REAL hpvv=hpvi*hpv1[j];
          REAL hpsv=hpsi*hpv1[j];
          REAL hpvs=hpvi*hps1[j];
          REAL hpss=hpsi*hps1[j];
          REAL hpssv=hpssi*hpv1[j];
          REAL hpvss=hpvi*hpss1[j];
          for (int k=2;k<n_pols_1d;k++) {
            Tensor<2,3> &ggv=gg[iv++];
            REAL hpvk=hpv2[k];
            REAL hpsk=hps2[k];
            ggv[0][0]=hpssv*hpvk;
            ggv[1][0]=hpss*hpvk;
            ggv[2][0]=hpsv*hpsk;
            ggv[1][1]=hpvss*hpvk;
            ggv[2][1]=hpvs*hpsk;
            ggv[2][2]=hpvv*hpss2[k];
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
            hpvv=-hpvv;
            hpsv=-hpsv;
            hpvs=-hpvs;
            hpss=-hpss;
            hpssv=-hpssv;
            hpvss=-hpvss;
          }
        }
      } else {
        REAL hpvi=hpv0[1];
        REAL hpsi=hps0[1];
        REAL hpssi=hpss0[1];
        for (int k=2;k<n_pols_1d;k++) {
          REAL hpvv=hpvi*hpv2[k];
          REAL hpsv=hpsi*hpv2[k];
          REAL hpvs=hpvi*hps2[k];
          REAL hpss=hpsi*hps2[k];
          REAL hpssv=hpssi*hpv2[k];
          REAL hpvss=hpvi*hpss2[k];
          for (int j=2;j<n_pols_1d;j++) {
            Tensor<2,3> &ggv=gg[iv++];
            REAL hpvj=hpv1[j];
            REAL hpsj=hps1[j];
            ggv[0][0]=hpssv*hpvj;
            ggv[1][0]=hpsv*hpsj;
            ggv[2][0]=hpss*hpvj;
            ggv[1][1]=hpvv*hpss1[j];
            ggv[2][1]=hpvs*hpsj;
            ggv[2][2]=hpvss*hpvj;
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
          }
          hpvi=-hpvi;
          hpsi=-hpsi;
          hpssi=-hpssi;
        }
      }
    } else {
      if (ev[3]<ev[5]) {
        REAL hpvi=hpv0[1];
        REAL hpsi=hps0[1];
        REAL hpssi=hpss0[1];
        for (int j=2;j<n_pols_1d;j++) {
          REAL hpvv=hpvi*hpv1[j];
          REAL hpsv=hpsi*hpv1[j];
          REAL hpvs=hpvi*hps1[j];
          REAL hpss=hpsi*hps1[j];
          REAL hpssv=hpssi*hpv1[j];
          REAL hpvss=hpvi*hpss1[j];
          for (int k=2;k<n_pols_1d;k++) {
            Tensor<2,3> &ggv=gg[iv++];
            REAL hpvk=hpv2[k];
            REAL hpsk=hps2[k];
            ggv[0][0]=hpssv*hpvk;
            ggv[1][0]=hpss*hpvk;
            ggv[2][0]=hpsv*hpsk;
            ggv[1][1]=hpvss*hpvk;
            ggv[2][1]=hpvs*hpsk;
            ggv[2][2]=hpvv*hpss2[k];
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
            hpvv=-hpvv;
            hpsv=-hpsv;
            hpvs=-hpvs;
            hpss=-hpss;
            hpssv=-hpssv;
            hpvss=-hpvss;
          }
          hpvi=-hpvi;
          hpsi=-hpsi;
          hpssi=-hpssi;
        }
      } else {
        REAL hpvi=hpv0[1];
        REAL hpsi=hps0[1];
        REAL hpssi=hpss0[1];
        for (int k=2;k<n_pols_1d;k++) {
          REAL hpvv=hpvi*hpv2[k];
          REAL hpsv=hpsi*hpv2[k];
          REAL hpvs=hpvi*hps2[k];
          REAL hpss=hpsi*hps2[k];
          REAL hpssv=hpssi*hpv2[k];
          REAL hpvss=hpvi*hpss2[k];
          for (int j=2;j<n_pols_1d;j++) {
            Tensor<2,3> &ggv=gg[iv++];
            REAL hpvj=hpv1[j];
            REAL hpsj=hps1[j];
            ggv[0][0]=hpssv*hpvj;
            ggv[1][0]=hpsv*hpsj;
            ggv[2][0]=hpss*hpvj;
            ggv[1][1]=hpvv*hpss1[j];
            ggv[2][1]=hpvs*hpsj;
            ggv[2][2]=hpvss*hpvj;
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
            hpvv=-hpvv;
            hpsv=-hpsv;
            hpvs=-hpvs;
            hpss=-hpss;
            hpssv=-hpssv;
            hpvss=-hpvss;
          }
          hpvi=-hpvi;
          hpsi=-hpsi;
          hpssi=-hpssi;
        }
      }
    }
  }
  {
    CHECK_TEST(ev[0]<min(ev[1],min(ev[4],ev[5])) && ev[1]<ev[4]);
    const REAL &hpvj=hpv1[0];
    for (int k=2;k<n_pols_1d;k++) {
      REAL hpvv=hpv1[0]*hpv2[k];
      REAL hpsv=hps1[0]*hpv2[k];
      REAL hpvs=hpv1[0]*hps2[k];
      REAL hpss=hps1[0]*hps2[k];
      REAL hpssv=hpss1[0]*hpv2[k];
      REAL hpvss=hpv1[0]*hpss2[k];
      for (int i=2;i<n_pols_1d;i++) {
        Tensor<2,3> &ggv=gg[iv++];
        REAL hpvi=hpv0[i];
        REAL hpsi=hps0[i];
        ggv[0][0]=hpss0[i]*hpvv;
        ggv[1][0]=hpsi*hpsv;
        ggv[2][0]=hpsi*hpvs;
        ggv[1][1]=hpvi*hpssv;
        ggv[2][1]=hpvi*hpss;
        ggv[2][2]=hpvi*hpvss;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    }
  }
  {
    REAL hpvj=hpv1[1];
    int vmin=min(min(ev[2],ev[3]),min(ev[6],ev[7]));
    if (ev[2]==vmin) {
      if (ev[3]<ev[6]) {
        for (int k=2;k<n_pols_1d;k++) {
          REAL hpvv=hpv1[1]*hpv2[k];
          REAL hpsv=hps1[1]*hpv2[k];
          REAL hpvs=hpv1[1]*hps2[k];
          REAL hpss=hps1[1]*hps2[k];
          REAL hpssv=hpss1[1]*hpv2[k];
          REAL hpvss=hpv1[1]*hpss2[k];
          for (int i=2;i<n_pols_1d;i++) {
            Tensor<2,3> &ggv=gg[iv++];
            REAL hpvi=hpv0[i];
            REAL hpsi=hps0[i];
            ggv[0][0]=hpss0[i]*hpvv;
            ggv[1][0]=hpsi*hpsv;
            ggv[2][0]=hpsi*hpvs;
            ggv[1][1]=hpvi*hpssv;
            ggv[2][1]=hpvi*hpss;
            ggv[2][2]=hpvi*hpvss;
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
          }
        }
      } else {
        for (int i=2;i<n_pols_1d;i++) {
          REAL hpvv=hpv0[i]*hpv1[1];
          REAL hpsv=hps0[i]*hpv1[1];
          REAL hpvs=hpv0[i]*hps1[1];
          REAL hpss=hps0[i]*hps1[1];
          REAL hpssv=hpss0[i]*hpv1[1];
          REAL hpvss=hpv0[i]*hpss1[1];
          for (int k=2;k<n_pols_1d;k++) {
            Tensor<2,3> &ggv=gg[iv++];
            REAL hpvk=hpv2[k];
            REAL hpsk=hps2[k];
            ggv[0][0]=hpssv*hpvk;
            ggv[1][0]=hpss*hpvk;
            ggv[2][0]=hpsv*hpsk;
            ggv[1][1]=hpvss*hpvk;
            ggv[2][1]=hpvs*hpsk;
            ggv[2][2]=hpvv*hpss2[k];
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
          }
        }
      }
    } else if (ev[3]==vmin) {
      if (ev[2]<ev[7]) {
        for (int k=2;k<n_pols_1d;k++) {
          REAL hpvv=hpv1[1]*hpv2[k];
          REAL hpsv=hps1[1]*hpv2[k];
          REAL hpvs=hpv1[1]*hps2[k];
          REAL hpss=hps1[1]*hps2[k];
          REAL hpssv=hpss1[1]*hpv2[k];
          REAL hpvss=hpv1[1]*hpss2[k];
          for (int i=2;i<n_pols_1d;i++) {
            Tensor<2,3> &ggv=gg[iv++];
            REAL hpvi=hpv0[i];
            REAL hpsi=hps0[i];
            ggv[0][0]=hpss0[i]*hpvv;
            ggv[1][0]=hpsi*hpsv;
            ggv[2][0]=hpsi*hpvs;
            ggv[1][1]=hpvi*hpssv;
            ggv[2][1]=hpvi*hpss;
            ggv[2][2]=hpvi*hpvss;
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
            hpvv=-hpvv;
            hpsv=-hpsv;
            hpvs=-hpvs;
            hpss=-hpss;
            hpssv=-hpssv;
            hpvss=-hpvss;
          }
        }
      } else {
        REAL hpvj=hpv1[1];
        REAL hpsj=hps1[1];
        REAL hpssj=hpss1[1];
        for (int i=2;i<n_pols_1d;i++) {
          REAL hpvv=hpv0[i]*hpvj;
          REAL hpsv=hps0[i]*hpvj;
          REAL hpvs=hpv0[i]*hpsj;
          REAL hpss=hps0[i]*hpsj;
          REAL hpssv=hpss0[i]*hpvj;
          REAL hpvss=hpv0[i]*hpssj;
          for (int k=2;k<n_pols_1d;k++) {
            Tensor<2,3> &ggv=gg[iv++];
            REAL hpvk=hpv2[k];
            REAL hpsk=hps2[k];
            ggv[0][0]=hpssv*hpvk;
            ggv[1][0]=hpss*hpvk;
            ggv[2][0]=hpsv*hpsk;
            ggv[1][1]=hpvss*hpvk;
            ggv[2][1]=hpvs*hpsk;
            ggv[2][2]=hpvv*hpss2[k];
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
          }
          hpvj=-hpvj;
          hpsj=-hpsj;
          hpssj=-hpssj;
        }
      }
    } else if (ev[6]==vmin) {
      if (ev[2]<ev[7]) {
        for (int i=2;i<n_pols_1d;i++) {
          REAL hpvv=hpv0[i]*hpv1[1];
          REAL hpsv=hps0[i]*hpv1[1];
          REAL hpvs=hpv0[i]*hps1[1];
          REAL hpss=hps0[i]*hps1[1];
          REAL hpssv=hpss0[i]*hpv1[1];
          REAL hpvss=hpv0[i]*hpss1[1];
          for (int k=2;k<n_pols_1d;k++) {
            Tensor<2,3> &ggv=gg[iv++];
            REAL hpvk=hpv2[k];
            REAL hpsk=hps2[k];
            ggv[0][0]=hpssv*hpvk;
            ggv[1][0]=hpss*hpvk;
            ggv[2][0]=hpsv*hpsk;
            ggv[1][1]=hpvss*hpvk;
            ggv[2][1]=hpvs*hpsk;
            ggv[2][2]=hpvv*hpss2[k];
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
            hpvv=-hpvv;
            hpsv=-hpsv;
            hpvs=-hpvs;
            hpss=-hpss;
            hpssv=-hpssv;
            hpvss=-hpvss;
          }
        }
      } else {
        REAL hpvj=hpv1[1];
        REAL hpsj=hps1[1];
        REAL hpssj=hpss1[1];
        for (int k=2;k<n_pols_1d;k++) {
          REAL hpvv=hpvj*hpv2[k];
          REAL hpsv=hpsj*hpv2[k];
          REAL hpvs=hpvj*hps2[k];
          REAL hpss=hpsj*hps2[k];
          REAL hpssv=hpssj*hpv2[k];
          REAL hpvss=hpvj*hpss2[k];
          for (int i=2;i<n_pols_1d;i++) {
            Tensor<2,3> &ggv=gg[iv++];
            REAL hpvi=hpv0[i];
            REAL hpsi=hps0[i];
            ggv[0][0]=hpss0[i]*hpvv;
            ggv[1][0]=hpsi*hpsv;
            ggv[2][0]=hpsi*hpvs;
            ggv[1][1]=hpvi*hpssv;
            ggv[2][1]=hpvi*hpss;
            ggv[2][2]=hpvi*hpvss;
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
          }
          hpvj=-hpvj;
          hpsj=-hpsj;
          hpssj=-hpssj;
        }
      }
    } else {
      if (ev[3]<ev[6]) {
        REAL hpvj=hpv1[1];
        REAL hpsj=hps1[1];
        REAL hpssj=hpss1[1];
        for (int i=2;i<n_pols_1d;i++) {
          REAL hpvv=hpv0[i]*hpvj;
          REAL hpsv=hps0[i]*hpvj;
          REAL hpvs=hpv0[i]*hpsj;
          REAL hpss=hps0[i]*hpsj;
          REAL hpssv=hpss0[i]*hpvj;
          REAL hpvss=hpv0[i]*hpssj;
          for (int k=2;k<n_pols_1d;k++) {
            Tensor<2,3> &ggv=gg[iv++];
            REAL hpvk=hpv2[k];
            REAL hpsk=hps2[k];
            ggv[0][0]=hpssv*hpvk;
            ggv[1][0]=hpss*hpvk;
            ggv[2][0]=hpsv*hpsk;
            ggv[1][1]=hpvss*hpvk;
            ggv[2][1]=hpvs*hpsk;
            ggv[2][2]=hpvv*hpss2[k];
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
            hpvv=-hpvv;
            hpsv=-hpsv;
            hpvs=-hpvs;
            hpss=-hpss;
            hpssv=-hpssv;
            hpvss=-hpvss;
          }
          hpvj=-hpvj;
          hpsj=-hpsj;
          hpssj=-hpssj;
        }
      } else {
        REAL hpvj=hpv1[1];
        REAL hpsj=hps1[1];
        REAL hpssj=hpss1[1];
        for (int k=2;k<n_pols_1d;k++) {
          REAL hpvv=hpvj*hpv2[k];
          REAL hpsv=hpsj*hpv2[k];
          REAL hpvs=hpvj*hps2[k];
          REAL hpss=hpsj*hps2[k];
          REAL hpssv=hpssj*hpv2[k];
          REAL hpvss=hpvj*hpss2[k];
          for (int i=2;i<n_pols_1d;i++) {
            Tensor<2,3> &ggv=gg[iv++];
            REAL hpvi=hpv0[i];
            REAL hpsi=hps0[i];
            ggv[0][0]=hpss0[i]*hpvv;
            ggv[1][0]=hpsi*hpsv;
            ggv[2][0]=hpsi*hpvs;
            ggv[1][1]=hpvi*hpssv;
            ggv[2][1]=hpvi*hpss;
            ggv[2][2]=hpvi*hpvss;
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
            hpvv=-hpvv;
            hpsv=-hpsv;
            hpvs=-hpvs;
            hpss=-hpss;
            hpssv=-hpssv;
            hpvss=-hpvss;
          }
          hpvj=-hpvj;
          hpsj=-hpsj;
          hpssj=-hpssj;
        }
      }
    }
  }
  {
    CHECK_TEST(ev[0]<min(ev[1],min(ev[2],ev[3])) && ev[1]<ev[2]);
    for (int j=2;j<n_pols_1d;j++) {
      REAL hpvv=hpv1[j]*hpv2[0];
      REAL hpsv=hps1[j]*hpv2[0];
      REAL hpvs=hpv1[j]*hps2[0];
      REAL hpss=hps1[j]*hps2[0];
      REAL hpssv=hpss1[j]*hpv2[0];
      REAL hpvss=hpv1[j]*hpss2[0];
      for (int i=2;i<n_pols_1d;i++) {
        Tensor<2,3> &ggv=gg[iv++];
        REAL hpvi=hpv0[i];
        REAL hpsi=hps0[i];
        ggv[0][0]=hpss0[i]*hpvv;
        ggv[1][0]=hpsi*hpsv;
        ggv[2][0]=hpsi*hpvs;
        ggv[1][1]=hpvi*hpssv;
        ggv[2][1]=hpvi*hpss;
        ggv[2][2]=hpvi*hpvss;
        ggv[0][1]=ggv[1][0];
        ggv[0][2]=ggv[2][0];
        ggv[1][2]=ggv[2][1];
      }
    }
  }
  {
    REAL hpvk=hpv2[1];
    int vmin=min(min(ev[4],ev[5]),min(ev[6],ev[7]));
    if (ev[4]==vmin) {
      if (ev[5]<ev[6]) {
        for (int j=2;j<n_pols_1d;j++) {
          REAL hpvv=hpv1[j]*hpv2[1];
          REAL hpsv=hps1[j]*hpv2[1];
          REAL hpvs=hpv1[j]*hps2[1];
          REAL hpss=hps1[j]*hps2[1];
          REAL hpssv=hpss1[j]*hpv2[1];
          REAL hpvss=hpv1[j]*hpss2[1];
          for (int i=2;i<n_pols_1d;i++) {
            Tensor<2,3> &ggv=gg[iv++];
            REAL hpvi=hpv0[i];
            REAL hpsi=hps0[i];
            ggv[0][0]=hpss0[i]*hpvv;
            ggv[1][0]=hpsi*hpsv;
            ggv[2][0]=hpsi*hpvs;
            ggv[1][1]=hpvi*hpssv;
            ggv[2][1]=hpvi*hpss;
            ggv[2][2]=hpvi*hpvss;
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
          }
        }
      } else {
        for (int i=2;i<n_pols_1d;i++) {
          REAL hpvv=hpv0[i]*hpv2[1];
          REAL hpsv=hps0[i]*hpv2[1];
          REAL hpvs=hpv0[i]*hps2[1];
          REAL hpss=hps0[i]*hps2[1];
          REAL hpssv=hpss0[i]*hpv2[1];
          REAL hpvss=hpv0[i]*hpss2[1];
          for (int j=2;j<n_pols_1d;j++) {
            Tensor<2,3> &ggv=gg[iv++];
            REAL hpvj=hpv1[j];
            REAL hpsj=hps1[j];
            ggv[0][0]=hpssv*hpvj;
            ggv[1][0]=hpsv*hpsj;
            ggv[2][0]=hpss*hpvj;
            ggv[1][1]=hpvv*hpss1[j];
            ggv[2][1]=hpvs*hpsj;
            ggv[2][2]=hpvss*hpvj;
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
          }
        }
      }
    } else if (ev[5]==vmin) {
      if (ev[4]<ev[7]) {
        for (int j=2;j<n_pols_1d;j++) {
          REAL hpvv=hpv1[j]*hpv2[1];
          REAL hpsv=hps1[j]*hpv2[1];
          REAL hpvs=hpv1[j]*hps2[1];
          REAL hpss=hps1[j]*hps2[1];
          REAL hpssv=hpss1[j]*hpv2[1];
          REAL hpvss=hpv1[j]*hpss2[1];
          for (int i=2;i<n_pols_1d;i++) {
            Tensor<2,3> &ggv=gg[iv++];
            REAL hpvi=hpv0[i];
            REAL hpsi=hps0[i];
            ggv[0][0]=hpss0[i]*hpvv;
            ggv[1][0]=hpsi*hpsv;
            ggv[2][0]=hpsi*hpvs;
            ggv[1][1]=hpvi*hpssv;
            ggv[2][1]=hpvi*hpss;
            ggv[2][2]=hpvi*hpvss;
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
            hpvv=-hpvv;
            hpsv=-hpsv;
            hpvs=-hpvs;
            hpss=-hpss;
            hpssv=-hpssv;
            hpvss=-hpvss;
          }
        }
      } else {
        REAL hpvk=hpv2[1];
        REAL hpsk=hps2[1];
        REAL hpssk=hpss2[1];
        for (int i=2;i<n_pols_1d;i++) {
          REAL hpvv=hpv0[i]*hpvk;
          REAL hpsv=hps0[i]*hpvk;
          REAL hpvs=hpv0[i]*hpsk;
          REAL hpss=hps0[i]*hpsk;
          REAL hpssv=hpss0[i]*hpvk;
          REAL hpvss=hpv0[i]*hpssk;
          for (int j=2;j<n_pols_1d;j++) {
            Tensor<2,3> &ggv=gg[iv++];
            REAL hpvj=hpv1[j];
            REAL hpsj=hps1[j];
            ggv[0][0]=hpssv*hpvj;
            ggv[1][0]=hpsv*hpsj;
            ggv[2][0]=hpss*hpvj;
            ggv[1][1]=hpvv*hpss1[j];
            ggv[2][1]=hpvs*hpsj;
            ggv[2][2]=hpvss*hpvj;
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
          }
          hpvk=-hpvk;
          hpsk=-hpsk;
          hpssk=-hpssk;
        }
      }
    } else if (ev[6]==vmin) {
      if (ev[4]<ev[7]) {
        for (int i=2;i<n_pols_1d;i++) {
          REAL hpvv=hpv0[i]*hpv2[1];
          REAL hpsv=hps0[i]*hpv2[1];
          REAL hpvs=hpv0[i]*hps2[1];
          REAL hpss=hps0[i]*hps2[1];
          REAL hpssv=hpss0[i]*hpv2[1];
          REAL hpvss=hpv0[i]*hpss2[1];
          for (int j=2;j<n_pols_1d;j++) {
            Tensor<2,3> &ggv=gg[iv++];
            REAL hpvj=hpv1[j];
            REAL hpsj=hps1[j];
            ggv[0][0]=hpssv*hpvj;
            ggv[1][0]=hpsv*hpsj;
            ggv[2][0]=hpss*hpvj;
            ggv[1][1]=hpvv*hpss1[j];
            ggv[2][1]=hpvs*hpsj;
            ggv[2][2]=hpvss*hpvj;
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
            hpvv=-hpvv;
            hpsv=-hpsv;
            hpvs=-hpvs;
            hpss=-hpss;
            hpssv=-hpssv;
            hpvss=-hpvss;
          }
        }
      } else {
        REAL hpvk=hpv2[1];
        REAL hpsk=hps2[1];
        REAL hpssk=hpss2[1];
        for (int j=2;j<n_pols_1d;j++) {
          REAL hpvv=hpv1[j]*hpvk;
          REAL hpsv=hps1[j]*hpvk;
          REAL hpvs=hpv1[j]*hpsk;
          REAL hpss=hps1[j]*hpsk;
          REAL hpssv=hpss1[j]*hpvk;
          REAL hpvss=hpv1[j]*hpssk;
          for (int i=2;i<n_pols_1d;i++) {
            Tensor<2,3> &ggv=gg[iv++];
            REAL hpvi=hpv0[i];
            REAL hpsi=hps0[i];
            ggv[0][0]=hpss0[i]*hpvv;
            ggv[1][0]=hpsi*hpsv;
            ggv[2][0]=hpsi*hpvs;
            ggv[1][1]=hpvi*hpssv;
            ggv[2][1]=hpvi*hpss;
            ggv[2][2]=hpvi*hpvss;
            ggv[0][1]=ggv[1][0];
            ggv[0][2]=ggv[2][0];
            ggv[1][2]=ggv[2][1];
          }
          hpvk=-hpvk;
          hpsk=-hpsk;
          hpssk=-hpssk;
        }
      }
    } else {
      if (ev[5]<ev[6]) {
        REAL hpvk=hpv2[1];
        REAL hpsk=hps2[1];
        REAL hpssk=hpss2[1];
        for (int j=2;j<n_pols_1d;j++) {
          REAL hpvv=hpv0[j]*hpvk;
          REAL hpsv=hps0[j]*hpvk;
          REAL hpvs=hpv0[j]*hpsk;
          REAL hpss=hps0[j]*hpsk;
          REAL hpssv=hpss0[j]*hpvk;
          REAL hpvss=hpv0[j]*hpssk;
          for (int i=2;i<n_pols_1d;i++) {
            Tensor<2,3> &ggv=gg[iv++];
            REAL hpvi=hpv1[i];
            REAL hpsi=hps1[i];
            ggv[1][1]=hpss1[i]*hpvv;
            ggv[0][1]=hpsi*hpsv;
            ggv[2][1]=hpsi*hpvs;
            ggv[0][0]=hpvi*hpssv;
            ggv[2][0]=hpvi*hpss;
            ggv[2][2]=hpvi*hpvss;
            ggv[1][0]=ggv[0][1];
            ggv[1][2]=ggv[2][1];
            ggv[0][2]=ggv[2][0];
            hpvv=-hpvv;
            hpsv=-hpsv;
            hpvs=-hpvs;
            hpss=-hpss;
            hpssv=-hpssv;
            hpvss=-hpvss;
          }
          hpvk=-hpvk;
          hpsk=-hpsk;
          hpssk=-hpssk;
        }
      } else {
        REAL hpvk=hpv2[1];
        REAL hpsk=hps2[1];
        REAL hpssk=hpss2[1];
        for (int i=2;i<n_pols_1d;i++) {
          REAL hpvv=hpv1[i]*hpvk;
          REAL hpsv=hps1[i]*hpvk;
          REAL hpvs=hpv1[i]*hpsk;
          REAL hpss=hps1[i]*hpsk;
          REAL hpssv=hpss1[i]*hpvk;
          REAL hpvss=hpv1[i]*hpssk;
          for (int j=2;j<n_pols_1d;j++) {
            Tensor<2,3> &ggv=gg[iv++];
            REAL hpvj=hpv0[j];
            REAL hpsj=hps0[j];
            ggv[1][1]=hpssv*hpvj;
            ggv[0][1]=hpsv*hpsj;
            ggv[2][1]=hpss*hpvj;
            ggv[0][0]=hpvv*hpss0[j];
            ggv[2][0]=hpvs*hpsj;
            ggv[2][2]=hpvss*hpvj;
            ggv[1][0]=ggv[0][1];
            ggv[1][2]=ggv[2][1];
            ggv[0][2]=ggv[2][0];
            hpvv=-hpvv;
            hpsv=-hpsv;
            hpvs=-hpvs;
            hpss=-hpss;
            hpssv=-hpssv;
            hpvss=-hpvss;
          }
          hpvk=-hpvk;
          hpsk=-hpsk;
          hpssk=-hpssk;
        }
      }
    }
  }
  {
    for (int k=2;k<n_pols_1d;k++) {
      const REAL &hpvk=hpv2[k];
      const REAL &hpsk=hps2[k];
      const REAL &hpssk=hpss2[k];
      for (int j=2;j<n_pols_1d;j++) {
        REAL hpvv=hpv1[j]*hpvk;
        REAL hpsv=hps1[j]*hpvk;
        REAL hpvs=hpv1[j]*hpsk;
        REAL hpss=hps1[j]*hpsk;
        REAL hpssv=hpss1[j]*hpvk;
        REAL hpvss=hpv1[j]*hpssk;
        for (int i=2;i<n_pols_1d;i++) {
          Tensor<2,3> &ggv=gg[iv++];
          REAL hpvi=hpv0[i];
          REAL hpsi=hps0[i];
          ggv[0][0]=hpss0[i]*hpvv;
          ggv[1][0]=hpsi*hpsv;
          ggv[2][0]=hpsi*hpvs;
          ggv[1][1]=hpvi*hpssv;
          ggv[2][1]=hpvi*hpss;
          ggv[2][2]=hpvi*hpvss;
          ggv[0][1]=ggv[1][0];
          ggv[0][2]=ggv[2][0];
          ggv[1][2]=ggv[2][1];
        }
      }
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalHexahedronPolynomials::shapeFunctionsOrder(
const Shape *element,NumPtr<int> &renumber) const {
  CHECK_POINTER(dynamic_cast<const Hexahedron*>(element));
//since requiresShapeToComputeValues=true, no reordering is needed
  renumber.allocate(getNumber());
  for (int i=0;i<renumber.getNumber();i++) renumber[i]=i;
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalHexahedronPolynomials::nonzeroShapesOnFace(
const Shape1 *e,int f,NumPtr<int> &indices) const {
  const Shape3 *e3=dynamic_cast<const Shape3*>(e);
  int order=n_pols_1d-1;
  int nv=e->numberVertices();
  int nvof=e3->numberVerticesOnFace(f);
  int lpf=linesPerFace(nv,f);
  int nieslpof=numberInteriorEquallySpacedLatticePointsOnFace(nv,f,order);
  indices.allocate(nvof+lpf*(order-1)+nieslpof);
  int iv=0;
  for (int v=0;v<nvof;v++) indices[iv++]=e->elementVertexFromFace(f,v);
  for (int l=0;l<lpf;l++) {
    int offset=nv+elementLineFromFace(nv,f,l)*(order-1)-1;
    for (int i=1;i<order;i++) indices[iv++]=offset+i;
  }
  int offset=nv+e3->numberCodimension2Shapes()*(order-1);
  for (int i=0;i<f;i++) {
    offset+=numberInteriorEquallySpacedLatticePointsOnFace(nv,i,order);
  }
  for (int i=0;i<nieslpof;i++) indices[iv++]=offset+i;
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void HierarchicalHexahedronPolynomials::printOn(ostream &os) const {
  os << "HierarchicalHexahedronPolynomials: n_pols_1d = " << n_pols_1d
     << "\n\tn_pols = " << n_pols << endl;
  ShapeFunction<3>::printOn(os);
}
#endif

template int __cmath_power<int>(int,unsigned int);
template class ShapeFunction<1>;
template class TensorProductPolynomials<1>;
#if (SPACEDIM>1)
template class ShapeFunction<2>;
template class TensorProductPolynomials<2>;
template class HomogeneousPolynomials<SPACEDIM>;
#endif
#if (SPACEDIM==3)
template class ShapeFunction<3>;
template class TensorProductPolynomials<3>;
template class HomogeneousPolynomials<2>;
#endif

#include "NumPtr.C"
INSTANTIATE_NUMPTR(const Polynomial*);
template class NumPtrBase<NumPtr<NumPtr<REAL> >,int>;
template class NumPtr<NumPtr<NumPtr<REAL> > >;
template class NumPtrBase<Tensor<1,SPACEDIM>,int>;
template class NumPtr<Tensor<1,SPACEDIM> >;
template class NumPtrBase<Tensor<2,SPACEDIM>,int>;
template class NumPtr<Tensor<2,SPACEDIM> >;

template class ArrayBase<2,NumPtr<REAL> >;
template class Array<2,NumPtr<REAL> >;
#if (SPACEDIM==3)
template class NumPtrBase<Tensor<1,2>,int>;
template class NumPtr<Tensor<1,2> >;
template class NumPtrBase<Tensor<2,2>,int>;
template class NumPtr<Tensor<2,2> >;
#endif
