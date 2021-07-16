#ifndef SHAPE_FUNCTION_H
#define SHAPE_FUNCTION_H "$Header:$"
//----------------------  tensor_product_polynomials.h  -------------
//   $Id: tensor_product_polynomials.h,v 1.24 2003/05/05 13:50:14 wolf Exp $
//   Version: $Name: Version-4-0-0 $
//
//   Copyright (C) 2000, 2001, 2002, 2003 by the deal.II authors
//
//   This file is subject to QPL and may not be  distributed
//   without copyright and license information. Please refer
//   to the file deal.II/doc/license.html for the  text  and
//   further information on this license.
//
//----------------------  tensor_product_polynomials.h  -------------
//
//modified from deal.II/base/include/base/tensor_product_polynomials.h
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

#include <DDecl.H>
#include <Point.H>
#include <Polynomial.H>
#include <Shape.H>
#include <Tensor.H>
#include <Types.H>

//not in DEAL
template<int dim> class ShapeFunction {
  public:
    ShapeFunction() {}
    virtual ~ShapeFunction() {}
    virtual bool requiresShapeToComputeValues() const { return false; }
    virtual void computeValues(const Point<dim>&,NumPtr<REAL>&) const=0;
    virtual void computeGrads(const Point<dim>&,NumPtr<Tensor<1,dim> >&) 
      const=0;
    virtual void computeGradGrads(const Point<dim>&,
      NumPtr<Tensor<2,dim> >&) const=0;
    virtual void computeValuesFor(const Shape*,const Point<dim> &p,
    NumPtr<REAL> &v) const {
      OBSOLETE("does not require shape");
    }
    virtual void computeGradsFor(const Shape*,const Point<dim> &p,
    NumPtr<Tensor<1,dim> > &g) const {
      OBSOLETE("does not require shape");
    }
    virtual void computeGradGradsFor(const Shape*,const Point<dim> &p,
    NumPtr<Tensor<2,dim> > &gg) const {
      OBSOLETE("does not require shape");
    }
//  renumber[fe_s] = dof_s:
    virtual void shapeFunctionsOrder(const Shape*,NumPtr<int>&) const=0;
    virtual void nonzeroShapesOnFace(const Shape1*,int,NumPtr<int>&) 
      const=0;
    virtual const Shape* reorientChildForParent(const Shape *parent,
    const Shape *child,int child_index,
    const NumPtr<int> &line_interior_child_vertices) const {
      return child->
        reorientFor(parent,child_index,line_interior_child_vertices);
    }
    virtual int getNumber() const=0;
    virtual void printOn(ostream&) const;
};

//compare to DEAL TensorProductPolynomials
template<int dim> class TensorProductPolynomials : 
public ShapeFunction<dim> {
  protected:
    const Polynomial *polynomials[dim];
    int n_tensor_pols;
    int n_pols[dim];
    bool delete_all;
  public:
    template <class Pol> TensorProductPolynomials(const Pol *pols,
    int degree) : delete_all(false),n_tensor_pols(1) {
      for (int d=0;d<dim;d++) { 
        polynomials[d]=pols;
        n_pols[d]=degree+1;
        n_tensor_pols*=n_pols[d];
      }
    }
    template <class Pol> TensorProductPolynomials(const Pol *pols0,
    const Pol *pols1,int d0,int d1) : delete_all(true) {
      CHECK_SAME(dim,2);
      polynomials[0]=pols0;
      polynomials[1]=pols1;
      n_pols[0]=d0+1;
      n_pols[1]=d1+1;
      n_tensor_pols=n_pols[0]*n_pols[1];
    }
    template <class Pol> TensorProductPolynomials(const Pol *pols0,
    const Pol *pols1,const Pol *pols2,int d0,int d1,int d2) : 
    delete_all(true) {
      CHECK_SAME(dim,3);
      polynomials[0]=pols0;
      polynomials[1]=pols1;
      polynomials[2]=pols2;
      n_pols[0]=d0+1;
      n_pols[1]=d1+1;
      n_pols[2]=d2+1;
      n_tensor_pols=n_pols[0]*n_pols[1]*n_pols[2];
    }
    ~TensorProductPolynomials() {
      if (delete_all) {
        for (int d=0;d<dim;d++) {
          delete polynomials[d]; polynomials[d]=0;
        }
      } else {
        delete polynomials[0]; polynomials[0]=0;
      }
    }
    void computeValues(const Point<dim>&,NumPtr<REAL>&) const;
    void computeGrads(const Point<dim>&,NumPtr<Tensor<1,dim> >&) const;
    void computeGradGrads(const Point<dim>&,NumPtr<Tensor<2,dim> >&) const;
    int getNumber() const { return n_tensor_pols; }
    void printOn(ostream&) const;
};

#if (SPACEDIM>1)
//compare to DEAL AnisotropicPolynomials
class AnisotropicPolynomials : public ShapeFunction<SPACEDIM> {
  private:
    NumPtr<const Polynomial*> polynomials;
    int n_tensor_pols;
    int n_pols[SPACEDIM];
  public:
    AnisotropicPolynomials(const NumPtr<const Polynomial*>&,
      const int (&)[SPACEDIM]);
    virtual ~AnisotropicPolynomials() {}
    void computeValues(const Point<SPACEDIM>&,NumPtr<REAL>&) const;
    void computeGrads(const Point<SPACEDIM>&,
      NumPtr<Tensor<1,SPACEDIM> >&) const;
    void computeGradGrads(const Point<SPACEDIM>&,
      NumPtr<Tensor<2,SPACEDIM> >&) const;
    virtual void shapeFunctionsOrder(const Shape*,NumPtr<int>&) const {
      OBSOLETE("not programmed");
    }
    void nonzeroShapesOnFace(const Shape1 *e,int f,NumPtr<int> &indices) 
    const {
      OBSOLETE("not programmed");
    }
    virtual const Shape* reorientChildForParent(const Shape *parent,
    const Shape *child,int child_index,
    const NumPtr<int> &line_interior_child_vertices) const {
      OBSOLETE("not programmed");
    }
    int getNumber() const { return n_tensor_pols; }
    void printOn(ostream&) const;
};
#endif

#if (SPACEDIM>1)
//not in DEAL
class TrianglePolynomials : public ShapeFunction<2> {
  protected:
    int order;
  public:
    TrianglePolynomials(int ord) : order(ord) {}
    ~TrianglePolynomials() {}
    int getNumber() const { return ((order+1)*(order+2))/2; }
    void printOn(ostream&) const;
};

//not in DEAL
template<int dim> class HomogeneousPolynomials : public ShapeFunction<dim>
{
  protected:
    int order;
    int n_pols;
  public:
    HomogeneousPolynomials(int degree) : order(degree),n_pols(1) {
      CHECK_TEST(dim>=2 && dim<=3);
      for (int d=1;d<dim;d++) n_pols*=degree+d;
      if (dim>2) n_pols/=2;
    }
    virtual ~HomogeneousPolynomials() { }
    void computeValues(const Point<dim>&,NumPtr<REAL>&) const;
    void computeGrads(const Point<dim>&,NumPtr<Tensor<1,dim> >&) const;
    void computeGradGrads(const Point<dim>&,NumPtr<Tensor<2,dim> >&) const;
    virtual void shapeFunctionsOrder(const Shape*,NumPtr<int>&) const {
      OBSOLETE("why is this called?");
    }
    virtual void nonzeroShapesOnFace(const Shape1*,int,NumPtr<int>&) 
    const {
      OBSOLETE("why is this called?");
    }
    virtual int getNumber() const { return n_pols; }
    virtual void printOn(ostream&) const;
};

//not in DEAL
class HomogeneousLagrangePolynomials2 : public ShapeFunction<2> {
  protected:
    NumPtr<REAL> zeros;
    int order;
    int n_pols;
  public: 
    HomogeneousLagrangePolynomials2(int);
    virtual ~HomogeneousLagrangePolynomials2() {}
    void computeValues(const Point<2>&,NumPtr<REAL>&) const;
    void computeGrads(const Point<2>&,NumPtr<Tensor<1,2> >&) const;
    void computeGradGrads(const Point<2>&,NumPtr<Tensor<2,2> >&) const;
    virtual void shapeFunctionsOrder(const Shape*,NumPtr<int>&) const {
      OBSOLETE("why is this called?");
    }
    virtual void nonzeroShapesOnFace(const Shape1*,int,NumPtr<int>&) 
    const {
      OBSOLETE("why is this called?");
    }
    virtual int getNumber() const { return n_pols; }
    virtual void printOn(ostream&) const;
};

# if (SPACEDIM==3)
//not in DEAL
class HomogeneousLagrangePolynomials3 : public ShapeFunction<3> {
  protected:
    NumPtr<REAL> zeros;
    int order;
    int n_pols;
  public: 
    HomogeneousLagrangePolynomials3(int);
    virtual ~HomogeneousLagrangePolynomials3() {}
    void computeValues(const Point<3>&,NumPtr<REAL>&) const;
    void computeGrads(const Point<3>&,NumPtr<Tensor<1,3> >&) const;
    void computeGradGrads(const Point<3>&,NumPtr<Tensor<2,3> >&) const;
    virtual void shapeFunctionsOrder(const Shape*,NumPtr<int>&) const {
      OBSOLETE("why is this called?");
    }
    virtual void nonzeroShapesOnFace(const Shape1*,int,NumPtr<int>&) 
    const {
      OBSOLETE("why is this called?");
    }
    virtual int getNumber() const { return n_pols; }
    virtual void printOn(ostream&) const;
};
# endif
#endif

#if (SPACEDIM==3)
//not in DEAL
class TetrahedronPolynomials : public ShapeFunction<3> {
  protected:
    int order;
  public:
    TetrahedronPolynomials(int ord) : order(ord) {}
    ~TetrahedronPolynomials() {}
    int getNumber() const { 
      return ((order+1)*(order+2)*(order+3))/6; 
    }
    void printOn(ostream&) const;
};
#endif

#if (SPACEDIM==3)
//not in DEAL
class PrismPolynomials : public ShapeFunction<3> {
  protected:
    ShapeFunction<2> *tp;
    Polynomial *lp;
    int order01;
    int order2;
  public:
    PrismPolynomials(ShapeFunction<2> *t,Polynomial *p,int ord01,
      int ord2) : tp(t),lp(p),order01(ord01),order2(ord2) {}
    ~PrismPolynomials() { delete tp; tp=0; delete lp; lp=0; }
    void computeValues(const Point<3>&,NumPtr<REAL>&) const;
    void computeGrads(const Point<3>&,NumPtr<Tensor<1,3> >&) const;
    void computeGradGrads(const Point<3>&,NumPtr<Tensor<2,3> >&) const;
    int getNumber() const { return tp->getNumber()*(order2+1); }
    void printOn(ostream&) const;
};
#endif

//not in DEAL
class C0LagrangeIntervalPolynomials : public TensorProductPolynomials<1> {
  public:
    C0LagrangeIntervalPolynomials(int ord) : 
      TensorProductPolynomials<1>(OPERATOR_NEW C0LagrangePolynomial(ord),
      ord) {}
    virtual ~C0LagrangeIntervalPolynomials() {}
    virtual void shapeFunctionsOrder(const Shape *element,
    NumPtr<int> &renumber) const {
      element->reorderEquallySpacedLatticePoints(n_pols[0]-1,renumber);
    }
    void nonzeroShapesOnFace(const Shape1 *e,int f,NumPtr<int> &indices)
    const {
      indices.allocate(1);
      indices[0]=e->equallySpacedLatticePointIndexOfVertex(f,n_pols[0]-1);
    }
    void printOn(ostream&) const;
};

#if (SPACEDIM>1)
//not in DEAL
class C0LagrangeTrianglePolynomials : public TrianglePolynomials {
  private:
    NumPtr<REAL> zeros;
  public:
    C0LagrangeTrianglePolynomials(int);
    virtual ~C0LagrangeTrianglePolynomials() {}
    void computeValues(const Point<2>&,NumPtr<REAL>&) const;
    void computeGrads(const Point<2>&,NumPtr<Tensor<1,2> >&) const;
    void computeGradGrads(const Point<2>&,NumPtr<Tensor<2,2> >&) const;
    void shapeFunctionsOrder(const Shape*,NumPtr<int>&) const;
    void nonzeroShapesOnFace(const Shape1 *e,int f,NumPtr<int> &indices) 
    const {
      dynamic_cast<const Shape2*>(e)->
        equallySpacedLatticePointIndicesOfFace(f,order,indices);
    }
    void printOn(ostream&) const;
};
#endif

#if (SPACEDIM>1)
//not in DEAL
class C0LagrangeQuadrilateralPolynomials : 
public TensorProductPolynomials<2> {
  public:
    C0LagrangeQuadrilateralPolynomials(int ord) : 
      TensorProductPolynomials<2>(OPERATOR_NEW C0LagrangePolynomial(ord),
      ord) {}
    C0LagrangeQuadrilateralPolynomials(int ord0,int ord1) : 
      TensorProductPolynomials<2>(OPERATOR_NEW C0LagrangePolynomial(ord0),
        OPERATOR_NEW C0LagrangePolynomial(ord1),ord0,ord1) {}
    virtual ~C0LagrangeQuadrilateralPolynomials() {}
    virtual void shapeFunctionsOrder(const Shape *element,
    NumPtr<int> &renumber) const {
      dynamic_cast<const Quadrilateral*>(element)->
      reorderEquallySpacedLatticePoints(n_pols[0]-1,n_pols[1]-1,renumber);
    }
    void nonzeroShapesOnFace(const Shape1 *e,int f,NumPtr<int> &indices)
    const {
      dynamic_cast<const Quadrilateral*>(e)->
        equallySpacedLatticePointIndicesOfFace(f,n_pols[0]-1,n_pols[1]-1,
        indices);
      return;
    }
    void printOn(ostream&) const;
};
#endif

#if (SPACEDIM==3)
//not in DEAL
class C0LagrangeTetrahedronPolynomials : public TetrahedronPolynomials {
  private:
    NumPtr<REAL> zeros;
  public:
    C0LagrangeTetrahedronPolynomials(int);
    virtual ~C0LagrangeTetrahedronPolynomials() {}
    void computeValues(const Point<3>&,NumPtr<REAL>&) const;
    void computeGrads(const Point<3>&,NumPtr<Tensor<1,3> >&) const;
    void computeGradGrads(const Point<3>&,NumPtr<Tensor<2,3> >&) const;
    void shapeFunctionsOrder(const Shape*,NumPtr<int>&) const;
    void nonzeroShapesOnFace(const Shape1 *e,int f,NumPtr<int> &indices) 
    const {
      dynamic_cast<const Shape2*>(e)->
        equallySpacedLatticePointIndicesOfFace(f,order,indices);
    }
    void printOn(ostream&) const;
};
#endif

#if (SPACEDIM==3)
//not in DEAL
class C0LagrangePrismPolynomials : public PrismPolynomials {
  public:
    C0LagrangePrismPolynomials(int ord) : 
      PrismPolynomials(OPERATOR_NEW C0LagrangeTrianglePolynomials(ord),
      OPERATOR_NEW C0LagrangePolynomial(ord),ord,ord) {}
    C0LagrangePrismPolynomials(int ord01,int ord2) : 
      PrismPolynomials(OPERATOR_NEW C0LagrangeTrianglePolynomials(ord01),
      OPERATOR_NEW C0LagrangePolynomial(ord2),ord01,ord2) {}
    virtual ~C0LagrangePrismPolynomials() {}
    virtual void shapeFunctionsOrder(const Shape*,NumPtr<int>&) const;
    void nonzeroShapesOnFace(const Shape1 *e,int f,NumPtr<int> &indices) 
    const {
      CHECK_SAME(order01,order2);
      dynamic_cast<const Shape2*>(e)->
        equallySpacedLatticePointIndicesOfFace(f,order01,indices);
    }
    void printOn(ostream&) const;
};
#endif

#if (SPACEDIM==3)
//not in DEAL
class C0LagrangeHexahedronPolynomials : 
public TensorProductPolynomials<3> {
  public:
    C0LagrangeHexahedronPolynomials(int ord) : 
      TensorProductPolynomials<3>(OPERATOR_NEW C0LagrangePolynomial(ord),
      ord) {}
    C0LagrangeHexahedronPolynomials(int ord0,int ord1,int ord2) : 
      TensorProductPolynomials<3>(OPERATOR_NEW C0LagrangePolynomial(ord0),
      OPERATOR_NEW C0LagrangePolynomial(ord1),
      OPERATOR_NEW C0LagrangePolynomial(ord2),ord0,ord1,ord2) {}
    virtual ~C0LagrangeHexahedronPolynomials() {}
    virtual void shapeFunctionsOrder(const Shape *element,
    NumPtr<int> &renumber) const {
      dynamic_cast<const Hexahedron*>(element)->
      reorderEquallySpacedLatticePoints(n_pols[0]-1,n_pols[1]-1,
        n_pols[2]-1,renumber);
    }
    void nonzeroShapesOnFace(const Shape1 *e,int f,NumPtr<int> &indices)
    const {
      dynamic_cast<const Hexahedron*>(e)->
        equallySpacedLatticePointIndicesOfFace(f,n_pols[0]-1,
        n_pols[1]-1,n_pols[2]-1,indices);
      return;
    }
    void printOn(ostream&) const;
};
#endif

//not in DEAL
class HierarchicalIntervalPolynomials : public ShapeFunction<1> {
  private:
    int n_pols;
  public:
    HierarchicalIntervalPolynomials(int ord) : n_pols(ord+1) {}
    ~HierarchicalIntervalPolynomials() {}
    void computeValues(const Point<1>&,NumPtr<REAL>&) const;
    void computeGrads(const Point<1>&,NumPtr<Tensor<1,1> >&) const;
    void computeGradGrads(const Point<1>&,NumPtr<Tensor<2,1> >&) const;
    void shapeFunctionsOrder(const Shape*,NumPtr<int>&) const;
    void nonzeroShapesOnFace(const Shape1*,int f,NumPtr<int> &indices)
    const {
      indices.allocate(1);
      indices[0]=f;
    }
    virtual const Shape* reorientChildForParent(const Shape *parent,
    const Shape *child,int child_index,
    const NumPtr<int> &line_interior_child_vertices) const {
      return child;
    }
    int getNumber() const { return n_pols; }
    void printOn(ostream&) const;
};

#if (SPACEDIM>1)
//not in DEAL
class HierarchicalTrianglePolynomials : public TrianglePolynomials {
  public:
    HierarchicalTrianglePolynomials(int ord) : TrianglePolynomials(ord) {}
    ~HierarchicalTrianglePolynomials() {}
    void computeValues(const Point<2>&,NumPtr<REAL>&) const;
    void computeGrads(const Point<2>&,NumPtr<Tensor<1,2> >&) const;
    void computeGradGrads(const Point<2>&,NumPtr<Tensor<2,2> >&) const;
    void shapeFunctionsOrder(const Shape*,NumPtr<int>&) const;
    void nonzeroShapesOnFace(const Shape1*,int,NumPtr<int>&) const;
    virtual const Shape* reorientChildForParent(const Shape *parent,
    const Shape *child,int child_index,
    const NumPtr<int> &line_interior_child_vertices) const {
      return child;
    }
    void printOn(ostream&) const;
};
#endif

#if (SPACEDIM>1)
//not in DEAL
class HierarchicalQuadrilateralPolynomials : public ShapeFunction<2> { 
  private:
    int n_pols_1d;
    int n_pols;
  public:
    HierarchicalQuadrilateralPolynomials(int ord) : n_pols_1d(ord+1),
      n_pols(__cmath_power(ord+1,2)) { }
    ~HierarchicalQuadrilateralPolynomials() {}
    void computeValues(const Point<2>&,NumPtr<REAL>&) const;
    void computeGrads(const Point<2>&,NumPtr<Tensor<1,2> >&) const;
    void computeGradGrads(const Point<2>&,NumPtr<Tensor<2,2> >&) const;
    void shapeFunctionsOrder(const Shape*,NumPtr<int>&) const;
    void nonzeroShapesOnFace(const Shape1*,int,NumPtr<int>&) const;
    virtual const Shape* reorientChildForParent(const Shape *parent,
    const Shape *child,int child_index,
    const NumPtr<int> &line_interior_child_vertices) const {
      return child;
    }
    int getNumber() const { return n_pols; }
    void printOn(ostream&) const;
};
#endif

#if (SPACEDIM==3)
//not in DEAL
class HierarchicalTetrahedronPolynomials : public TetrahedronPolynomials {
  public:
    HierarchicalTetrahedronPolynomials(int ord) :
      TetrahedronPolynomials(ord) {}
    virtual ~HierarchicalTetrahedronPolynomials() {}
    virtual bool requiresShapeToComputeValues() const { 
      return true;
    }
    void computeValues(const Point<3>&,NumPtr<REAL>&) const;
    void computeGrads(const Point<3>&,NumPtr<Tensor<1,3> >&) const;
    void computeGradGrads(const Point<3>&,NumPtr<Tensor<2,3> >&) const;
    void computeValuesFor(const Shape*,const Point<SPACEDIM>&,
      NumPtr<REAL>&) const;
    void computeGradsFor(const Shape*,const Point<SPACEDIM>&,
      NumPtr<Tensor<1,SPACEDIM> >&) const;
    void computeGradGradsFor(const Shape*,const Point<SPACEDIM>&,
      NumPtr<Tensor<2,SPACEDIM> >&) const;
    void shapeFunctionsOrder(const Shape*,NumPtr<int>&) const;
    void nonzeroShapesOnFace(const Shape1*,int,NumPtr<int>&) const;
    virtual const Shape* reorientChildForParent(const Shape *parent,
    const Shape *child,int child_index,
    const NumPtr<int> &line_interior_child_vertices) const {
      return child;
    }
    void printOn(ostream&) const;
};
#endif

#if (SPACEDIM==3)
//not in DEAL
class HierarchicalPrismPolynomials : public PrismPolynomials {
  public:
    HierarchicalPrismPolynomials(int ord) : 
      PrismPolynomials(OPERATOR_NEW HierarchicalTrianglePolynomials(ord),
      OPERATOR_NEW HierarchicalPolynomial,ord,ord) {}
    virtual ~HierarchicalPrismPolynomials() {}
    virtual bool requiresShapeToComputeValues() const { 
      return true;
    }
    void computeValues(const Point<3>&,NumPtr<REAL>&) const;
    void computeGrads(const Point<3>&,NumPtr<Tensor<1,3> >&) const;
    void computeGradGrads(const Point<3>&,NumPtr<Tensor<2,3> >&) const;
    void computeValuesFor(const Shape*,const Point<SPACEDIM>&,
      NumPtr<REAL>&) const;
    void computeGradsFor(const Shape*,const Point<SPACEDIM>&,
      NumPtr<Tensor<1,SPACEDIM> >&) const;
    void computeGradGradsFor(const Shape*,const Point<SPACEDIM>&,
      NumPtr<Tensor<2,SPACEDIM> >&) const;
    void shapeFunctionsOrder(const Shape*,NumPtr<int>&) const;
    void nonzeroShapesOnFace(const Shape1*,int,NumPtr<int>&) const;
    virtual const Shape* reorientChildForParent(const Shape *parent,
    const Shape *child,int child_index,
    const NumPtr<int> &line_interior_child_vertices) const {
      return child;
    }
    void printOn(ostream&) const;
};
#endif

#if (SPACEDIM==3)
//not in DEAL
class HierarchicalHexahedronPolynomials : public ShapeFunction<3> {
  private:
    int n_pols_1d;
    int n_pols;
  public:
    HierarchicalHexahedronPolynomials(int ord) : n_pols_1d(ord+1), 
      n_pols(__cmath_power(ord+1,3)) { }
    virtual ~HierarchicalHexahedronPolynomials() {}
    virtual bool requiresShapeToComputeValues() const { 
      return true;
    }
    void computeValues(const Point<3>&,NumPtr<REAL>&) const;
    void computeGrads(const Point<3>&,NumPtr<Tensor<1,3> >&) const;
    void computeGradGrads(const Point<3>&,NumPtr<Tensor<2,3> >&) const;
    void computeValuesFor(const Shape*,const Point<SPACEDIM>&,
      NumPtr<REAL>&) const;
    void computeGradsFor(const Shape*,const Point<SPACEDIM>&,
      NumPtr<Tensor<1,SPACEDIM> >&) const;
    void computeGradGradsFor(const Shape*,const Point<SPACEDIM>&,
      NumPtr<Tensor<2,SPACEDIM> >&) const;
    void shapeFunctionsOrder(const Shape*,NumPtr<int>&) const;
    void nonzeroShapesOnFace(const Shape1*,int,NumPtr<int>&) const;

    virtual const Shape* reorientChildForParent(const Shape *parent,
    const Shape *child,int child_index,
    const NumPtr<int> &line_interior_child_vertices) const {
      return child;
    }
    int getNumber() const { return n_pols; }
    void printOn(ostream&) const;
};
#endif

#endif
