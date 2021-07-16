// "$Header:$"
//----------------------------  solver_control.cc  ------------------------
//  $Id: solver_control.cc,v 1.30 2003/01/08 17:58:16 wolf Exp $
//  Version: $Name: Version-4-0-0 $
//
//  Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//  This file is subject to QPL and may not be  distributed
//  without copyright and license information. Please refer
//  to the file deal.II/doc/license.html for the  text  and
//  further information on this license.
//
//----------------------------  solver_control.cc  ------------------------
//----------------------------  solver.h  ---------------------------
//  $Id: solver.h,v 1.29 2003/02/10 09:10:08 guido Exp $
//  Version: $Name: Version-4-0-0 $
//
//  Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//  This file is subject to QPL and may not be  distributed
//  without copyright and license information. Please refer
//  to the file deal.II/doc/license.html for the  text  and
//  further information on this license.
//
//----------------------------  solver.h  ---------------------------
//-------------------------  solver_bicgstab.h  ---------------------------
//  $Id: solver_bicgstab.h,v 1.47 2003/01/08 17:58:14 wolf Exp $
//  Version: $Name: Version-4-0-0 $
//
//  Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//  This file is subject to QPL and may not be  distributed
//  without copyright and license information. Please refer
//  to the file deal.II/doc/license.html for the  text  and
//  further information on this license.
//
//-------------------------  solver_bicgstab.h  ---------------------------
//----------------------------  solver_cg.h  ---------------------------
//  $Id: solver_cg.h,v 1.48 2003/02/14 16:35:45 wolf Exp $
//  Version: $Name: Version-4-0-0 $
//
//  Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//  This file is subject to QPL and may not be  distributed
//  without copyright and license information. Please refer
//  to the file deal.II/doc/license.html for the  text  and
//  further information on this license.
//
//----------------------------  solver_cg.h  ---------------------------
//----------------------------  solver_gmres.h  ---------------------------
//  $Id: solver_gmres.h,v 1.64 2003/05/18 02:26:47 wolf Exp $
//  Version: $Name: Version-4-0-0 $
//
//  Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//  This file is subject to QPL and may not be  distributed
//  without copyright and license information. Please refer
//  to the file deal.II/doc/license.html for the  text  and
//  further information on this license.
//
//----------------------------  solver_gmres.h  ---------------------------
//---------------------------  solver_minres.h  ---------------------------
//  $Id: solver_minres.h,v 1.22 2003/04/26 16:10:31 wolf Exp $
//  Version: $Name: Version-4-0-0 $
//
//  Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//  This file is subject to QPL and may not be  distributed
//  without copyright and license information. Please refer
//  to the file deal.II/doc/license.html for the  text  and
//  further information on this license.
//
//----------------------------  solver_minres.h  -----------------------
//----------------------------  solver_qmrs.h  ---------------------------
//  $Id: solver_qmrs.h,v 1.31 2003/04/30 23:07:00 wolf Exp $
//  Version: $Name: Version-4-0-0 $
//
//  Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//  This file is subject to QPL and may not be  distributed
//  without copyright and license information. Please refer
//  to the file deal.II/doc/license.html for the  text  and
//  further information on this license.
//
//----------------------------  solver_qmrs.h  ---------------------------
//-----------------------  solver_richardson.h  ---------------------------
//   $Id: solver_richardson.h,v 1.35 2003/04/26 16:11:11 wolf Exp $
//   Version: $Name: Version-4-0-0 $
//
//   Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//   This file is subject to QPL and may not be  distributed
//   without copyright and license information. Please refer
//   to the file deal.II/doc/license.html for the  text  and
//   further information on this license.
//
//------------------------  solver_richardson.h  --------------------------
//
//modified from deal.II/lac/source/solver_control.h and
//  deal.II/lac/include/lac/solver_control, solver.h,
//  solver_bicgstab.h, solver_cg.h, solver_gmres.h, solver_minres.h,
//  solver_qmrs.h and solver_richardson.h
//  by John Trangenstein, August 2009
#include "Solver.H"
#include <cmath>
#include <sstream>
#include "Precondition.H"

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
extern "C" {
void F77NAME(daxpy)(const int &n,const double &da,const double *dx,
  const int &incx,double *dy,const int &incy);
void F77NAME(dcopy)(const int &n,const double *sx,const int &incx,
                    double *sy,const  int &incy);
double F77NAME(ddot)(const int &n,const double *dx,const int &incx,
  const double *dy,const int &incy);
void F77NAME(dhels)(const double *a,const int &lda,const int &n,
  const double *q,double *b);
void F77NAME(dheqr)(double *a,const int &lda,const int &n,double *q,
  int &info,const int &ijob);
double F77NAME(dnrm2)(const int &n,const double *x,const int &incx);
void F77NAME(dorth)(double *vnew,const double *v,double *hes,
  const int &n,const int &ll,const int &ldhes,const int &kmp,
  double &snormw);
void F77NAME(drot)(int &n,double *sx,int &incx,double *sy,int &incy,
                   const double &c,const double &s);
void F77NAME(drotg)(const double *sa,const double *sb,double &c,
  double &s);
void F77NAME(dscal)(const int &n,const double &da,double *dx,
  const int &incs);
void F77NAME(dtrsm)(const char &side,const char &uplo,
  const char &transa,const char &diag,int *m,int *n,double *alpha,
  const double *a,int *lda,double *b,int *ldb);
};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ostream& operator<<(ostream &os,SOLVER_STATE ss) {
  switch (ss) {
    case SOLVER_ITERATE:
      os << "SOLVER_ITERATE";
      break;
    case SOLVER_SUCCESS:
      os << "SOLVER_SUCCESS";
      break;
    case SOLVER_FAILURE:
      os << "SOLVER_FAILURE";
      break;
    default:
      os << "unknown";
      break;
  }
  return os;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
SOLVER_STATE SolverControl::check(int step,double check_value) {
  if (m_log_history && ((step % m_log_frequency)==0)) {
    cout << "   Check " << step << "\t" << check_value << endl;
  }
  if (step==0) {
    if (check_failure) {
      failure_residual=relative_failure_residual*check_value;
    }
    if (m_log_result) {
      cout << "   Starting value " << check_value << endl;
    }
  }
  if (check_value<=tol) {
    if (m_log_result) {
      cout << "   Convergence step " << step
           << " value " << check_value << endl;
    }
    return SOLVER_SUCCESS;
  }
  if (step>=maxsteps || isnan(check_value) ||
  (check_failure && (check_value > failure_residual))) {
    if (m_log_result) {
      cout << "   Failure step " << step << " value " << check_value 
           << " goal = " << tol << endl;
    }
//  lcheck=SOLVER_FAILURE;
    return SOLVER_FAILURE;
  }
//lcheck=SOLVER_ITERATE;
  return SOLVER_ITERATE;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SolverControl::printOn(ostream &os) const {
  os << "SolverControl:" << endl;
  os << "\tmaxsteps = " << maxsteps
     << "\n\ttol = " << tol
//   << "\n\tlcheck = " << lcheck
//   << "\n\tlvalue = " << lvalue
//   << "\n\tlstep = " << lstep
     << "\n\tcheck_failure = " << check_failure
     << "\n\trelative_failure_residual = " << relative_failure_residual
     << "\n\tfailure_residual = " << failure_residual
//   << "\n\tm_log_history = " << m_log_history
//   << "\n\tm_log_frequency = " << m_log_frequency
//   << "\n\tm_log_result = " << m_log_result
     << endl;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
SOLVER_STATE ReductionControl::check(int step,double check_value) {
  if (step==0) {
    if (m_log_result) { 
      cout << "   Starting value " << check_value << endl;
    }
    initial_val=check_value;
    return SOLVER_ITERATE;
  } else {
    if (m_log_history && ((step % m_log_frequency)==0)) {
      cout << "   Check " << step << "\t" << check_value << endl;
    }
    if (check_value<=initial_val*tolerance()) {
      if (m_log_result) {
        cout << "Convergence step " << step << " value " << check_value 
             << endl;
      }
      return SOLVER_SUCCESS;
    } 
    if (step>=maxsteps || isnan(check_value)) {
      if (m_log_result) {
        cout << "   Failure step " << step << " value " << check_value
             << " goal = " << initial_val*tolerance() << endl;
      }
      return SOLVER_FAILURE;
    }
    return SOLVER_ITERATE;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void ReductionControl::printOn(ostream &os) const {
  os << "ReductionControl:" << endl;
  os << "\tinitial_val = " << initial_val
     << endl;
  SolverControl::printOn(os);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
SOLVER_STATE RelativeErrorControl::check(int step,double check_value) {
  if (step==0) {
    if (m_log_result) {
      cout << "   Starting value " << check_value << endl;
    }
//  rhs_norm is set in the solvers by calling setRhsNorm
  } else {
    if (m_log_history && ((step % m_log_frequency)==0)) {
      cout << "   Check " << step << "\t" << check_value << endl;
    }
    if (check_value<=tol*rhs_norm) {
      if (m_log_result) {
        cout << "   Convergence step " << step
             << " value " << check_value << endl;
      }
      return SOLVER_SUCCESS;
    }
    if (step>=maxsteps || isnan(check_value)) {
      if (m_log_result) {
        cout << "   Failure step " << step << " value " << check_value 
             << " goal " << tol*rhs_norm << endl;
      }
      return SOLVER_FAILURE;
    }
  }
  return SOLVER_ITERATE;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void RelativeErrorControl::printOn(ostream &os) const {
  os << "RelativeErrorControl:" << endl;
  os << "\trhs_norm = " << rhs_norm << endl;
  SolverControl::printOn(os);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SolverPreconditioner::printOn(ostream &os) const {
  os << "SolverPreconditioner:" << endl;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Solver::printOn(ostream &os) const {
  os << "Solver:" << endl;
  os << "\tcntrl:" << endl;
  cntrl.printOn(os);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void BicgstabSolver::solve(const SolverMatrix &A,Vector &x,
const Vector &b,const SolverPreconditioner &precondition) {
//TRACER_CALL(t,"BicgstabSolver::solve");
  RelativeErrorControl *rec=
    dynamic_cast<RelativeErrorControl*>(&control());
  if (rec!=0) rec->setRhsNorm(b.linftyNorm());
  Vector *Vr=x.clone();
  Vector *Vrbar=x.clone();
  Vector *Vp=x.clone();
  Vector *Vy=x.clone();
  Vector *Vz=x.clone();
  Vector *Vs=x.clone();
  Vector *Vt=x.clone();
  Vector *Vv=x.clone();

  int it=0;
  bool keep_going=true;
  SOLVER_STATE state=SOLVER_ITERATE;
  do {
    A.vmult(*Vr,x); 
    Vr->scale(-1.);
    Vr->add(1.,b);
    res=Vr->linftyNorm();
    *Vp=0.;
    *Vv=0;
    *Vrbar=*Vr;
    if (control().check(it,res)==SOLVER_SUCCESS) break;

    state=SOLVER_ITERATE;
    alpha=omega=rho=1.;
    keep_going=false;
    do {
      it++;
      rhobar=(*Vr)*(*Vrbar);
      beta=rhobar*alpha/(rho*omega);
      rho=rhobar;
      Vp->scale(beta);
      Vp->add(1.,*Vr);
      Vp->add(-beta*omega,*Vv);
      precondition.vmult(*Vy,*Vp);
      A.vmult(*Vv,*Vy);
      rhobar=(*Vrbar)*(*Vv);
      alpha=rho/rhobar;
      if (abs(alpha)>breakdown) {
        keep_going=true;
        break;
      }
      Vs->equ(1.,*Vr);
      Vs->add(-alpha,*Vv);
      precondition.vmult(*Vz,*Vs);
      A.vmult(*Vt,*Vz);
      rhobar=(*Vt)*(*Vs);
      double vtnorm=(*Vt)*(*Vt);
      omega=(vtnorm>0. ? rhobar/vtnorm : 1.);
      x.add(alpha,*Vy);
      x.add(omega,*Vz);
      Vr->equ(1.,*Vs);
      Vr->add(-omega,*Vt);
      if (exact_residual) res=criterion(A,x,b,*Vt);
      else res=Vr->linftyNorm();
      state=control().check(it,res);
    } while (state==SOLVER_ITERATE);
  } while (keep_going);
  if (Vr) delete Vr; Vr=0;
  if (Vrbar) delete Vrbar; Vrbar=0;
  if (Vp) delete Vp; Vp=0;
  if (Vy) delete Vy; Vy=0;
  if (Vz) delete Vz; Vz=0;
  if (Vs) delete Vs; Vs=0;
  if (Vt) delete Vt; Vt=0;
  if (Vv) delete Vv; Vv=0;
  CHECK_TEST(state==SOLVER_SUCCESS);
//cout << "\tnumber iterations = " << it << endl;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void BicgstabSolver::printOn(ostream &os) const {
  os << "BicgstabSolver:" << endl;
  os << "\talpha = " << alpha
     << "\n\tbeta = " << beta
     << "\n\tomega = " << omega
     << "\n\trho = " << rho
     << "\n\trhobar = " << rhobar
     << "\n\tres = " << res
     << "\n\texact_residual = " << exact_residual
     << "\n\tbreakdown = " << breakdown
     << endl;
  Solver::printOn(os);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void CGSolver::solve(const SolverMatrix &A,Vector &x,
const Vector &b,const SolverPreconditioner &precondition) {
//TRACER_CALL(t,"CGSolver::solve");
  RelativeErrorControl *rec=
    dynamic_cast<RelativeErrorControl*>(&control());
  if (rec!=0) rec->setRhsNorm(b.linftyNorm());
  Vector *g=x.clone();
  Vector *h=x.clone();
  Vector *d=x.clone();
  Vector *Ad=x.clone();
  if (!x.allZero()) {
    A.vmult(*g,x);
    g->scale(-1.);
    g->add(1.,b);
  } else *g=b;
  double res=g->linftyNorm();
//cout << "\tinitial residual norm = " << res << endl;
  SOLVER_STATE conv=control().check(0,res);
  int it=0;
  if (conv==SOLVER_ITERATE) {
    g->scale(-1.);
    precondition.vmult(*h,*g);
    d->equ(-1.,*h);
    double gh = (*g)*(*h);
    while (conv==SOLVER_ITERATE) {
      it++;
      A.vmult(*Ad,*d);
      double alpha=(*d)*(*Ad);
      alpha=gh/alpha;
      g->add(alpha,*Ad);
      x.add(alpha,*d);
      res=g->linftyNorm();
      conv=control().check(it,res);
//    cout << "\tresidual norm[" << it << "] = " << res << endl;
      if (conv!=SOLVER_ITERATE) break;
      precondition.vmult(*h,*g);
      double beta=gh;
      gh=(*g)*(*h);
      beta=gh/beta;
      d->scale(beta);
      d->add(-1.,*h);
    }
  }
  delete g;
  delete h;
  delete d;
  delete Ad;
  CHECK_TEST(conv==SOLVER_SUCCESS);
//cout << "\tnumber iterations = " << it << endl;
//x.printOn(cout);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void CGSolver::printOn(ostream &os) const {
  os << "CGSolver:" << endl;
  os << "\tres2 = " << res2 << endl;
  Solver::printOn(os);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GMRESSolver::givensRotation(Vector &h,Vector &b,Vector &ci,
Vector &si,int col) const {
  int n=1;
  int incx=1;
  for (int i=0;i<col;i++) {
    F77NAME(drot)(n,&h[i],incx,&h[i+1],incx,ci[i],si[i]);
  };
  double db=h[col+1];
  F77NAME(drotg)(&h[col],&db,ci[col],si[col]);
  b[col+1]=-si[col]*b[col];
  b[col]*=ci[col];
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
template<class T> SOLVER_STATE GMRESSolver::templateSolve(const T &A,
void (T::*vmult)(Vector&,const Vector&) const,
Vector &x,const Vector &b,const SolverPreconditioner &precondition) {
//TRACER_CALL(t,"GMRESSolver::templateSolve");
  RelativeErrorControl *rec=
    dynamic_cast<RelativeErrorControl*>(&control());
  if (rec!=0) rec->setRhsNorm(b.linftyNorm());

  NumPtr<Vector*> tmp_vectors(max_n_tmp_vectors);
  tmp_vectors.initialize(0);

  int accumulated_iterations=0;
  Matrix H(max_n_tmp_vectors,max_n_tmp_vectors-1);
  Vector gamma(max_n_tmp_vectors),ci(max_n_tmp_vectors-1),
    si(max_n_tmp_vectors-1);
  SOLVER_STATE iteration_state=SOLVER_ITERATE;
  Vector *&v=tmp_vectors[0];
  if (!v) v=x.clone();
  Vector *&p=tmp_vectors[max_n_tmp_vectors-1];
  if (!p) p=x.clone();
  Vector *r=0;
  Vector *x_=0;
  Vector *gamma_=0;
  if (!use_default_residual) {
    r=x.clone();
    x_=x.clone();
    gamma_=gamma.clone();
  }
  do {
    Vector h(max_n_tmp_vectors-1);
    if (!right_preconditioning) {
      (A.*vmult)(*p,x);
      p->scale(-1.);
      p->add(1.,b);
      precondition.vmult(*v,*p);
    } else {
      (A.*vmult)(*v,x);
      v->scale(-1.);
      v->add(1.,b);
    }
    double rho=v->l2Norm();
    if (use_default_residual) {
      iteration_state=
        control().check(accumulated_iterations,v->linftyNorm());
      if (iteration_state!=SOLVER_ITERATE) break;
    } else {
      if (!right_preconditioning) {
        (A.*vmult)(*r,x);
        r->scale(-1.);
        r->add(1.,b);
      } else precondition.vmult(*r,*v);
      double res=r->linftyNorm();
      iteration_state=control().check(accumulated_iterations,res);
      if (iteration_state!=SOLVER_ITERATE) {
        if (r) delete r; r=0;
        if (x_) delete x_; x_=0;
        if (gamma_) delete gamma_; gamma_=0;
        break;
      }
    }
    gamma[0]=rho;
    v->scale(1./rho);
    int dim=0;
    for (int inner_iteration=0;
    ((inner_iteration<max_n_tmp_vectors-2) &&
    (iteration_state==SOLVER_ITERATE));inner_iteration++) {
      ++accumulated_iterations;
      Vector *&vv=tmp_vectors[inner_iteration+1];
      if (!vv) vv=x.clone();
      if (!right_preconditioning) {
        (A.*vmult)(*p,*tmp_vectors[inner_iteration]);
        precondition.vmult(*vv,*p);
      } else {
        precondition.vmult(*p,*tmp_vectors[inner_iteration]);
        (A.*vmult)(*vv,*p);
      }
      dim=inner_iteration+1;
      for (int i=0;i<dim;i++) {
        h[i]=(*vv)*(*tmp_vectors[i]);
        vv->add(-h[i],*tmp_vectors[i]);
      }
      double s=vv->l2Norm();
      h[inner_iteration+1]=s;
      for (int i=0;i<dim;i++) {
        double htmp=(*vv)*(*tmp_vectors[i]);
        h[i]+=htmp;
        vv->add(-htmp,* tmp_vectors[i]);
      };
      s=vv->l2Norm();
      h[inner_iteration+1]=s;
      vv->scale(1./s);
      givensRotation(h,gamma,ci,si,inner_iteration);
      for (int i=0;i<dim;i++) {
        H(i,inner_iteration)=h[i];
      }
      rho=abs(gamma[dim]);
      if (use_default_residual) {
        iteration_state=control().check(accumulated_iterations,rho);
      } else {
        Vector h_(dim);
        *x_=x;
        *gamma_=gamma;
        Matrix H1(dim+1,dim);
        H.copyInto(H1);

        h_=*gamma_;
        int arg_m=h_.size();
        int arg_n=1;
        double arg_alpha=1.;
        int arg_lda=H1.size(0);
        int arg_ldb=h_.size();
        F77NAME(dtrsm)('L','U','N','N',&arg_m,&arg_n,&arg_alpha,
          H1.addr(),&arg_lda,h_.addr(),&arg_ldb);

        if (!right_preconditioning) {
          for (int i=0;i<dim;i++) {
            x_->add(h_[i],*tmp_vectors[i]);
          }
        } else {
          *p=0.;
          for (int i=0;i<dim;i++) {
            p->add(h_[i],*tmp_vectors[i]);
          }
          precondition.vmult(*r,*p);
          x_->add(1.,*r);
        }
        (A.*vmult)(*r,*x_);
        r->scale(-1.);
        r->add(1.,b);
        if (!right_preconditioning) {
          double res=r->linftyNorm();
          iteration_state=control().check(accumulated_iterations,res);
        } else {
          precondition.vmult(*x_, *r);
          double preconditioned_res=x_->linftyNorm();
          iteration_state=
            control().check(accumulated_iterations,preconditioned_res);
        }
      }
    }
    Vector hh(dim);
    Matrix H1(dim+1,dim);
    H.copyInto(H1);
    int inc=1;
    int arg_m=hh.size();
    F77NAME(dcopy)(arg_m,gamma.addr(),inc,hh.addr(),inc);
    int arg_n=1;
    double arg_alpha=1.;
    int arg_lda=H1.size(0);
    int arg_ldb=hh.size();
    F77NAME(dtrsm)('L','U','N','N',&arg_m,&arg_n,&arg_alpha,
      H1.addr(),&arg_lda,hh.addr(),&arg_ldb);

    if (!right_preconditioning) {
      for (int i=0;i<dim;i++) x.add(hh[i],*tmp_vectors[i]);
    } else {
      *p=0.;
      for (int i=0;i<dim;i++) p->add(hh[i],*tmp_vectors[i]);
      precondition.vmult(*v,*p);
      x.add(1.,*v);
    }
  } while (iteration_state==SOLVER_ITERATE);
  if (!use_default_residual) {
    if (r) delete r; r=0;
    if (x_) delete x_; x_=0;
    if (gamma_) delete gamma_; gamma_=0;
  }
  for (int i=0;i<tmp_vectors.getNumber();i++) {
    if (tmp_vectors[i]) delete tmp_vectors[i]; tmp_vectors[i]=0;
  }
//cout << "\tnumber iterations = " << accumulated_iterations << endl;
  return iteration_state;
}
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//adapted from DSPIGM in DASKR ode package
template<class T> SOLVER_STATE GMRESSolver::daskrSolve(const T &A,
void (T::*vmult)(Vector&,const Vector&) const,
Vector &x,const Vector &b,const SolverPreconditioner &precondition) {
  int neq=b.size();
  int max_steps=cntrl.maxSteps();

//maxl = max number of Arnoldi iterations before restart
//int maxl=min(5,neq); // default in DASKR
  int maxl=neq;
  maxl=min(maxl,max_steps);
  max_n_tmp_vectors=min(max_n_tmp_vectors,maxl);

  int nrmax=max_steps/maxl; // max restarts
  if (nrmax*maxl<max_steps) nrmax++;

  double bnrm=F77NAME(dnrm2)(neq,b.addr(),1);
  RelativeErrorControl *rec=
    dynamic_cast<RelativeErrorControl*>(&control());
  if (rec!=0) rec->setRhsNorm(bnrm);

  int maxlp1=maxl+1;
  int ldhes=maxl+2;
  Vector *dl=b.clone();
  Vector *r=b.clone();
  Vector *dx=x.clone();
  Vector *Av=b.clone();
  Vector q(2*maxl);
  Vector y(maxlp1);
  Matrix hes(ldhes,maxlp1);
  Matrix orthonormal_vectors(neq,maxlp1);

  (*r)=b;
  x=0.;
  SOLVER_STATE iteration_state=SOLVER_ITERATE;
  int its=0;
  int iflag=1;
  for (int nrsts=0;iflag==1 && nrsts<nrmax;nrsts++) {
    if (nrsts>0) (*r)=(*dl);
    iflag=0;
    int lgmr=0;
    (*dx)=0.;
    if (nrsts==0) {
      Vector Pb(orthonormal_vectors.addr(),neq);
      precondition.vmult(Pb,b);
    } else {
      memcpy(orthonormal_vectors.addr(),r->addr(),neq*sizeof(double));
    }
    double rnrm=F77NAME(dnrm2)(neq,orthonormal_vectors.addr(),1);
//  if (rnrm>eplin) 
    {
      double da=1./rnrm;
      F77NAME(dscal)(neq,da,orthonormal_vectors.addr(),1);
      hes=0.;
      double prod=1.;
      bool failed=false;
      int ll=0;
      double rho=numeric_limits<double>::infinity();
      double snormw=numeric_limits<double>::infinity();
      double *v_ll=orthonormal_vectors.addr();
      for (;ll<maxl;ll++,its++) {
        lgmr=ll;
        int llp1=ll+1;
        double *v_llp1=v_ll+neq;
        Vector v_in(v_ll,neq);
        (A.*vmult)(*Av,v_in);
        Vector PAv(v_llp1,neq);
        precondition.vmult(PAv,*Av); //v_llp1=preconditioner*jacobian*v_ll
//      Gram-Schmidt orthogonalization of v_llp1 
//        against cols max(0,ll-max_n_tmp_vectors+1),...,ll 
//        of orthonormal_vectors
//      snormw = || v_llp1 ||
        F77NAME(dorth)(v_llp1,orthonormal_vectors.addr(),hes.addr(),neq,
          llp1,ldhes,max_n_tmp_vectors,snormw);
        hes(llp1,ll)=snormw;
//      update QR factorization of hes
        int info=INT_MAX;
        F77NAME(dheqr)(hes.addr(),ldhes,llp1,q.addr(),info,llp1);
        CHECK_TEST(info==0);
        if (info==llp1) {
          failed=true;
          break;
        }
        prod*=q[2*ll+1];
        rho=abs(prod)*rnrm;
//      if (llp1>max_n_tmp_vectors then orthonormal_vectors 0...ll+1
//      are not necessarily orthogonal for ll>max_n_tmp_vectors;
//      the vector P_{ll+1} Q_{ll+1}^T e_{ll+1} must be computed, and
//      its norm used in the calculation of rho
        if (llp1>max_n_tmp_vectors && max_n_tmp_vectors<maxl) {
          if (ll==max_n_tmp_vectors) {
            memcpy(dl->addr(),orthonormal_vectors.addr(),
              neq*sizeof(double));
            for (int i=0;i<max_n_tmp_vectors;i++) {
              int i2=i*2;
              double s=q[i2+1];
              double c=q[i2];
              double *v_ip1=&orthonormal_vectors(0,i+1);
              for (int k=0;k<neq;k++) (*dl)[k]=s*(*dl)[k]+c*v_ip1[k];
            }
          } 
          double s=q[2*ll+1];
          double c=q[2*ll]/snormw;
          for (int k=0;k<neq;k++) (*dl)[k]=s*(*dl)[k]+c*v_llp1[k];
          // dl = P_{k+1} Q_k^T e_{k+1}
          rho*=dl->l2Norm();
        }
        iteration_state=cntrl.check(its,rho);
        if (iteration_state!=SOLVER_ITERATE) break;
        if (ll+1<maxl) F77NAME(dscal)(neq,1./snormw,v_llp1,1);
        v_ll=v_llp1;
      } // end loop over ll
      if (rho>=rnrm) failed=true;
      if (failed) { // statement 120
        iflag=2;
        (*dx)=0.;
      } else if (iteration_state==SOLVER_ITERATE) { // statement 150
        iflag=1;
        if (max_n_tmp_vectors==maxl) {
          memcpy(dl->addr(),orthonormal_vectors.addr(),neq*sizeof(double));
          for (int i=0;i<maxl-1;i++) {
            int i2=i*2;
            double s=q[i2+1];
            double c=q[i2];
            double *v_ip1=&orthonormal_vectors(0,i+1);
            for (int k=0;k<neq;k++) (*dl)[k]=s*(*dl)[k]+c*v_ip1[k];
          }
          double s=q[2*maxl-1];
          double c=q[2*maxl-2];
          double *v_maxl=&orthonormal_vectors(0,maxl);
          for (int k=0;k<neq;k++) (*dl)[k]=s*(*dl)[k]+c*v_maxl[k];
        }
        dl->scale(rnrm*prod);
        //-residual=dl = P_{k+1} Q_k^T e_{k+1} e_{k+1}^T Q_k e_0 || r_0 ||
      }
      if (!failed) { // statement 200
        ll=lgmr;
        int llp1=ll+1;
        for (int k=0;k<=llp1;k++) y[k]=0.;
        y[0]=rnrm; // e_0 * || r_0 ||
//      solve least squares problem for y
        F77NAME(dhels)(hes.addr(),ldhes,llp1,q.addr(),y.addr()); 
        (*dx)=0.;
        double *vi=orthonormal_vectors.addr();
        for (int i=0;i<=ll;i++,vi+=neq) {
          F77NAME(daxpy)(neq,y[i],vi,1,dx->addr(),1);
        } // dx = P y
        x+=(*dx);
      }
    }
  }
  delete dl; dl=0;
  delete r; r=0;
  delete dx; dx=0;
  delete Av;
  return iteration_state;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GMRESSolver::printOn(ostream &os) const {
  os << "GMRESSolver:" << endl;
  os << "\tmax_n_tmp_vectors = " << max_n_tmp_vectors 
     << "\n\tright_preconditioning = " << right_preconditioning
     << "\n\tuse_default_residual = " << use_default_residual
     << endl;
  Solver::printOn(os);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void MinResSolver::solve(const SolverMatrix &A,Vector &x,
const Vector &b,const SolverPreconditioner &precondition) {
//TRACER_CALL(t,"MinResSolver::solve");
  SOLVER_STATE conv=SOLVER_ITERATE;
  Vector* u[3]={b.clone(),b.clone(),b.clone()};
  Vector* m[3]={b.clone(),b.clone(),b.clone()};
  Vector* v=b.clone();
  double delta[3];
  double f[2];
  double e[2]; 
  double tau=0.;
  double c=0.;
  double s=0.;
  double d_=0.;
  int it=1;
  A.vmult(*m[0],x);
  *u[1]=b;
  *u[1]-=(*m[0]);
  precondition.vmult(*v,*u[1]);
  delta[1]=(*v)*(*u[1]);
  CHECK_TEST(delta[1]>=0);
  double r0=sqrt(delta[1]);
  double r_l2=r0;
  *u[0]=0.;
  delta[0]=1.;
  *m[0]=0.;
  *m[1]=0.;
  *m[2]=0.;
  conv=control().check(0,r_l2);
  while (conv==SOLVER_ITERATE) {      
    if (delta[1]>0.) v->scale(1./sqrt(delta[1]));
    else *v=0.;
    A.vmult(*u[2],*v);
    u[2]->add(-sqrt(delta[1]/delta[0]),*u[0]);
    double gamma=(*u[2])*(*v);
    u[2]->add(-gamma/sqrt(delta[1]),*u[1]);
    *m[0]=(*v);
    precondition.vmult(*v,*u[2]);
    delta[2]=(*v)*(*u[2]);
    CHECK_TEST(delta[2]>=0.);
    if (it==1) {
      d_=gamma;
      e[1]=sqrt(delta[2]);
    }
    if (it>1) {
      d_=s*e[0]-c*gamma;
      e[0]=c*e[0]+s*gamma;
      f[1]=s*sqrt(delta[2]);
      e[1]=-c*sqrt(delta[2]);
    }
    double d=sqrt(d_*d_+delta[2]);
    if (it>1) tau*=s/c;
    c=d_/d;
    tau*=c;
    s=sqrt(delta[2])/d;
    if (it==1) tau=r0*c;
    m[0]->add(-e[0],*m[1]);
    if (it>1) m[0]->add(-f[0],*m[2]);
    m[0]->scale(1./d);
    x.add(tau,*m[0]);
    r_l2*=abs(s);
    conv=control().check(it,r_l2);
    it++;
    m[2]->swapWith(*m[1]);
    m[1]->swapWith(*m[0]);
    u[0]->swapWith(*u[1]);
    u[1]->swapWith(*u[2]);
    f[0]=f[1];
    e[0]=e[1];
    delta[0]=delta[1];
    delta[1]=delta[2];
  }
  if (u[0]) delete u[0]; u[0]=0;
  if (u[1]) delete u[1]; u[1]=0;
  if (u[2]) delete u[2]; u[2]=0;
  if (m[0]) delete m[0]; m[0]=0;
  if (m[1]) delete m[1]; m[1]=0;
  if (m[2]) delete m[2]; m[2]=0;
  if (v) delete v; v=0;
//cout << "\tnumber iterations = " << it << endl;
  CHECK_TEST(conv==SOLVER_SUCCESS);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void MinResSolver::printOn(ostream &os) const {
  os << "MinResSolver:" << endl;
  os << "\tres2 = " << res2 << endl;
  Solver::printOn(os);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void QMRSSolver::solve(const SolverMatrix &A,Vector &x,
const Vector &b,const SolverPreconditioner &precondition) {
  Vector *Vv=x.clone();
  Vector *Vp=x.clone();
  Vector *Vq=x.clone();
  Vector *Vt=x.clone();
  Vector *Vd=x.clone();
  int step=0;
  bool keep_going=true;
  SOLVER_STATE state=SOLVER_ITERATE;
  do {
    state=SOLVER_ITERATE;
    int it=0;
    double theta=0.;
    *Vd=0.;
    precondition.vmult(*Vq,x);  
    A.vmult(*Vv,*Vq);
    Vv->scale(-1.);
    Vv->add(1.,b);
    double res=Vv->l2Norm();
    if (control().check(step,Vv->linftyNorm())==SOLVER_SUCCESS) {
      keep_going=false;
      break;
    }
    *Vp=*Vv;
    precondition.vmult(*Vq,*Vp);
    double tau=(*Vv)*(*Vv);
    double rho=(*Vq)*(*Vv);
    keep_going=false;
    while (state==SOLVER_ITERATE) {
      step++;
      it++;
      A.vmult(*Vt,*Vq);
      double sigma=(*Vq)*(*Vt);
      if (abs(sigma/rho)<breakdown) {
        keep_going=true;
        break;
      }
      double alpha=rho/sigma;
      Vv->add(-alpha,*Vt);
      double theta_old=theta;
      theta=(*Vv)*(*Vv)/tau;
      double psi=1./(1.+theta);
      tau*=theta*psi;
      Vd->scale(psi*theta_old);
      Vd->add(psi*alpha,*Vp);
      x.add(1.,*Vd);
      if (exact_residual) {
        A.vmult(*Vq,x);
        Vq->scale(-1.);
        Vq->add(1.,b);
        res=Vq->l2Norm();
      } else res=sqrt(static_cast<double>(it+1)*tau);
      state=control().check(step,res);
      if ((state==SOLVER_SUCCESS) || (state==SOLVER_FAILURE)) {
        A.vmult(*Vq,x);
        keep_going=false;
        break;
      }
      if (abs(rho)<breakdown) {
        keep_going=true;
        break;
      }
      double rho_old=rho;
      precondition.vmult(*Vq,*Vv);
      rho=(*Vq)*(*Vv);
      double beta=rho/rho_old;
      Vp->scale(beta);
      Vp->add(1.,*Vv);
      precondition.vmult(*Vq,*Vp);
    }
  } while (keep_going);
  if (Vv) delete Vv; Vv=0;
  if (Vp) delete Vp; Vp=0;
  if (Vq) delete Vq; Vq=0;
  if (Vt) delete Vt; Vt=0;
  if (Vd) delete Vd; Vd=0;
  CHECK_TEST(state==SOLVER_SUCCESS);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void QMRSSolver::printOn(ostream &os) const {
  os << "QMRSSolver:" << endl;
  os << "\tres2 = " << res2 
     << "\n\texact_residual = " << exact_residual
     << "\n\tbreakdown = " << breakdown
     << endl;
  Solver::printOn(os);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void RichardsonSolver::solve(const SolverMatrix &A,Vector &x,
const Vector &b,const SolverPreconditioner &precondition) {
  SOLVER_STATE conv=SOLVER_ITERATE;
  Vector *r=x.clone();
  Vector *d=x.clone();
  for (int iter=0;conv==SOLVER_ITERATE;iter++) {
    A.vmult(*r,x);
    r->scale(-1.);
    r->add(1.,b);
    if (!use_preconditioned_residual) {
      conv=control().check(iter,r->linftyNorm());
      if (conv!=SOLVER_ITERATE) break;
    }
    precondition.vmult(*d,*r);
    if (use_preconditioned_residual) {
      conv=control().check(iter,d->linftyNorm());
      if (conv!=SOLVER_ITERATE) break;
    }
    x.add(omega,*d);
  }
  delete r;
  delete d;
  CHECK_TEST(conv==SOLVER_SUCCESS);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void RichardsonSolver::transposeSolve(const SolverMatrix &A,
Vector &x,const Vector &b,
const SolverPreconditioner &precondition) {
  SOLVER_STATE conv=SOLVER_ITERATE;
  Vector *r=x.clone();
  Vector *d=x.clone();
  for(int iter=0;conv==SOLVER_ITERATE;iter++) {
    A.transposeVmult(*r,x);
    r->scale(-1.);
    r->add(1.,b);
    conv=control().check(iter,r->linftyNorm());
    if (conv!=SOLVER_ITERATE) break;
    precondition.transposeVmult(*d,*r);
    x.add(omega,*d);
  }
  delete r;
  delete d;
  CHECK_TEST(conv==SOLVER_SUCCESS);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void RichardsonSolver::printOn(ostream &os) const {
  os << "RichardsonSolver:" << endl;
  os << "\tomega = " << omega 
     << "\n\tuse_preconditioned_residual = " 
     << use_preconditioned_residual << endl;
  Solver::printOn(os);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void HdivCGSolver::solve(const SolverMatrix &A,const SolverMatrix &B,
Vector &x,const Vector &b,const SolverPreconditioner &precondition) {
  CHECK_SAME(A.size(0),A.size(1));
  CHECK_SAME(B.size(1),A.size(1));
  CHECK_SAME(x.size(),A.size(0)+B.size(0));
  CHECK_SAME(x.size(),b.size());

  RelativeErrorControl cg_cntrl(A.size(0),cntrl.tolerance(),
    cntrl.logHistory(),cntrl.logResult());
  CGSolver cg(cg_cntrl);
  RelativeErrorControl *rec=
    dynamic_cast<RelativeErrorControl*>(&control());
  if (rec!=0) rec->setRhsNorm(b.linftyNorm());

  Vector x_A(&x[0],A.size(0));
  Vector x_B(x.addr()+A.size(0),B.size(0));
  const Vector b_B(const_cast<double*>(b.addr())+A.size(0),B.size(0));
  Vector rhs_A(A.size(0));
  Vector r_B(B.size(0));
  Vector p_B(B.size(0));
  Vector z_B(B.size(0));
  Vector p_A(A.size(0));
  Vector v_A(A.size(0));

  B.transposeVmult(rhs_A,x_B);
  rhs_A.scale(-1.);
  rhs_A.add(1.,b); // b_A - B^T x_B
  cg.solve(A,x_A,rhs_A,precondition); // x_A = A^{-1} ( b_A - B^T x_B )

  B.vmult(r_B,x_A);
  r_B.add(-1.,b_B);
  p_B=r_B;
  double qq=p_B*r_B;

  double res=p_B.linftyNorm();
  SOLVER_STATE conv=control().check(0,res);
  int it=0;
  while (conv==SOLVER_ITERATE) {
    it++;
    B.transposeVmult(v_A,p_B);
    p_A=0.;
    cg.solve(A,p_A,v_A,precondition);
    B.vmult(z_B,p_A);

    double alpha=p_B*z_B;
    alpha=qq/alpha;
    x_B.add(alpha,p_B);
    x_A.add(-alpha,p_A);
    r_B.add(-alpha,z_B);
    res=r_B.linftyNorm();
    conv=control().check(it,res);
    if (conv!=SOLVER_ITERATE) break;
    double beta=qq;
    qq=r_B*r_B;
    beta=qq/beta;
    p_B.scale(beta);
    p_B.add(1.,r_B);
  }
  CHECK_TEST(conv==SOLVER_SUCCESS);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void HdivCGSolver::printOn(ostream &os) const {
  os << "HdivCGSolver:" << endl;
  Solver::printOn(os);
}

#include "NumPtr.C"
INSTANTIATE_NUMPTR(Vector*);
template SOLVER_STATE GMRESSolver::daskrSolve<SolverMatrix>(
  SolverMatrix const&, 
  void (SolverMatrix::*)(Vector&, Vector const&) const, 
  Vector&, 
  Vector const&, 
  SolverPreconditioner const&
);
template SOLVER_STATE GMRESSolver::templateSolve<SolverMatrix>(
  SolverMatrix const&, 
  void (SolverMatrix::*)(Vector&, Vector const&) const, 
  Vector&, 
  Vector const&, 
  SolverPreconditioner const&
);
