//adapted from DSPIGM in DASKR ode package
template<class T> void GMRESSolver::daskrSolve(T &A,
void (T::*vmult)(Vector&,const Vector&) const,
Vector &x,const Vector &b) {
//int n_orthogonalization_vectors,int maxl,
//TRACER_CALL(tr,"dslvk");
#ifdef DEBUG
  for (int i=0;i<x.size();i++) {
    cout << "\tx[" << i << "] = " << x[i] << endl;
  }   
  for (int i=0;i<x.size();i++) {
    cout << "\tb[" << i << "] = " << b[i] << endl;
  }
#endif
  int neq=b.size();
  int max_steps=control.maxSteps();

  int maxl=min(5,neq); // default in DASKR
  maxl=min(maxl,max_steps);
  max_n_tmp_vectors=min(max_n_tmp_vectors,maxl);

  int nrmax=max_steps/maxl; // max restarts
  if (nrmax*maxl<max_steps) nrmax++;

  RelativeErrorControl *rec=
    dynamic_cast<RelativeErrorControl*>(&control());
  if (rec!=0) rec->setRhsNorm(b.linftyNorm());

  int maxlp1=maxl+1;
  Vector *dl=b.clone();
  Vector *r=b.clone();
  Vector *dx=x.clone();
  Vector *Av=b.clone();
  double *q=OPERATOR_NEW_BRACKET(double,2*neq);
  double *hes=OPERATOR_NEW_BRACKET(double,maxlp1*maxl);
  double *krylov_vectors=OPERATOR_NEW_BRACKET(double,neq*maxlp1);
  double *work=OPERATOR_NEW_BRACKET(double,maxlp1);

  (*r)=b;
  x=0.;
  int iflag=1;
  for (int nrsts=0;iflag==1 && nrsts<nrmax;nrsts++) {
    if (nrsts>0) (*r)=(*dl);
    iflag=0;
    int lgmr=0;
    x=0.;
    if (nrsts==0) precondition.vmult(krylov_vectors,b);
    else memcpy(krylov_vectors,b.addr(),neq*sizeof(double));
    double rnrm=F77NAME(dnrm2)(neq,krylov_vectors,1);
//  if (rnrm>eplin) 
    {
      double da=1./rnrm;
      F77NAME(dscal)(neq,da,krylov_vectors,1);
      for (int i=0;i<maxlp1*maxl;i++) hes[i]=0.;
      double prod=1.;
      bool failed=false;
      int ll=0;
      double rho=numeric_limits<double>::infinity();
      double snormw=numeric_limits<double>::infinity();
      for (;ll<maxl;ll++) {
        lgmr=ll;
        int llp1=ll+1;
        double *v_ll=krylov_vectors+ll*neq;
        double *v_llp1=krylov_vectors+llp1*neq;
//      v_llp1 = preconditioner * jacobian * v_ll
        Vector v_in(v_ll,neq);
        A.vmult(*Av,v_in);
        Vector PAv(v_llp1,neq);
        precondition.vmult(PAv,*Av);
//      Gram-Schmidt orthogonalization of v_llp1 
//        against cols 0,...,ll of krylov_vectors
//      snormw = || v_llp1 ||
        F77NAME(dorth)(v_llp1,krylov_vectors,hes,neq,llp1,maxlp1,
          max_n_tmp_vectors,snormw);
        hes[llp1+ll*maxlp1]=snormw;
//      update QR factorization
        int info=INT_MAX;
        F77NAME(dheqr)(hes,maxlp1,llp1,q,info,llp1);
        if (info==llp1) {
          failed=true;
          break;
        }
        prod*=q[2*ll+1];
        rho=abs(prod*rnrm);
        if (llp1>max_n_tmp_vectors && max_n_tmp_vectors<maxl) {
          TRACER_CALL(tr,"dspigm");
          if (ll==max_n_tmp_vectors) {
            memcpy(dl.addr(),krylov_vectors,neq*sizeof(double));
            for (int i=0;i<max_n_tmp_vectors;i++) {
              int ip1=i+1;
              int i2=i*2;
              double s=q[i2+1];
              double c=q[i2];
              double *v_ip1=krylov_vectors+ip1*neq;
              for (int k=0;k<neq;k++) dl[k]=s*dl[k]+c*v_ip1[k];
            }
          }
          double s=q[2*ll+1];
          double c=q[2*ll]/snormw;
          for (int k=0;k<neq;k++) dl[k]=s*dl[k]+c*v_llp1[k];
          rho*=dl.l2Norm();
        }
        if (rho<=eplin*rnrm || llp1==maxl) break;
        F77NAME(dscal)(neq,1./snormw,v_llp1,1);
      }
      if (rho>=rnrm) failed=true;
      if (failed) {
        TRACER_CALL(tr,"dspigm failed");
        iflag=2;
        x=0.;
      } else if (rho>eplin*rnrm) {
        iflag=1;
        if (max_n_tmp_vectors==maxl) {
          memcpy(dl.addr(),krylov_vectors,neq*sizeof(double));
          for (int i=0;i<maxl-1;i++) {
            int ip1=i+1;
            int i2=i*2;
            double s=q[i2+1];
            double c=q[i2];
            double *v_ip1=krylov_vectors+ip1*neq;
            for (int k=0;k<neq;k++) {
              dl[k]=s*dl[k]+c*v_ip1[k];
            }
          }
          double s=q[2*maxl-1];
          double c=q[2*maxl-2]/snormw;
          double *v_maxl=krylov_vectors+maxl*neq;
          for (int k=0;k<neq;k++) {
            dl[k]=s*dl[k]+c*v_maxl[k];
          }
        }
        dl.scale(rnrm*prod);
      }
      if (!failed) {
        ll=lgmr;
        int llp1=ll+1;
        for (int k=0;k<=llp1;k++) work[k]=0.;
        work[0]=rnrm;
        F77NAME(dhels)(hes,maxlp1,llp1,q,work);
        x=0.;
        double *vi=krylov_vectors;
        for (int i=0;i<=ll;i++,vi+=neq) {
          F77NAME(daxpy)(neq,work[i],vi,1,x,1);
        }
      }
    }
    x+=(*dx);
  }
  delete dl; dl=0;
  delete r; r=0;
  delete dx; dx=0;
  delete Av;
  delete [] work;
  delete [] q;
  delete [] hes;
  delete [] krylov_vectors;
}
