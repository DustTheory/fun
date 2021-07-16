#include <float.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include "Arch.H"
#include "Debug.H"
#include "SetTraps.H"
#include "InputParameter.H"
#include "Palette.H"
#include "TimedObject.H"
#include "Tracer.H"
#include "VGT.H"
#include "XColormap.H"
#include "XYGraphTool.H"
#include <unistd.h>

#define ifirst 0

#define left 0.
#define right 1.

#define ncells 64*16
#define tmax 0.1
#define ilast ncells
#define tol 1.e-8

int multigrid(double *r, double *d, int n) {
  int iter=1;
  double dx=1./double(n); 
  double dt=1./double(ncells);
  double a1=1./(6.*dt)-1./(2.*dx*dx);
  double a2=2.*(1./(3.*dt))+1./(dx*dx);
  double a3=a1;
  for(int i=0; i<n; i++) d[i]=0.;
  
  double *b = new double [n];
  if (n%4==0) {
    int cn=n/2;
    for(int k=0; k<iter; k++) {
      for (int i=1; i<n; i++) {
        b[i]=r[i];
	d[i]=(-b[i]-a1*d[i-1]-a3*d[i+1])/a2;
      } 
      for (int i=n-1; i>0; i--) {
	d[i]=(-b[i]-a1*d[i-1]-a3*d[i+1])/a2; 
      }	 
    }
    for (int i=1; i<n; i++) { 
      r[i]=b[i]+(2.*d[i]/dt+(-d[i-1]+2.*d[i]-d[i+1])/dx/dx);
    }
     
    double *cr = new double [cn]; 
    double *cd = new double [cn]; 		 
    for (int j=2; j<n; j+=2){ 
      cr[j/2]=0.5*r[j]+0.25*(r[j-1]+r[j+1]);
    }
    
    multigrid(cr, cd, cn); 
	
    for (int j=1; j<cn; j++) {
      d[2*j-1]+=0.5*cd[j];
      d[2*j+1]+=0.5*cd[j];
    }
    for (int j=1; j<cn; j++) d[2*j]+=cd[j]; 	

    delete [] cr;
    delete [] cd; 
    for (int i=1; i<n; i++) {	
      double dd=-d[i-1]+2.*d[i]-d[i+1];	       	
      r[i]=2.*(d[i]/dt+dd/dx/dx)+b[i]; 		
    }
    for(int k=0; k<iter; k++) {
      for (int i=1;i<n; i++) d[i]=(-b[i]-a1*d[i-1]-a3*d[i+1])/a2; 
      for (int i=n-1;i>0; i--) d[i]=(-b[i]-a1*d[i-1]-a3*d[i+1])/a2;
    }
  } else { 
    for (int i=1;i<n; i++) d[i]=-r[i]/a2;
  }
  delete [] b;
  return EXIT_SUCCESS; 
}

int main(int argc,char* argv[]) {
  ASSERT(argc>1);

//since ncells is determined dynamically, 
//must use dynamic memory allocation
  double *u=new double[ncells+1];    // -2 <= cell <= ncells+1
  double *u1=new double[ncells+1];
  double *du=new double[ncells+1];
  double *x=new double[ncells+1];    //  0 <= edge <= ncells
  double *a=new double[ncells+1];
  double *b=new double[ncells+1];
  double *c=new double[ncells+1];
  double *f=new double[ncells+1];
  double *r=new double[ncells+1];
  double *a1=new double[ncells+1];
  double *b1=new double[ncells+1];
  double *c1=new double[ncells+1];

// initialization
  double dx = 1./double(ncells);
  for (int i=0;i<=ncells;i++) {
    x[i]=left+(double)(i)*dx;
  }
  for (int i=0;i<=ncells;i++) {
//  u1[i]=u[i]=sin(2.*M_PI*x[i]);
//  u1[i]=u[i]=0.5-fabs(0.5-x[i]);
    if (i<=int(ncells/3.) || i> int(2.*ncells/3.)) u1[i]=u[i]=0.;
    else  u1[i]=u[i]=1.;
  }
//find min,max data values
  double xlo=x[ifirst];
  double xhi=x[ifirst];
  for (int ie=ifirst+1;ie<=ilast;ie++) {
    double xi=x[ie];
    if (xi<xlo) xlo=xi;
    if (xi>xhi) xhi=xi;
  }
  double ulo=u[ifirst];
  double uhi=u[ifirst];
    
  for (int ic=ifirst+1;ic<=ilast;ic++) {
    double ui=u[ic];
    if (ui<ulo) ulo=ui;
    if (ui>uhi) uhi=ui;
  }
  double first_uhi=uhi;
    
//setup interactive graphics
  Palette pal;
  XColormap cmap(&pal);
  double winsize=0.5;
  XYGraphTool gt("u vs x","x","u",xlo,xhi,ulo,uhi,&cmap,NULL,winsize);

  for (int n=0;n<10;n++) {
//  initialize graphics display
    gt.setbgColor("white");
    gt.setfgColor("blue");
    gt.drawAxes();

//  draw numerical solution
    int ic=ifirst;
    double width=0.5*(x[ic+1]-x[ic]);
//  double center=0.5*(x[ic]+x[ic+1]);
    double center=x[ic];
    gt.setfgColor("blue");
    gt.drawPlus(center,u[ic],width);

    for (ic=ifirst+1;ic<ilast;ic++) {
      double width=0.5*(x[ic+1]-x[ic]);
//    double center=0.5*(x[ic]+x[ic+1]);
      double center=x[ic];
      gt.drawPlus(center,u[ic],width);
    }

//  force X to perform the requests
    gt.flush();
    XYGraphTool::WINDOW_TYPE::QuitButton qb;

    double dt=dx;
    for (int ie=1; ie<ncells; ie++) {	    
      a[ie]=1./(6.*dt)-1./(2.*dx*dx); 
      b[ie]=2.*(1./(3.*dt))+1./(dx*dx);
      c[ie]=a[ie];	
      a1[ie]=1./(6.*dt)+1./(2.*dx*dx);
      b1[ie]=2.*(1./(3.*dt))-1./(dx*dx);
      c1[ie]=a1[ie];
      f[ie]=u[ie-1]*a1[ie]+u[ie]*b1[ie]+u[ie+1]*c1[ie];
    }

/*
    for (int ie=1; ie<ncells; ie++) {		
      u1[ie]=(f[ie]-a[ie]*u1[ie-1]+c[ie]*u1[ie+1])/b[ie];
    }

    for (int ie=ncells-1; ie>0; ie--) {
      u1[ie]=(f[ie]-a[ie]*u1[ie-1]+c[ie]*u1[ie+1])/b[ie];
    }
*/
    int it=0;
    double error=0.;      	
    for (int ie=1; ie<ncells; ie++) {	
      r[ie]=(a[ie]*u1[ie-1]+b[ie]*u1[ie]+c[ie]*u1[ie+1])-f[ie];    
      if (fabs(r[ie])>error) error=fabs(r[ie]);
    }
    r[0]=r[ncells]=0.;
//  cout<<"it="<<it<<"	error="<<error<<"\n";
    cout<<n<<"\t"<<log(error)/log(10.)<<"\n";

    it++;
    multigrid(r, du, ncells);
    for (int ie=1; ie<ncells; ie++)  u1[ie]+=du[ie]; 
  }

  {
 // find min,max data values
    double ulo=u1[0];
    double uhi=u1[0];
    
    for (int ic=ifirst;ic<=ilast;ic++) {
      double ui=u1[ic];
      if (ui<ulo) ulo=ui;
      if (ui>uhi) uhi=ui;
    }

//  initialize graphics display

//  gt.rescale(xlo,xhi,ulo,uhi);
    gt.rescale(0.,1.,ulo,uhi);

    gt.newPage();
    gt.setbgColor("white");
    gt.setfgColor("black");
    gt.drawAxes();

//  draw numerical solution
    int ic=ifirst;

    double width=0.5*(x[ic+1]-x[ic]);
    double center=0.5*(x[ic]+x[ic+1]);
//  double center=x[ic];
    gt.setfgColor("blue");
    gt.drawPlus(center,u1[ic],width);

    for (ic=ifirst+1;ic<ilast;ic++) {
      double width=0.5*(x[ic+1]-x[ic]);
      double center=0.5*(x[ic]+x[ic+1]);
//    double center=x[ic];
      gt.drawPlus(center,u1[ic],width);
    }
            
//  force X to perform the requests
    gt.flush();
  }
  XYGraphTool::WINDOW_TYPE::QuitButton qb;

//since x, u, flux were created with operator new, must delete
  delete [] u;
  delete [] du;
  delete [] x;
  delete [] u1;
  delete [] b;
  delete [] a;
  delete [] c;
  delete [] f;delete [] r;
  delete [] b1;
  delete [] a1;
  delete [] c1;
  
  return EXIT_SUCCESS;
}
