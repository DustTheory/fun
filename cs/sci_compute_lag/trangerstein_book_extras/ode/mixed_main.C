#include <fstream.h>
#include <iostream.h>
//#include <sys/resource.h>
#include <stdlib.h>
#include <sys/times.h>
#include <sys/param.h>
#define MAX_STEPS 1000000

extern "C" {
  void integrate_(const int nsteps,const float &tmax, 
                  float *soln,float *time);
}
struct odepar_common {
  float rate;
};
extern odepar_common odepar_;

int main(int /*argc*/,char* /*argv[]*/) {
  int nsteps;
  float tmax;
  float t_array[MAX_STEPS+1],x_array[MAX_STEPS+1];

  cout << "enter rate" << endl;
  cin >> odepar_.rate;
  cout << "enter initial x" << endl;
  cin >> x_array[0];
  t_array[0]=0.;
  cout << "enter tmax" << endl;
  cin >> tmax;
  cout << "enter number steps" << endl;
  cin >> nsteps;
  if (nsteps<0) nsteps=0;
  if (nsteps>MAX_STEPS) nsteps=MAX_STEPS;
//cout << "number steps = " << nsteps << endl;
//cout << "rate,x = " << odepar_.rate << " " << x_array[0] << endl;
//cout << "tmax,nsteps = " << tmax << " " << nsteps << endl;

  struct tms buffer;
  times(&buffer);
  clock_t times_start=buffer.tms_utime+buffer.tms_stime;
  
//solve ode
  integrate_(nsteps,tmax, x_array,t_array);

//stop timing
  times(&buffer);
  clock_t times_stop=buffer.tms_utime+buffer.tms_stime;

  cout << "\ntimes: integrate took " 
       << double(times_stop-times_start)/double(HZ) << " seconds" 
       << endl;

  ofstream out_file;
  out_file.open("mixed_main_output",ios::out | ios::noreplace);
  for (int i=0;i<=nsteps;i++) {
    out_file << t_array[i] << " " << x_array[i] << endl;
  }

  return EXIT_SUCCESS;
}

#ifdef __GNUC__
extern "C" {
void MAIN__(void) {;}
}
#endif
