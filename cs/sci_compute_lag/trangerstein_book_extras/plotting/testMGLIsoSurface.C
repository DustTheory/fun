#include <iostream>
#include <mgl2/mgl.h>
#include <MemoryDebugger.H>
using namespace std;
int main(int argc,char* argv[]) {
  cout << boolalpha;
  MemoryDebugger md(1);
  {
    mglGraph gr;
    gr.Light(true);
    gr.Alpha(true);
    gr.SetOrigin(0.,0.,0.);
    gr.Title("IsoSurface x*x+y*y+z*z");
    gr.SetRanges(-1.,1.,-1.,1.,-1.,1.);

    mglData a(50,50,50);
    a.Fill(gr.Self(),"x*x+y*y+z*z");
    gr.Rotate(40,60);
    gr.Surf3(0.5,a);
    gr.Surf3(1.,a);
    gr.Box();
    gr.Axis();
    gr.Grid();
    gr.Label('x',"xaxis");
    gr.Label('y',"yaxis");
    gr.Label('z',"zaxis");
    gr.ShowImage("xv");
  }
  return 0;
}
