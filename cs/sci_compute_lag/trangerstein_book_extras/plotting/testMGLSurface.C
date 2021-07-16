#include <iostream>
#include <mgl2/mgl.h>
#include <MemoryDebugger.H>
using namespace std;
int main(int argc,char* argv[]) {
  cout << boolalpha;
  MemoryDebugger md(1);
  {
    mglGraph gr;
    gr.SetOrigin(0.,0.,0.);
    gr.Title("surface x*x+y*y in [0,1]x[0,1]");
    gr.SetRanges(0.,1.,0.,1.,0.,2.);

    mglData a(50,50);
    a.Fill(gr.Self(),"x*x+y*y");
    gr.Rotate(40,60);
    gr.Surf(a,"BbcwyrR#"); // hot-cold           
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
