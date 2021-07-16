#include <iostream>
#include <mgl2/mgl.h>
#include <MemoryDebugger.H>
using namespace std;
int main(int argc,char* argv[]) {
  cout << boolalpha;
  MemoryDebugger md(1);
  {
    mglGraph gr;
    gr.Title("contour x*x+y*y in [0,1]x[0,1]");
    gr.SetRanges(0.,1.,0.,1.);
    gr.Axis();
    gr.Label('x',"xaxis");
    gr.Label('y',"yaxis");
    gr.Box();
//  gr.SetSize(,); // in pixels; should be fraction of screen width
//  gr.SetPalette(char *colors);
//  gr.SetAxisStl(...);

    mglData a(50,50);
    a.Fill(gr.Self(),"x*x+y*y");
    mglData v(21);
    for (int i=0;i<21;i++) v.a[i]=static_cast<float>(i)*0.1;

    const char *color_scheme="BbcwyrR"; // hot-cold
    gr.Cont(v,a,"t");
    gr.Finish();
    gr.ShowImage("xv");
  }
  return 0;
}
