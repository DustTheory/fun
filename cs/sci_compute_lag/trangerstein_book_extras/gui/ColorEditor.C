//**********************************************************************
// Copyright 2006 John A. Trangenstein
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
// "$Header: /home/faculty/johnt/cvs/deal_new/gui/ColorEditor.C,v 1.1 2009/08/20 17:32:34 johnt Exp $"
#include "ColorEditor.H"
#ifndef USE_GTK
extern "C" {
#include <Xm/DrawnB.h>
#include <Xm/DrawingA.h>
#include <Xm/Form.h>
#include <Xm/Frame.h>
#include <Xm/Label.h>
#include <Xm/MessageB.h>
#include <Xm/RowColumn.h>
#include <Xm/XmStrDefs.h>
#include <Xm/Scale.h>
#include <Xm/ToggleB.h>
}
#include "MemoryDebugger.H"
#include "Palette.H"
#include "Tracer.H"
#include "XWindow.H"
#include "ISLList.C"
Widget XColorEditor::dialog=0;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//XColormap::XColormap uses the DefaultColormap to avoid "flashing"
//  and checks that the Colormap is not StaticGray or GrayScale
//  sometimes, this Colormap is read only (StaticColor or TrueColor)
//  if so, then XColorEditor::xxxSliderMoved will cause an error
//    in calling XStoreColor
//Thus, we install a new visual if necessary
XColorEditor::XColorEditor(XColormap *xmap,Widget p) :
xcolormap(xmap),parent(p),current_toggle(0),form(0),selected(-1) {
  display=XtDisplay(parent);
  screen=DefaultScreen(display);

  XWindowAttributes a;
  Status status=XGetWindowAttributes(display,XtWindow(parent),&a);
  Visual *visual=a.visual;
  if ((a.visual->c_class!=DirectColor) && 
  (a.visual->c_class!=PseudoColor)) {
    XVisualInfo vi;
    int depth=a.depth;
    if (!XMatchVisualInfo(display,screen,depth,PseudoColor,&vi)) {
      ASSERT(XMatchVisualInfo(display,screen,depth,DirectColor,&vi));
    }
    xcolormap->installNewVisual(XtWindow(parent),vi.visual);
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
bool XColorEditor::warnUserNoColor(Widget dialog_parent) {
  if (current_toggle || selected!=-1) return FALSE;
  if (dialog==0)  {
    dialog=XmCreateWarningDialog(dialog_parent,
      const_cast<char*>("noColorWarningDialog"),0,0);
  }
  XtManageChild(dialog);
  return TRUE;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void XColorEditor::redSliderMovedCallback(Widget w,
XtPointer client_data,XtPointer callData) {
  XColorEditor *ce=reinterpret_cast<XColorEditor*>(client_data);
  XmScaleCallbackStruct *cb = 
    reinterpret_cast<XmScaleCallbackStruct*>(callData);
  if (ce->warnUserNoColor(w)) return;

  Pixel pixel;
  XtVaGetValues(ce->swatch,XmNbackground,&pixel,0);

  XColor color;
  color.red  =cb->value*256;
  color.pixel=pixel;
  color.flags=DoRed;
  XStoreColor(XtDisplay(w),ce->xcolormap->getColormap(),&color);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void XColorEditor::blueSliderMovedCallback(Widget w,
XtPointer client_data,XtPointer callData ) {
  XColorEditor *ce=reinterpret_cast<XColorEditor*>(client_data);
  XmScaleCallbackStruct *cb = 
  reinterpret_cast<XmScaleCallbackStruct*>(callData);
  if (ce->warnUserNoColor(w)) return;

  Pixel pixel;
  XtVaGetValues(ce->swatch,XmNbackground,&pixel,0);    

  XColor color;
  color.blue =cb->value*256;
  color.pixel=pixel;
  color.flags=DoBlue;

  XStoreColor(XtDisplay(w),ce->xcolormap->getColormap(),&color );
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void XColorEditor::greenSliderMovedCallback(Widget w,
XtPointer client_data,XtPointer callData ) {
  XColorEditor *ce=reinterpret_cast<XColorEditor*>(client_data);
  XmScaleCallbackStruct *cb = 
    reinterpret_cast<XmScaleCallbackStruct*>(callData);
  if (ce->warnUserNoColor(w)) return;

  Pixel pixel;
  XtVaGetValues(ce->swatch,XmNbackground,&pixel,0);    

  XColor color;
  color.green=cb->value*256;
  color.pixel=pixel;
  color.flags=DoGreen;

  XStoreColor(XtDisplay(w),ce->xcolormap->getColormap(),&color);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void XColorEditor::selectColorCallback(Widget w,XtPointer client_data,
XtPointer callData) {
  XColorEditor *ce=reinterpret_cast<XColorEditor*>(client_data);
  XmToggleButtonCallbackStruct *cbs =
    reinterpret_cast<XmToggleButtonCallbackStruct*>(callData);
  if (!cbs->set) return;

  ce->current_toggle=w;

  Pixel bg;
  XtVaGetValues(ce->current_toggle,XmNbackground,&bg,0);
  XtVaSetValues(ce->swatch,XmNbackground,bg,0);

  XColor color;
  color.flags=DoRed | DoGreen | DoBlue;
  color.pixel=bg;

  XQueryColor(XtDisplay(w),ce->xcolormap->getColormap(),&color);

  XtVaSetValues(ce->red_slider,  XmNvalue,color.red  /256,0);
  XtVaSetValues(ce->green_slider,XmNvalue,color.green/256,0);
  XtVaSetValues(ce->blue_slider, XmNvalue,color.blue /256,0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void XColorEditor::resizeCallback(Widget w,XtPointer client_data,
XtPointer callData ) {
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Pixel XColorEditor::allocNamedColor(char* colorname,Pixel default_color)
{
//allocate read-only color
  XColor hardwarecolor,exactcolor;
  return (
    XAllocNamedColor(display,xcolormap->getColormap(),colorname,
      &hardwarecolor,&exactcolor) ?
    hardwarecolor.pixel : default_color);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void XColorEditor::exposeCallback(Widget w,XtPointer client_data,
XtPointer callData) {
  if (reinterpret_cast<XmDrawingAreaCallbackStruct*>(callData)->event->
      xexpose.count) return;

  reinterpret_cast<XColorEditor*>(client_data)->expose(w,callData);
}                                
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void XColorEditor::expose(Widget w,XtPointer) {
  Dimension width,height;
  XtVaGetValues(w,XmNwidth,&width,XmNheight,&height,0);

  int num_cmap_colors=xcolormap->getNumberColors();
  int space=1;//space between color rects
  float bar_widthf=float(width-(num_cmap_colors+1)*space)
         /float(num_cmap_colors);
  int bar_width=static_cast<int>(bar_widthf);
  int bar_height=height/3;
  int i, x, y=height/2, wb,xp=0;
  int is=-1,xs,ys,ws;

  Display *w_display=XtDisplay(w);
  Window w_window=XtWindow(w);

  Pixel bg;
  XtVaGetValues ( w,XmNbackground,&bg,0 );
  XSetForeground(w_display,gc,bg);
//erase all
  XFillRectangle(w_display,w_window,gc,0,0,width,height); 
  
  XColor *colors=OPERATOR_NEW_BRACKET(XColor,num_cmap_colors);
  ASSERT(colors!=0);
  Dimension *xc=OPERATOR_NEW_BRACKET(Dimension,num_cmap_colors);
  ASSERT(xc!=0);

  for (i=0,x=space;i<num_cmap_colors;i++,
  x=space+(int)(bar_widthf*i+space*i)) {
    wb=bar_width; 
    if (x-xp>space) {x--;wb++;}
    xp=x+wb;
    colors[i].pixel=xcolormap->getColor(i);
    xc[i]=x+wb/2;
    if (i == selected){
      xs=x;ys=y;ws=wb;is=i;
      continue;
    }
    XSetForeground(w_display,gc,xcolormap->getColor(i));
    XFillRectangle(w_display,w_window,gc,x,y,wb,bar_height);
  }
  if (is!=-1) {
    XSetForeground(w_display,gc,
              allocNamedColor(const_cast<char*>("White"),
              WhitePixel(w_display,DefaultScreen(w_display))));
    XFillRectangle(w_display,w_window,gc,
            xs-3,ys-3,bar_width+6,bar_height+6);
    XSetForeground(w_display,gc,xcolormap->getColor(is));
    XFillRectangle(w_display,w_window,gc,xs,ys,ws,bar_height);
  }
  // draw the RGB curves above the color strip
  XQueryColors(w_display,xcolormap->getColormap(),colors,
    num_cmap_colors);
  int h=height/2-10;
  float r=float(h)/float(65535);

#define ycr(i) 2+h-(int)(r*colors[i].red)
  XSetForeground(w_display,gc,allocNamedColor(const_cast<char*>("red"),
                WhitePixel(w_display,DefaultScreen(w_display))));
  for (i=1;i<num_cmap_colors;i++) {
   XDrawLine(w_display,w_window,gc,xc[i-1],ycr(i-1),xc[i],ycr(i));
  }

#define ycg(i) 2+h-(int)(r*colors[i].green)
  XSetForeground(w_display,gc,allocNamedColor(const_cast<char*>("green"),
                WhitePixel(w_display,DefaultScreen(w_display))));
  for (i=1;i<num_cmap_colors;i++) {
   XDrawLine(w_display,w_window,gc,xc[i-1],ycg(i-1),xc[i],ycg(i));
  }

#define ycb(i) 2+h-(int)(r*colors[i].blue)
  XSetForeground(w_display,gc,allocNamedColor(const_cast<char*>("blue"),
                WhitePixel(w_display,DefaultScreen(w_display))));
  for (i=1;i<num_cmap_colors;i++) {
   XDrawLine(w_display,w_window,gc,xc[i-1],ycb(i-1),xc[i],ycb(i));
  }
  delete []colors; delete []xc;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void XColorEditor::inputEventHandler(Widget w,XtPointer client_data,
XEvent* event,Boolean* /* continue_to_dispatch */ ) {
  reinterpret_cast<XColorEditor*>(client_data)->
    input(w,reinterpret_cast<XtPointer>(event));
}            
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void XColorEditor::inputCallback(Widget w,XtPointer client_data,
XtPointer callData) {
  reinterpret_cast<XColorEditor*>(client_data)->input(w,callData);
}                                
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void XColorEditor::input(Widget w,XtPointer callData ) {
  static unsigned int button=0;// button number
  
  XEvent* event=reinterpret_cast<XEvent*>(callData);
  if (!event) return;
  
  Dimension width,height;
  XtVaGetValues(w,XmNwidth,&width,XmNheight,&height,0);

  int mx, my; //mouse coords
  int i, x, y=height/2, wb,xp=0;
  int space=1;//space between color rects
  int num_cmap_colors=xcolormap->getNumberColors();
  float bar_widthf=float(width-(num_cmap_colors+1)*space)
         /float(num_cmap_colors);
  int bar_width=static_cast<int>(bar_widthf);
  int bar_height=height/3;

  if (event->type == ButtonPress){
    mx=event->xbutton.x; my=event->xbutton.y; 
    button=event->xbutton.button;
    if(my<y || my>y+bar_height) return; // not in any rects
    if(event->xbutton.button != Button1) return;
  } else if (event->type == MotionNotify) {
    mx=event->xmotion.x; my=event->xmotion.y;
    // event->xbutton.button is 0 if cursor is outside the window
    if (my>y) return;
    Colormap cmap=xcolormap->getColormap();
    for (i=0,x=space;i<num_cmap_colors;i++,
    x=space+(int)(bar_widthf*i+space*i)) {

      wb=bar_width; 
      if(x-xp>space){x--;wb++;}
      xp=x+wb;

      if (mx>=x && mx <=x+wb) {
        XColor color; color.pixel=xcolormap->getColor(i);
        XQueryColor(XtDisplay(w), cmap, &color);
        int h=height/2-10;
        float r=float(h)/float(65535);
        int mc=static_cast<int>(float(2+h-my)/r);
        if(mc>65535) mc=65535;
        if(mc<0)mc=0;
        if (button == Button1) { //red
          color.red=mc;
          color.flags =  DoRed;
          XStoreColor(XtDisplay(w),cmap,&color);
        } else if (button == Button2){ //green
          color.green=mc;
          color.flags =  DoGreen;
          XStoreColor ( XtDisplay(w), cmap, &color );
        } else if (button == Button3){ //blue
          color.blue=mc;
          color.flags =  DoBlue;
          XStoreColor ( XtDisplay(w), cmap, &color );
        }
        break;
      }
    }
    expose(w,0); // could only update changed porttion
//  or redraw the line seg's using bg color to erase first.
    return;
  } else if (event->type == ButtonRelease) {
    mx=event->xbutton.x; my=event->xbutton.y;
    return; // do nothing here
  }

  for (i=0,x=space;i<num_cmap_colors;i++,
  x=space+(int)(bar_widthf*i+space*i)) {

    wb=bar_width; 
    if(x-xp>space){x--;wb++;}
    xp=x+wb;

    if (mx>=x && mx <=x+wb) {
      if (selected == i) {  
      } else {
        selected=i;
        XColor color;
        XtVaSetValues(swatch,XmNbackground,xcolormap->getColor(i),0);
        color.flags = DoRed | DoGreen | DoBlue;
        color.pixel = xcolormap->getColor(i);
                XQueryColor ( XtDisplay(w), cmap, &color );
        XtVaSetValues ( red_slider,XmNvalue,color.red/256,0 );
        XtVaSetValues ( green_slider,XmNvalue,color.green/256,0 );
        XtVaSetValues ( blue_slider,XmNvalue, color.blue/256,0 );
      }
//    Instead of calling expose() to redrawl all, we can just redraw 
//    previously selected cell and the newly selected cell and at most
//     2 additional cells around them.
      expose(w,0);
      break;
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Widget XColorEditor::makeSlider(char *name,Widget slider_parent, 
XtCallbackProc callback) {
  Widget w=XtVaCreateManagedWidget(name,xmScaleWidgetClass,
    slider_parent,
    XmNminimum,     0,
    XmNmaximum,     255,
    XmNshowValue,   TRUE,
    XmNorientation, XmHORIZONTAL,
    0 );
  XtAddCallback(w,XmNvalueChangedCallback,callback,this);
  XtAddCallback(w,XmNdragCallback,callback,this);
  return w;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void XColorEditor::colorEditorCallback(Widget w,XtPointer client_data,
XtPointer) {
#if (SPACEDIM>=3)
  return;
#endif
  XColorEditor *ce=reinterpret_cast<XColorEditor*>(client_data);
  if (ce->form) { 
    XtManageChild(ce->form); 
    return;
  }
  
  Arg args[10];
  Cardinal n=0;
  XmString xm_string; 

  n=0;
  char name[80];
  snprintf(name,80,"XColorEditor: %s",
    ce->xcolormap->getPalette()->getName());
  xm_string=XmStringCreateSimple(name);
  XtSetArg(args[n],XmNallowResize,True);n++;
  XtSetArg(args[n],XmNautoUnmanage,False);n++;
  XtSetArg(args[n],XmNdialogTitle,xm_string);n++;
//set colormap for form's shell parent

  Visual *visual=ce->xcolormap->getVisual();
  int depth=ce->xcolormap->getDepth();

  XtSetArg(args[n],XmNvisual,visual);n++; //page 513
  XtSetArg(args[n],XmNdepth,depth);n++;
  XtSetArg(args[n],XmNcolormap,ce->xcolormap->getColormap());n++;
  ce->form=
    XmCreateFormDialog(ce->parent,const_cast<char*>("ColorEdit"),args,n);
  XmStringFree(xm_string);

// create a colorstrip
  Widget frame=
    XtVaCreateManagedWidget("frame",xmFrameWidgetClass,ce->form,
                            XmNtopAttachment,   XmATTACH_FORM,
                            XmNbottomAttachment,XmATTACH_NONE,
                            XmNleftAttachment,  XmATTACH_FORM,
                            XmNrightAttachment, XmATTACH_FORM,
                            XmNleftOffset,      10,
                            XmNtopOffset,       10,
                            XmNbottomOffset,    10,
                            XmNrightOffset,     10,
                            0);

  XtVaCreateManagedWidget("selectorLabel", xmLabelWidgetClass,
                          frame, 
                          XmNchildType, XmFRAME_TITLE_CHILD,
                          0 );

  Widget drawing_area=
    XtVaCreateManagedWidget("drawingArea",
                            xmDrawingAreaWidgetClass,frame, 
                            XmNresizePolicy,XmRESIZE_ANY,
                            XmNwidth, 300,
                            XmNheight, 100,
                            0 ); 

  XtAddCallback(drawing_area,
    XmNresizeCallback,&XColorEditor::resizeCallback,client_data);

  XtAddCallback(drawing_area,
    XmNexposeCallback,&XColorEditor::exposeCallback,client_data);

//  to ask for motion events
  XtAddEventHandler(drawing_area,
    ButtonMotionMask | ButtonPressMask | ButtonReleaseMask, False,
    reinterpret_cast<XtEventHandler>(XColorEditor::inputEventHandler),
    client_data);
                  

  Widget colors=
    XtCreateManagedWidget("frame",xmFrameWidgetClass,ce->form,0,0);

  XtVaCreateManagedWidget("selectorLabel",xmLabelWidgetClass,colors, 
                          XmNchildType, XmFRAME_TITLE_CHILD,
                          0 );

  Widget panel=XtVaCreateManagedWidget("colorpanel",
                                       xmRowColumnWidgetClass, colors, 
                                       XmNradioBehavior, TRUE,
                                       XmNnumColumns,    6,
                                       XmNpacking,      XmPACK_COLUMN,
                                       XmNadjustLast,    FALSE,
                                       0 );
  
  int num_cmap_colors=ce->xcolormap->getNumberColors();
  Colormap cmap=ce->xcolormap->getColormap();
  for (int i=0;i<num_cmap_colors;i++) {
    snprintf(name,80,"%d",i);
    Widget toggle=
      XtVaCreateManagedWidget(name,xmToggleButtonWidgetClass,panel, 
        XmNcolormap,cmap,
        XmNbackground,ce->xcolormap->getColor(i),
        0);
    XtAddCallback(toggle,XmNvalueChangedCallback,selectColorCallback, 
      client_data);
  }

  XtVaSetValues(colors, 
                XmNtopAttachment,    XmATTACH_WIDGET,
                XmNtopWidget,        frame,
                XmNbottomAttachment, XmATTACH_NONE,
                XmNleftAttachment,   XmATTACH_FORM,
                XmNrightAttachment,  XmATTACH_FORM,
                XmNleftOffset,       10,
                XmNtopOffset,        10,
                XmNbottomOffset,     10,
                XmNrightOffset,      10,
                0 );


  Widget controls=XtCreateManagedWidget("controls",xmFormWidgetClass,
                                         ce->form, 0, 0 );

  Widget controls_frame= 
     XtVaCreateManagedWidget("frame", xmFrameWidgetClass, controls,
                             XmNtopAttachment,    XmATTACH_FORM,
                             XmNbottomAttachment, XmATTACH_FORM,
                             XmNleftAttachment,   XmATTACH_NONE,
                             XmNrightAttachment,  XmATTACH_FORM,
                             XmNbottomOffset,     20,
                             XmNrightOffset,      20,
                             0);

 XtVaCreateManagedWidget("swatchLabel", xmLabelWidgetClass, 
                         controls_frame,
                         XmNchildType, XmFRAME_TITLE_CHILD,
                         0 );

 ce->swatch=XtVaCreateManagedWidget("display",
                                xmDrawnButtonWidgetClass,controls_frame,
                                XmNwidth,            100,
                                XmNheight,           100,
                                0 );
  Widget sliders=XtVaCreateManagedWidget("sliderpanel", 
      xmRowColumnWidgetClass, controls, 
      XmNtopAttachment,    XmATTACH_FORM,
      XmNbottomAttachment, XmATTACH_NONE,
      XmNleftAttachment,   XmATTACH_FORM,
      XmNrightAttachment,  XmATTACH_WIDGET,
      XmNrightWidget,      ce->swatch,
      XmNrightOffset,      20,
      0 );

  ce->red_slider=ce->makeSlider(const_cast<char*>("red_slider"),sliders,
    redSliderMovedCallback);
  ce->green_slider=ce->makeSlider(const_cast<char*>("green_slider"),
    sliders,greenSliderMovedCallback);
  ce->blue_slider=ce->makeSlider(const_cast<char*>("blue_slider"),
    sliders,blueSliderMovedCallback);

  XtVaSetValues(controls, 
                XmNtopAttachment,    XmATTACH_WIDGET,
                XmNbottomAttachment, XmATTACH_FORM,
                XmNleftAttachment,   XmATTACH_FORM,
                XmNrightAttachment,  XmATTACH_FORM,
                XmNtopWidget, colors,
                XmNleftOffset,       10,
                XmNtopOffset,        10,
                XmNbottomOffset,     10,
                XmNrightOffset,      10,
                0 );

  XtManageChild(ce->form);

//  set up the colormap //XtParent(w) XtIsShell(w) XtIsTopLevelShell
//  form has to be managed to get a window ID
  XGCValues gv;
  gv.line_width=2;
  unsigned long gmsk=GCLineWidth;
  ce->gc=XCreateGC(XtDisplay(ce->form),XtWindow(drawing_area),
                              gmsk,reinterpret_cast<XGCValues*>(&gv));
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void XColorEditor::printOn(ostream &os) const {
  os << "XColorEditor: display = " << display
     << "screen = " << screen
     << "xcolormap = " << xcolormap
     << "cmap = " << xcolormap->getColormap()
//   << "num_cmap_colors = " << num_cmap_colors
     << endl;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void XColorEditorList::printOn(ostream &os) const {
  os << "\tXColorEditorList:" << endl;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

template class ISLList<XColorEditor>;
#endif
