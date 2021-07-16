// Copyright (C) 2004 Prof. John A. Trangenstein, Duke University
// Modified by Wenjun Ying (Tue May  4 11:30:09 EDT 2004)

#include "GUIInputParameter.H"
#include <Xm/RowColumn.h>
#include <Xm/TextFP.h>
#include <Xm/ToggleB.h>
#include "Tracer.H"

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Widget GUIInputParameter<bool>::createWidget(Widget parent) {
  Widget radiobox = XtVaCreateWidget("radiobox",
                                xmRowColumnWidgetClass,
                                parent,
                                XmNorientation, XmHORIZONTAL,
                                XmNpacking, XmPACK_COLUMN,
                                XmNradioBehavior, True,
                                XmNradioAlwaysOne, True,
                                0 );

  XmString xm_string = XmStringCreateSimple("True");
  Widget toggle1 = XtVaCreateManagedWidget ( "toggle1",
                               xmToggleButtonWidgetClass,
                               radiobox,
                               XmNlabelType, XmSTRING,
                               XmNlabelString, xm_string,
                               0 );
  XmStringFree(xm_string);

  XtAddCallback(toggle1,
                XmNvalueChangedCallback,
                &GUIVirtualInput::valueChangedCallback,
                reinterpret_cast<XtPointer>(this) );

  xm_string = XmStringCreateSimple("False");
  Widget toggle2 = XtVaCreateManagedWidget  ( "toggle2",
                               xmToggleButtonWidgetClass,
                               radiobox,
                               XmNlabelType, XmSTRING,
                               XmNlabelString, xm_string,
                               0 );
  XmStringFree(xm_string);

  XtAddCallback(toggle2,
                XmNvalueChangedCallback,
                &GUIVirtualInput::valueChangedCallback,
                reinterpret_cast<XtPointer>(this) );

  if (*ptr_data) {
    XmToggleButtonSetState(toggle1, True, False);
  } else {
    XmToggleButtonSetState(toggle2, True, False);
  }

  return radiobox;
}

void GUIInputParameter<bool>::valueChanged(Widget w)
{
  value_changed = TRUE;
  Widget toggle = w;
  if (!XmToggleButtonGetState(toggle)) {
    return;
  }
  XmString xm_string;
  char *label;
  XtVaGetValues(toggle, XmNlabelString,&xm_string, 0);
  XmStringGetLtoR(xm_string, XmFONTLIST_DEFAULT_TAG,&label);
  XmStringFree(xm_string);
  if (!strcmp(label, "True")) {
    *ptr_data = 1;
  } else {
    *ptr_data = 0;
  }
  XtFree(label);
  return;
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GUIInputParameter<int>::valueChanged(Widget w)
{
  value_changed = TRUE;
  char *text = XmTextFieldGetString(w);
  int t = atoi(text);
  if (t >= lower_bound && t <= upper_bound) {
    *ptr_data = t;
  } else {
    XBell(XtDisplay(w), 40);
  }
  writeTo(w);
  XtFree(text);
}

Widget GUIInputParameter<int>::createWidget(Widget parent)
{
  return GUIVirtualInput::createWidget(parent);
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GUIInputParameter<double>::valueChanged(Widget w)
{
  value_changed = TRUE;
  char *text = XmTextFieldGetString(w);
  double t = atof(text);
  if (t >= lower_bound && t <= upper_bound) {
    *ptr_data = t;
  } else {
    XBell(XtDisplay(w), 40);
  }
  writeTo(w);
  XtFree(text);
}

Widget GUIInputParameter<double>::createWidget(Widget parent)
{
  return GUIVirtualInput::createWidget(parent);
}

template class GUIInputParameter<bool>;
template class GUIInputParameter<int>;
template class GUIInputParameter<double>;
