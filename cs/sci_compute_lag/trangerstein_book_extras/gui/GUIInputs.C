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
// "$Header: /home/faculty/johnt/cvs/deal_new/gui/GUIInputs.C,v 1.1 2009/08/20 17:32:35 johnt Exp $"
#include "GUIInputs.H"
#include <Xm/TextF.h>
#include "MemoryDebugger.H"
#include "Tracer.H"
//:::::::::::::::::::::::: InputString :::::::::::::::::::::::::::::::::
void GUIInputString::valueChanged(GUI_WIDGET_POINTER text_widget) {
  GUI_INPUTS_BASE::valueChanged(text_widget);
#ifdef USE_GTK
  const gchar *text=gtk_entry_get_text(GTK_ENTRY(text_widget));
#else
  char *text=XmTextFieldGetString(text_widget);  
#endif
  if (strlen(text)<maxlen) {
    if (ptr_string[0]) delete ptr_string[0];
    length_string=strlen(text)+1;
    ptr_string[0]=OPERATOR_NEW_BRACKET(char,length_string);
    strcpy(ptr_string[0],text);
  }
  writeTo(text_widget);
#ifndef USE_GTK
  XtFree(text);
#endif
}
//::::::::::::::: InputIFStream ::::::::::::::::::::::::::::::::::::::::
void GUIInputIFStream::valueChanged(GUI_WIDGET_POINTER text_widget) {
  GUI_INPUTS_BASE::valueChanged(text_widget);
#ifdef USE_GTK
  const gchar *text=gtk_entry_get_text(GTK_ENTRY(text_widget));
#else
  char *text=XmTextFieldGetString(text_widget);
#endif
  changeString(text);
  writeTo(text_widget);
#ifndef USE_GTK
  XtFree(text);
#endif
}

//::::::::::::::::::::::::: InputOFStream ::::::::::::::::::::::::::::::
void GUIInputOFStream::valueChanged(GUI_WIDGET_POINTER text_widget) {
  GUI_INPUTS_BASE::valueChanged(text_widget);
#ifdef USE_GTK
  const gchar *text=gtk_entry_get_text(GTK_ENTRY(text_widget));
#else
  char *text=XmTextFieldGetString(text_widget);
#endif
  changeString(text);
  writeTo(text_widget);
#ifndef USE_GTK
  XtFree(text);
#endif
}
