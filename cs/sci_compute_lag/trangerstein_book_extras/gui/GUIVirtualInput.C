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
// "$Header: /home/faculty/johnt/cvs/deal_new/gui/GUIVirtualInput.C,v 1.1 2009/08/20 17:32:35 johnt Exp $"
#include "GUIVirtualInput.H"
#include <Xm/Xm.h>
#include <Xm/ArrowB.h>
#include <Xm/Form.h>
#include <Xm/Label.h>
#include <Xm/MessageB.h>
#include <Xm/PushB.h>
#include <Xm/RowColumn.h>
#include <Xm/ScrolledW.h>
#include <Xm/Separator.h>
#include <Xm/TextF.h>
#include <Xm/ToggleB.h>
#include <Xm/XmStrDefs.h>
#include "MemoryDebugger.H"
#include "Tracer.H"
#include "NISLList.C"
//::::::::::::::::::::::: VirtualInput :::::::::::::::::::::::::::::::::
GUIVirtualInput::GUIVirtualInput(const char *g) : group_name(0),
group_arrow(0),group_label(0),label(0),textfield(0),value_changed(FALSE)
{
  group_name=OPERATOR_NEW_BRACKET(char,strlen(g)+1);
  strcpy(group_name,g);
}

GUIVirtualInput::~GUIVirtualInput() {
  delete group_name; group_name=0;
}

void GUIVirtualInput::writeTo(Widget tf) const {
  XtVaSetValues(tf,XmNvalue,getValue(),0);
}

void GUIVirtualInput::valueChangedCallback(Widget w,
XtPointer client_data,XtPointer) {
  reinterpret_cast<GUIVirtualInput*>(client_data)->valueChanged(w);
}

void printConstraintResources(ostream &os,const Widget constraint) {
  Widget topWidget;
  XtVaGetValues(constraint,
    XmNtopWidget,&topWidget,
    0);
  os << "Constraint Resources for " << constraint
     << " XmNtopWidget = " << topWidget
     << endl;
}

void printCoreResources(ostream &os,const Widget core) {
  XtAccelerators accelerators;
  Boolean ancestor_sensitive,initial_resources_persistent=True,
          mapped_when_managed=True,sensitive=True;
  Pixel background,border_color;
  Pixmap background_pixmap=XmUNSPECIFIED_PIXMAP,
         border_pixmap=XmUNSPECIFIED_PIXMAP;
  Dimension border_width=1,height,width;
  Colormap colormap;
  int depth;
  XtCallbackList destroy_callback=0;
  Screen *screen;
  XtTranslations translations;
  Position x=0,y=0;

  XtVaGetValues(core, 
    XmCAccelerators, &accelerators,
    XmNancestorSensitive, &ancestor_sensitive,
    XmNbackground, &background,
    XmNbackgroundPixmap, &background_pixmap,
    XmNborderColor, &border_color,
    XmNborderPixmap, &border_pixmap,
    XmNborderWidth, &border_width,
    XmNcolormap, &colormap,
    XmNdepth, &depth,
    XmNdestroyCallback, &destroy_callback,
    XmNheight, &height,
    XmNinitialResourcesPersistent, &initial_resources_persistent,
    XmNmappedWhenManaged, &mapped_when_managed,
    XmNscreen, &screen,
    XmNsensitive, &sensitive,
    XmNtranslations, &translations,
    XmNwidth, &width,
    XmNx, &x,
    XmNy, &y,
    0);

  os << "Core Resources for " << core 
     << " isWidget = " << XtIsWidget(core) 
     << " isManaged = " << XtIsManaged(core) 
     << endl;
  if (XtIsWidget(core)) {
    os 
//     << "\tancestor_sensitive = " << ancestor_sensitive << "\n"
//     << "\tbackground = " << background << "\n"
//     << "\tborder_color = " << border_color << "\n"
//     << "\tborder_width = " << border_width << "\n"
//     << "\tdepth = " << depth << "\n"
//     << "\tdestroy_callback = " << destroy_callback << "\n"
//     << "\theight = " << height << "\n"
//     << "\tinitial_resources_persistent = " 
//        << initial_resources_persistent << "\n"
//     << "\tmapped_when_managed = " 
//        << mapped_when_managed << "\n"
//     << "\tscreen = " << screen << "\n"
//     << "\ttranslations = " << translations << "\n"
//     << "\twidth = " << width << "\n"
       << "\tx = " << x << "\n"
       << "\ty = " << y 
       << endl;
  }

}

void printPrimitiveResources(ostream &os,const Widget primitive) {
  Pixel bottom_shadow_color,foreground,highlight_color,top_shadow_color;
  Pixmap bottom_shadow_pixmap=XmUNSPECIFIED_PIXMAP,highlight_pixmap,
    top_shadow_pixmap;
  Boolean highlight_on_enter=False,traversal_on=True;
  Dimension highlight_thickness=2,shadow_thickness=2;
  XmNavigationType navigation_type=XmNONE;
  unsigned char unit_type;
  XtPointer user_data=0;

  XtVaGetValues(primitive, 
    XmNbottomShadowColor, &bottom_shadow_color,
    XmNbottomShadowPixmap, &bottom_shadow_pixmap,
    XmNforeground, &foreground,
    XmNhighlightColor, &highlight_color,
    XmNhighlightOnEnter, &highlight_on_enter,
    XmNhighlightPixmap, &highlight_pixmap,
    XmNhighlightThickness, &highlight_thickness,
    XmNnavigationType, &navigation_type,
    XmNshadowThickness, &shadow_thickness,
    XmNtopShadowColor, &top_shadow_color,
    XmNtopShadowPixmap, &top_shadow_pixmap,
    XmNtraversalOn, &traversal_on,
    XmNunitType, &unit_type,
    XmNuserData, &user_data,
    0);

  os << "Primitive Resources for " << primitive 
     << " isPrimitive = " << XmIsPrimitive(primitive)
     << endl; 
  if (XmIsPrimitive(primitive)) {
    os << "\n\tbottom_shadow_color = " << bottom_shadow_color
       << "\n\tbottom_shadow_pixmap = " << bottom_shadow_pixmap
       << "\n\tforeground = " << foreground
       << "\n\thighlight_color = " << highlight_color
       << "\n\thighlight_on_enter = " << highlight_on_enter
       << "\n\thighlight_pixmap = " << highlight_pixmap
       << "\n\thighlight_thickness = " << highlight_thickness
       << "\n\tnavigation_type = " << navigation_type
       << "\n\tshadow_thickness = " << shadow_thickness
       << "\n\ttop_shadow_color = " << top_shadow_color
       << "\n\ttop_shadow_pixmap = " << top_shadow_pixmap
       << "\n\ttraversal_on = " << traversal_on
       << "\n\tunit_type = " << unit_type
       << "\n\tuser_data = " << user_data
       << endl;
  }
}

void printLabelResources(ostream &os,const Widget label) {
  String accelerator=0,mnemonic_char_set=0;
  XmString accelerator_text=0,label_string=0;
  unsigned char alignment,label_type;
  XmFontList font_list;
  Pixmap label_insensitive_pixmap=XmUNSPECIFIED_PIXMAP,
         label_pixmap=XmUNSPECIFIED_PIXMAP;
  Dimension margin_bottom,margin_height,margin_left,margin_right,
    margin_top,margin_width;
  KeySym mnemonic;
  Boolean recompute_size=True;
  XmStringDirection string_direction;

  XtVaGetValues(label, 
    XmNaccelerator, &accelerator,
    XmNacceleratorText, &accelerator_text,
    XmNalignment, &alignment,
    XmNfontList, &font_list,
    XmNlabelInsensitivePixmap, &label_insensitive_pixmap,
    XmNlabelPixmap, &label_pixmap,
    XmNlabelString, &label_string,
    XmNlabelType, &label_type,
    XmNmarginBottom, &margin_bottom,
    XmNmarginHeight, &margin_height,
    XmNmarginLeft, &margin_left,
    XmNmarginRight, &margin_right,
    XmNmarginTop, &margin_top,
    XmNmarginWidth, &margin_width,
    XmNmnemonic, &mnemonic,
    XmNmnemonicCharSet, &mnemonic_char_set,
    XmNrecomputeSize, &recompute_size,
    XmNstringDirection, &string_direction,
    0);

  os << "Label Resources for " << label 
     << " isLabel = " << XmIsLabel(label)
     << endl; 
  if (XmIsLabel(label)) {
    os << "\n\taccelerator = "<< accelerator
       << "\n\taccelerator_text = " << accelerator_text
       << "\n\talignment = " << alignment
       << "\n\tlabel_insensitive_pixmap = " << label_insensitive_pixmap
       << "\n\tlabel_pixmap = " << label_pixmap
       << "\n\tlabel_string = " << label_string
       << "\n\tlabel_type = " << label_type
       << "\n\tmargin_bottom = " << margin_bottom 
       << "\n\tmargin_height = " << margin_height 
       << "\n\tmargin_left = " << margin_left 
       << "\n\tmargin_right = " << margin_right 
       << "\n\tmargin_top = " << margin_top 
       << "\n\tmargin_width = " << margin_width 
       << "\n\tmnemonic = " << mnemonic 
       << "\n\tmnemonic_char_set = " << mnemonic_char_set 
       << "\n\trecompute_size = " << recompute_size
       << "\n\tstring_direction = " << string_direction
       << endl;
  }
}

void printTextFieldResources(ostream &os,const Widget text_field) {
  int blink_rate=500,max_length,select_threshold=5;
  short columns;
  XmTextPosition cursor_position=0;
  Boolean cursor_position_visible=True,editable=True,
          pending_delete=True,resize_width=False,verify_bell;
  XmFontList font_list;
  Dimension margin_height=5,margin_width=5;
  String value;
  wchar_t *value_wcs;

  XtVaGetValues(text_field, 
    XmNblinkRate, &blink_rate,
    XmNcolumns, &columns,
    XmNcursorPosition, &cursor_position,
    XmNcursorPositionVisible, &cursor_position_visible,
    XmNeditable, &editable,
    XmNfontList, &font_list,
    XmNmarginHeight, &margin_height,
    XmNmarginWidth, &margin_width,
    XmNmaxLength, &max_length,
    XmNpendingDelete, &pending_delete,
    XmNresizeWidth, &resize_width,
    XmNselectThreshold, &select_threshold,
    XmNvalue, &value,
    XmNvalueWcs, &value_wcs,
    XmNverifyBell, &verify_bell,
    0);

  os << "TextField Resources for " << text_field 
     << " isTextField = " << XmIsTextField(text_field)
     << endl; 
  if (XmIsTextField(text_field)) {
    os << "\n\tblink_rate = "<< blink_rate
       << "\n\tcolumns = " << columns
       << "\n\tcursor_position = " << cursor_position
       << "\n\tcursor_position_visible = " 
          << cursor_position_visible
       << "\n\teditable = " << editable
       << "\n\tfont_list = " << font_list
       << "\n\tmargin_height = " << margin_height
       << "\n\tmargin_width = " << margin_width
       << "\n\tmax_length = " << max_length 
       << "\n\tpending_delete = " << pending_delete
       << "\n\tresize_width = " << resize_width
       << "\n\tselect_threshold = " << select_threshold 
       << "\n\tvalue = " << value 
       << "\n\tvalue_wcs = " << value_wcs 
       << "\n\tverify_bell = " << verify_bell
       << endl;
  }
}

void printManaged(ostream &os,Widget widget) {
  os << "\twidget = " << widget
     << " isManaged = " << XtIsManaged(widget)
     << " isRealized = " << XtIsRealized(widget) << endl;
  Widget parent=XtParent(widget);
  os << "\t\tparent = " << parent  
     << " isManaged = " << XtIsManaged(parent)
     << " isRealized = " << XtIsRealized(parent)
     << endl;
}

void printFamilyManaged(ostream &os,Widget widget) {
  Widget w=widget;
  int count=0;
  while (w!=0 && count<10) {
    os << "\tWidget " << w
       << (XtIsManaged(w) ? " is managed and " : " is NOT managed and ")
       << (XtIsRealized(w) ? " is realized" : " is NOT realized") << endl;
    count++;
    w=XtParent(w);
  }
}

void GUIVirtualInput::printOn(ostream &os) const {
  os << "GUIVirtualInput: name = " << getName()
     << "\n\tgroup_name = " << group_name << endl;
  if (label!=0) {
    os << "\tlabel:" << endl;
    printManaged(os,label);
//  printLabelResources(os,label);
//  printPrimitiveResources(os,label);
    printCoreResources(os,label);
  }

  if (textfield!=0) {
    os << "\n\ttextfield:" << endl;
    printManaged(os,textfield);
//  printTextFieldResources(os,textfield);
//  printPrimitiveResources(os,textfield);
    printCoreResources(os,textfield);
  }

  if (group_arrow!=0) {
    os << "\n\tgroup_arrow:" << endl;
    printManaged(os,group_arrow);
    printCoreResources(os,group_arrow);
  }

  if (group_label!=0) {
    os << "\n\tgroup_label:" << endl;
    printManaged(os,group_label);
    printCoreResources(os,group_label);
  }
}

Widget GUIVirtualInput::createWidget(Widget parent) {
  Widget textfield=XtVaCreateWidget("textfield",
                                    xmTextFieldWidgetClass,
                                    parent, 
                                    XmNwidth, 251, 
                                    0); 
  XtAddCallback(textfield,XmNactivateCallback,
                &GUIVirtualInput::valueChangedCallback,
                reinterpret_cast<XtPointer>(this));
  XtAddCallback(textfield,XmNlosingFocusCallback,
                &GUIVirtualInput::valueChangedCallback,
                reinterpret_cast<XtPointer>(this));
  XtVaSetValues(textfield,XmNvalue,getValue(),0);
  return textfield;
}

ostream& operator<<(ostream &os,const GUIVirtualInput& vi) {
  os << "GUIVirtualInput: group_name = " << vi.getGroupName() 
     << ", value = " << vi.getValue();
  return os;
}

void GUIVirtualInput::reload() {
  XtVaSetValues(textfield, XmNvalue, getValue(), 0 );
}

void GUIVirtualInput::manage(unsigned char adir) {
  if (adir == XmARROW_RIGHT) { // show all params in group
    XtManageChild(label);
    XtManageChild(textfield);
  } else {
    XtUnmanageChild(label);
    XtUnmanageChild(textfield);
  }
}

void GUIVirtualInput::setTopWidget(Widget w) {
  if (group_arrow==0 || group_label==0) return;
  XtUnmanageChild(group_arrow);
  XtUnmanageChild(group_label);
  XtVaSetValues(group_arrow, XmNtopWidget, w, 0);
  XtVaSetValues(group_label, XmNtopWidget, w, 0);
  XtManageChild(group_arrow);
  XtManageChild(group_label);
}

//:::::::::::::: InputParameterList ::::::::::::::::::::::::::::::::::::
bool GUIInputParameterList::badparams=FALSE;
GUIInputParameterList::GUIInputParameterList(const char *n) :
InputParameterList(n), form_dialog_widget(0), warn(0), //needreload(FALSE),
warn_up(FALSE), verify(0), verify_data(0) {
}

//called by done_button
void GUIInputParameterList::hideCallback(Widget,XtPointer client_data,
XtPointer) {
  reinterpret_cast<GUIInputParameterList*>(client_data)->hide();
}

void GUIInputParameterList::hide() {
  bool bad=badparams;
  badparams=FALSE;

  if (verify) verify(0,verify_data,0);

  if (!badparams){
    XtUnmanageChild(form_dialog_widget);
    badparams=bad;
  } else if (warn_up) showWarningDialog(0);
}

void GUIInputParameterList::arrowActCallback(Widget group_arrow,
XtPointer client_data,XtPointer) {
  GUIInputParameterList* list=
    reinterpret_cast<GUIInputParameterList*>(client_data);
  Widget form=XtParent(group_arrow);
  XtUnmanageChild(form);

  unsigned char adir=XmARROW_RIGHT; 
  XtVaGetValues(group_arrow,XmNarrowDirection,&adir,0);
  if (adir == XmARROW_RIGHT) { // to show all params in the group
    XtVaSetValues(group_arrow, XmNarrowDirection, XmARROW_DOWN, 0);
  } else {
    XtVaSetValues(group_arrow, XmNarrowDirection, XmARROW_RIGHT, 0);
  }
  
  NISLListNode<VirtualInput> *p=0;
  GUIVirtualInput *gp=0;
  for (p=list->first();p;p=list->next(p)) {
    gp=dynamic_cast<GUIVirtualInput*>(p->getData());
    if (gp->isGroupArrowWidget(group_arrow)) break;
  }
  CHECK_POINTER(p);

//manage or unmanage all label and textfield in this group
  GUIVirtualInput *gl=gp;
  for (NISLListNode<VirtualInput> *q=p;q;q=list->next(q)) {
    GUIVirtualInput *gq=dynamic_cast<GUIVirtualInput*>(q->getData());
    if (gq->sameGroupNameAs(gp)) {
      gq->manage(adir);
      gl=gq;
    }
  }
//find top attachment for next group
  Widget tw=(adir==XmARROW_RIGHT) ? gl->getTextField()
                                  : gp->getGroupLabel();

//set top attachment for next group arrow and label
  for (NISLListNode<VirtualInput> *q=list->next(p);q;q=list->next(q)) {
    GUIVirtualInput *gq=dynamic_cast<GUIVirtualInput*>(q->getData());
    if (!gq->sameGroupNameAs(gp)) {
      bool q_in_previous_group=FALSE;
      for (NISLListNode<VirtualInput> *r=list->first();r!=q;
      r=list->next(r)) {
        GUIVirtualInput *gr=
          dynamic_cast<GUIVirtualInput*>(r->getData());
        if (gr->sameGroupNameAs(gq)) { 
          q_in_previous_group=TRUE; 
          break; 
        }
      }
      if (q_in_previous_group) continue;
      gq->setTopWidget(tw);
      break;
    }
  }
  XtManageChild(form);
}

void GUIInputParameterList::createListWin(Widget parent,
void (*verify_callback)(Widget,XtPointer,XtPointer),XtPointer data,
const char *verify_button_label) {
  CHECK_TEST(form_dialog_widget==0);

//FormDialog widget hold everything for the InputParameterList 
  XmString xm_string=XmStringCreateSimple(getName());
  Arg args[5];
  Cardinal n=0;
  XtSetArg(args[n],XmNallowResize,TRUE); n++;
  XtSetArg(args[n],XmNautoUnmanage,FALSE); n++;
  XtSetArg(args[n],XmNdialogTitle,xm_string); n++;
  form_dialog_widget=XmCreateFormDialog(parent,getName(),args,n);
  XmStringFree(xm_string);

//PushButton done_button hides the list window
  Widget done_button=XmCreatePushButton(form_dialog_widget,
    const_cast<char*>("Done"),0,0);
  XtAddCallback(done_button,XmNactivateCallback,
                 &GUIInputParameterList::hideCallback, 
                 reinterpret_cast<XtPointer>(this));
  XtManageChild(done_button);  
  XtVaSetValues(done_button,
                XmNtopAttachment, XmATTACH_FORM,
                0);

  if (verify_callback) {
//  PushButton verify_button calls routine to check input parameters
    verify=verify_callback; 
    verify_data=data;
    Widget verify_button=0;
    if (verify_button_label!=0) {
      verify_button=
        XmCreatePushButton(form_dialog_widget,
          const_cast<char*>(verify_button_label),0,0);
    } else {
      verify_button=XmCreatePushButton(form_dialog_widget,
        const_cast<char*>("Verify"),0,0);
    }
    XtAddCallback(verify_button,XmNactivateCallback,verify_callback,
                  data);
    XtManageChild(verify_button);  
    XtVaSetValues(verify_button,
                  XmNtopAttachment, XmATTACH_FORM,
                  XmNleftAttachment,XmATTACH_WIDGET,
                  XmNleftWidget, done_button,
                  0);
  }

//SeparatorWidget separator separates done_button from input parameters
  Widget separator=XtVaCreateManagedWidget("separator",
                     xmSeparatorWidgetClass,
                     form_dialog_widget, 
                     XmNresizePolicy, XmRESIZE_NONE, 
                     XmNtopAttachment, XmATTACH_WIDGET, 
                     XmNbottomAttachment, XmATTACH_NONE, 
                     XmNleftAttachment, XmATTACH_FORM, 
                     XmNrightAttachment, XmATTACH_FORM, 
                     XmNtopWidget, done_button, 
                     0); 

//ScrolledWindow scrolled_window contains form holding the input parameters
  Widget scrolled_window=XtVaCreateManagedWidget("scrolled",
                           xmScrolledWindowWidgetClass,
                           form_dialog_widget,
                           XmNresizePolicy, XmRESIZE_GROW, 
                           XmNscrollingPolicy, XmAUTOMATIC,
                           XmNtopAttachment, XmATTACH_WIDGET,
                           XmNtopWidget,separator,
                           XmNbottomAttachment, XmATTACH_FORM,
                           XmNleftAttachment,XmATTACH_FORM,
                           XmNrightAttachment, XmATTACH_FORM,
                           XmNwidth,600,XmNheight,300,
                           0);

//Form form contains the input parameter groups
  Widget form=XtVaCreateWidget("form",
                               xmFormWidgetClass,
                               scrolled_window, 
                               XmNresizePolicy, XmRESIZE_GROW, 
                               0 ); 

  XmFontList fontlist;
  XtVaGetValues(form,XmNlabelFontList,&fontlist,0);

  NISLList<Widget> arrow_widget_list;
  Widget previous_group_arrow=0,previous_group_label=0;
  for (NISLListNode<VirtualInput> *p=first();p;p=next(p)) {
    GUIVirtualInput *gp=dynamic_cast<GUIVirtualInput*>(p->getData());
    CHECK_POINTER(gp);
//  look for VirtualInput before p in same group as p
    bool p_group_previously_processed=FALSE;
    for (NISLListNode<VirtualInput> *q=first();q!=p;q=next(q)) {
      GUIVirtualInput *gq=dynamic_cast<GUIVirtualInput*>(q->getData());
      if (gq->sameGroupNameAs(gp)) { 
        p_group_previously_processed=TRUE; 
        break; 
      }
    }
    if (p_group_previously_processed) continue;
//  new group: set group label
    char *group_name=gp->getGroupName();
    char group_label[80];
    if (strcmp(group_name,"")==0) strcpy(group_label,"Misc Params");
    else strcpy(group_label,group_name);
//  ArrowButton current_group_arrow is managed from the beginning
    Widget current_group_arrow=XtVaCreateManagedWidget(group_name,
                          xmArrowButtonWidgetClass,
                          form, 
                          XmNarrowDirection,XmARROW_RIGHT, 
                          XmNwidth,30, 
                          XmNheight,30, 
                          0); 
    arrow_widget_list.append(OPERATOR_NEW Widget(current_group_arrow));
//  Label current_group_label is managed from the beginning
    xm_string=XmStringCreateSimple(group_label);
    Widget current_group_label=XtVaCreateManagedWidget(group_name,
                          xmLabelWidgetClass,
                          form, 
                          XmNalignment,XmALIGNMENT_BEGINNING, 
                          XmNlabelType,XmSTRING, 
                          XmNlabelString,xm_string,
                          XmNwidth,300, 
                          XmNheight,30, 
                          0); 
    XmStringFree(xm_string);

    if (previous_group_arrow==0) {
//    first group to be processed; attach to form
      XtVaSetValues(current_group_arrow,
              XmNtopAttachment, XmATTACH_FORM, 
              XmNbottomAttachment, XmATTACH_NONE, 
              XmNleftAttachment, XmATTACH_FORM, 
              XmNrightAttachment, XmATTACH_NONE, 
              0 );
      XtVaSetValues(current_group_label,
              XmNtopAttachment, XmATTACH_FORM, 
              XmNbottomAttachment, XmATTACH_NONE, 
              XmNleftAttachment, XmATTACH_WIDGET, 
              XmNrightAttachment, XmATTACH_FORM, 
              XmNleftWidget, current_group_arrow,
              0 );
    } else {
//    not first group to be processed; attach to previous group arrow, label
      XtVaSetValues(current_group_arrow,
              XmNtopAttachment, XmATTACH_WIDGET, 
              XmNbottomAttachment, XmATTACH_NONE, 
              XmNleftAttachment, XmATTACH_FORM, 
              XmNrightAttachment, XmATTACH_OPPOSITE_WIDGET, 
              XmNtopWidget, previous_group_arrow, 
              XmNrightWidget, previous_group_arrow, 
              0 );
      XtVaSetValues(current_group_label,
              XmNtopAttachment, XmATTACH_WIDGET, 
              XmNbottomAttachment, XmATTACH_NONE, 
              XmNleftAttachment, XmATTACH_WIDGET, 
              XmNrightAttachment, XmATTACH_FORM, 
              XmNtopWidget, previous_group_arrow,
              XmNleftWidget, current_group_arrow, 
              0 );
    }
//  arrowActCallback expands or contracts parameter list
    XtAddCallback(current_group_arrow,
                  XmNactivateCallback,
                  &GUIInputParameterList::arrowActCallback,
                  reinterpret_cast<XtPointer>(this) );  
//  find max width of labels in this group
    Dimension labelwidth=100;
    for (NISLListNode<VirtualInput> *q=p;q;q=next(q)) {
      GUIVirtualInput *gq=dynamic_cast<GUIVirtualInput*>(q->getData());
      CHECK_POINTER(gq);
      if (gq->sameGroupNameAs(gp)) {
        xm_string=XmStringCreateSimple(gq->getName());
        Dimension w=XmStringWidth(fontlist,xm_string);
        XmStringFree(xm_string);
        if (w>labelwidth) labelwidth=w;
      }
    }

//  put all input parameters widgets in this group into form widget
//    (managed!)
    Widget previous_label_in_group=0,previous_textfield_in_group=0;
    for (NISLListNode<VirtualInput> *q=p;q;q=next(q)) {
      GUIVirtualInput *gq=dynamic_cast<GUIVirtualInput*>(q->getData());
      CHECK_POINTER(gq);
      if (!gq->sameGroupNameAs(gp)) continue;
      xm_string=XmStringCreateSimple(gq->getName());

      Widget gq_textfield=gq->createWidget(form); 
      Widget gq_label=0;
      if (q==p) { // The first one
        gq_label=XtVaCreateWidget("label01",
                                  xmLabelWidgetClass,
                                  form, 
                                  XmNalignment,XmALIGNMENT_END,
                                  XmNlabelType, XmSTRING, 
                                  XmNlabelString,xm_string,
                                  XmNtopAttachment, XmATTACH_WIDGET, 
                                  XmNbottomAttachment, XmATTACH_NONE, 
                                  XmNleftAttachment, XmATTACH_WIDGET,
                                  XmNrightAttachment, XmATTACH_NONE, 
                                  XmNtopWidget, current_group_label, 
                                  XmNleftWidget, current_group_arrow, 
                                  XmNtopOffset, 5, 
                                  XmNwidth, labelwidth, 
                                  0); 
        XtVaSetValues(gq_textfield,
                      XmNtopAttachment, XmATTACH_WIDGET, 
                      XmNbottomAttachment, XmATTACH_NONE, 
                      XmNleftAttachment, XmATTACH_WIDGET, 
                      XmNrightAttachment, XmATTACH_FORM, 
                      XmNtopWidget, current_group_label, 
                      XmNleftWidget,gq_label,
                      XmNwidth, 251, 
                      0 ); 
      } else {
        gq_label=XtVaCreateWidget("label02",
                                  xmLabelWidgetClass,
                                  form, 
                                  XmNalignment,XmALIGNMENT_END,
                                  XmNlabelType, XmSTRING, 
                                  XmNlabelString,xm_string,
                                  XmNtopAttachment, XmATTACH_WIDGET, 
                                  XmNbottomAttachment, XmATTACH_NONE, 
                                  XmNleftAttachment, XmATTACH_WIDGET, 
                                  XmNrightAttachment,XmATTACH_NONE,
                                  XmNtopWidget,previous_textfield_in_group, 
                                  XmNleftWidget, current_group_arrow,
                                  XmNtopOffset, 5, 
                                  XmNwidth, labelwidth, 
                                  0 );   
        XtVaSetValues(gq_textfield,
                      XmNtopAttachment, XmATTACH_WIDGET, 
                      XmNbottomAttachment, XmATTACH_NONE, 
                      XmNleftAttachment, XmATTACH_WIDGET, 
                      XmNrightAttachment, XmATTACH_FORM, 
                      XmNtopWidget,previous_textfield_in_group,
                      XmNleftWidget,gq_label,
                      XmNwidth, 251, 
                      0 ); 
      }
      XmStringFree(xm_string);
      gq->setLabel(gq_label);
      gq->setTextField(gq_textfield);
      gq->setGroupArrow(current_group_arrow);
      gq->setGroupLabel(current_group_label);
      previous_label_in_group=gq_label;
      previous_textfield_in_group=gq_textfield;                      
    }
    previous_group_arrow=current_group_arrow;
    previous_group_label=current_group_label;
  }
  XtManageChild(form);

  while (arrow_widget_list.notEmpty()) delete arrow_widget_list.delAfter(0);
}

void GUIInputParameterList::showWarningDialog(const char* str) {
  if (str) badparams=TRUE;
  if (form_dialog_widget==0) return;

  XmString xstr;
  if (str) {
    const char *def=XmSTRING_DEFAULT_CHARSET;
    xstr=XmStringCreateLtoR(const_cast<char*>(str),const_cast<char*>(def));
  }
    
  if (warn) {
    if (str) XtVaSetValues(warn,XmNmessageString, xstr, 0);
    XtManageChild(form_dialog_widget);
    XtManageChild(warn); 
  } else {
    Arg args[10];
    int n=0;
    if (str) { XtSetArg(args[n], XmNmessageString, xstr); n++; }
    warn=XmCreateWarningDialog(form_dialog_widget,
      const_cast<char*>("warning"),args,n);
  }
  if(str) XmStringFree(xstr);
  XtManageChild(form_dialog_widget);
  XtManageChild(warn);
  warn_up=TRUE;
}

void GUIInputParameterList::reloadListWin() {
  CHECK_TEST(form_dialog_widget!=0);
  for (NISLListNode<VirtualInput> *p=first();p;p=next(p)) {
    dynamic_cast<GUIVirtualInput*>(p->getData())->reload();
  }
  XtManageChild(form_dialog_widget);
}

bool GUIInputParameterList::getValueChanged() {
  if (form_dialog_widget==0) { return TRUE; }
  
  for (NISLListNode<VirtualInput> *p=first();p;p=next(p)) {
    if (dynamic_cast<GUIVirtualInput*>(p->getData())->valueChanged())
      return TRUE;
  }
  return FALSE;
}

void GUIInputParameterList::unsetValueChanged() {
  if (form_dialog_widget) return;
  for (NISLListNode<VirtualInput> *p=first();p;p=next(p)) {
    dynamic_cast<GUIVirtualInput*>(p->getData())->unsetValueChanged();
  }
}

void GUIInputParameterList::printOn(ostream &os) const {
  os << "GUIInputParameterList :"
//   << "\tneedreload = " << needreload
     << "\tbadparams = " << badparams
     << "\twarn_up = " << warn_up << endl;
  InputParameterList::printOn(os);
}

template class NISLListNode<Widget>;
template ostream& operator<<(ostream&,const NISLListNode<Widget>&);
template class NISLList<Widget>;
template ostream& operator<<(ostream&,const NISLList<Widget>&);
