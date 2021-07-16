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
// "$Header: /home/faculty/johnt/cvs/deal_new/gui/VirtualInput.C,v 1.1 2009/08/20 17:32:36 johnt Exp $"
#include "InputParameter.H"
#include <iomanip>
#include <stdlib.h>
#include "Tracer.H"
#include "NISLList.C"
//::::::::::::::::::::::: VirtualInput :::::::::::::::::::::::::::::::::
VirtualInput::VirtualInput(const char *n,const char*) : name(0) { 
  int length_name=strlen(n)+1;
  name=OPERATOR_NEW_BRACKET(char,length_name);
  strcpy(name,n);
}

VirtualInput::~VirtualInput() {
  delete name; name=0;
}

void VirtualInput::printOn(ostream &os) const {
  os << "VirtualInput: name = " << name << endl;
}

ostream& operator<<(ostream &os,const VirtualInput& vi) {
  os << "VirtualInput: name = " << vi.getName() << ", value = "
     << vi.getValue() << endl;
  return os;
}
//:::::::::::::::::::::::: InputString :::::::::::::::::::::::::::::::::
InputString::InputString(char *&p,const char *n,int m,const char*) : 
VirtualInput(n),ptr_string(&p),length_string(0),maxlen(m) {
  if (p) length_string=strlen(p)+1;
  copy=OPERATOR_NEW_BRACKET(char,maxlen); 
  strcpy(copy,"");
}

InputString::~InputString() {
  delete copy; copy=0;
  maxlen=0;
}

ifstream& InputString::formattedRead(ifstream &is) {
  is >> setw(maxlen) >> copy;
  length_string=strlen(copy)+1;
  if (ptr_string[0]) delete ptr_string[0];
  ptr_string[0]=OPERATOR_NEW_BRACKET(char,length_string);
  strcpy(ptr_string[0],copy);
  return is;
}

ofstream& InputString::unformattedWrite(ofstream &os) {
  os.write(reinterpret_cast<const char*>(&length_string),
    sizeof(length_string));
  os.write(ptr_string[0],length_string);
  return os;
}

ifstream& InputString::unformattedRead(ifstream &is) {
  is.read(reinterpret_cast<char*>(&length_string),
    sizeof(length_string));
  is.read(copy,length_string);
  if (ptr_string[0]) delete ptr_string[0];
  ptr_string[0]=OPERATOR_NEW_BRACKET(char,length_string);
  strcpy(ptr_string[0],copy);
  return is;
}

void InputString::printOn(ostream &os) const {
  os << "InputString: ptr_string = " 
     << reinterpret_cast<void*>(ptr_string)
     << "\n\tcopy = " << copy 
     << "\n\tlength_string = " << length_string
     << "\n\tmaxlen = " << maxlen << endl;
  VirtualInput::printOn(os);
}
//::::::::::::::: InputIFStream ::::::::::::::::::::::::::::::::::::::::
InputIFStream::InputIFStream(ifstream &fs,const char *n,int m,
const char*) : VirtualInput(n),stream(&fs),string(0),maxlen(m) {
  string=OPERATOR_NEW_BRACKET(char,maxlen);
  strcpy(string,"");
}

InputIFStream::~InputIFStream() {
  delete string;
  string=0;
  maxlen=0;
}

void InputIFStream::changeString(const char *text) {
  if (strlen(text)<maxlen) {
    stream->open(text,ios::in);
    if (*stream) strcpy(string,text);
    else reOpen();
  }
}

ifstream& InputIFStream::formattedRead(ifstream &is) {
  is >> setw(maxlen) >> string;
  stream->clear(ios::goodbit);
  stream->open(string,ios::in);
  ASSERT(!stream->fail());
  return is;
}

ofstream& InputIFStream::unformattedWrite(ofstream &os) {
  int len=strlen(string)+1;
  os.write(reinterpret_cast<const char*>(&len),sizeof(len));
  os.write(string,len);
  return os;
}

ifstream& InputIFStream::unformattedRead(ifstream &is) {
  int len;
  is.read(reinterpret_cast<char*>(&len),sizeof(len));
  is.read(string,len);
  stream->clear(ios::goodbit);
  stream->open(string,ios::in);
  ASSERT(!stream->fail());
  return is;
}

ifstream& InputIFStream::reOpen() {
  stream->clear(ios::goodbit);
  stream->open(string,ios::in);
  ASSERT(!stream->fail());
  return *stream;
}

void InputIFStream::printOn(ostream &os) const {
  os << "InputIFStream: stream = " << stream
     << "\n\tstring = " << string 
     << "\n\tmaxlen = " << maxlen << endl;
  VirtualInput::printOn(os);
}
//::::::::::::::::::::::::: InputOFStream ::::::::::::::::::::::::::::::
InputOFStream::InputOFStream(ofstream &fs,const char *n,int m,
const char*) : VirtualInput(n),stream(&fs),string(0),maxlen(m) {
  string=OPERATOR_NEW_BRACKET(char,maxlen);
  strcpy(string,"");
}

InputOFStream::~InputOFStream() {
  delete string;
  string=0;
  maxlen=0;
}

void InputOFStream::changeString(const char *text) {
  if (strlen(text)<maxlen) {
    stream->open(text);
    if (*stream) strcpy(string,text);
    else reOpen();
  }
}

ifstream& InputOFStream::formattedRead(ifstream &is) {
  is >> setw(maxlen) >> string;
  stream->clear(ios::goodbit);
#ifdef __DECCXX
  stream->open(string,ios::out);
#else
  stream->open(string,ios::out|ios::binary);
#endif
  ASSERT(!stream->fail());
  return is;
}

ofstream& InputOFStream::unformattedWrite(ofstream &os) {
  int len=strlen(string)+1;
  os.write(reinterpret_cast<const char*>(&len),sizeof(len));
  os.write(string,len);
  return os;
}

ifstream& InputOFStream::unformattedRead(ifstream &is) {
  int len;
  is.read(reinterpret_cast<char*>(&len),sizeof(len));
  is.read(string,len);
  stream->clear(ios::goodbit);
#ifdef __DECCXX
  stream->open(string,ios::out);
#else
  stream->open(string,ios::out|ios::binary);
#endif
  ASSERT(!stream->fail());
  return is;
}

ofstream& InputOFStream::reOpen() {
#ifndef __DECCXX
  if (!stream->is_open()) {
    stream->open(string,ios::out);
  } else 
#endif
  {
    stream->seekp(0,ios::beg);
    stream->clear(ios::goodbit);
  }
  ASSERT(!stream->fail());
  return *stream;
}

void InputOFStream::printOn(ostream &os) const {
  os << "InputOFStream: stream = " << stream
     << "\n\tstring = " << string 
     << "\n\tmaxlen = " << maxlen << endl;
  VirtualInput::printOn(os);
}
//:::::::::::::: InputParameterList ::::::::::::::::::::::::::::::::::::
InputParameterList::InputParameterList(const char *n) :
NISLList<VirtualInput>(),name(0) {
  int length_name=strlen(n)+1;
  name=OPERATOR_NEW_BRACKET(char,length_name);
  strcpy(name,n);
}

InputParameterList::~InputParameterList() {
  while (notEmpty()) delete delAfter(0);
  delete name; name=0;
}

ifstream& InputParameterList::formattedRead(ifstream &is,const char *n){
  for (NISLListNode<VirtualInput> *p=first();p;p=next(p)) {
    if (p->getData()->isNamed(n)) { 
      p->getData()->formattedRead(is); 
      break; 
    } 
  }
  return is;
}

ofstream& InputParameterList::unformattedWrite(ofstream &os) {
  for (NISLListNode<VirtualInput> *p=first();p;p=next(p)) {
    p->getData()->unformattedWrite(os); 
  }
  return os;
}

ifstream& InputParameterList::unformattedRead(ifstream &is) {
  for (NISLListNode<VirtualInput> *p=first();p;p=next(p)) {
    p->getData()->unformattedRead(is); 
  }
  return is;
}

ostream& operator<<(ostream &s,const InputParameterList &l) {
   s << *static_cast<const NISLList<VirtualInput> *>(&l);
   return s;
}

void InputParameterList::printOn(ostream &os) const {
  os << "InputParameterList:" << endl;
  for (NISLListNode<VirtualInput> *p=first();p;p=next(p)) {
    os << "\n\tp = " << p << endl;
    os << "\tp->getData = " << p->getData() << endl;
    p->getData()->printOn(os);
  }
//NISLList<VirtualInput>::printOn(os);
}

template class NISLListNode<VirtualInput>;
template ostream& operator<<(ostream&,const NISLListNode<VirtualInput>&);
template class NISLList<VirtualInput>;
template ostream& operator<<(ostream&,const NISLList<VirtualInput>&);
