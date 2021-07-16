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
// "$Header: /home/faculty/johnt/cvs/deal_new/memdebug/NumPtr.C,v 1.1 2009/08/20 17:33:32 johnt Exp $"
#include "NumPtr.H"
#include <algorithm>
#include "Tracer.H"

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> void NumPtr<T>::initialize(T t) const {
  T *d=getData();
  for (unsigned int i=0;i<getNumber();i++) {
    d[i]=t;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> void NumPtr<T>::copyFrom(const NumPtr<T> &p)
const {
  int min_number=min(getNumber(),p.getNumber());
  T *d=getData();
  T *pd=p.getData();
  for (unsigned int i=0;i<min_number;i++) d[i]=pd[i];
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> void NumPtrC<T>::heapify(int index,int heap_size) {
  int l=2*index+1;
  int r=l+1;
  int s=index;
  if (l<heap_size && this->operator[](s)<this->operator[](l)) s=l;
  if (r<heap_size && this->operator[](s)<this->operator[](r)) s=r;
  if (s!=index) {
    swap(s,index);
    heapify(s,heap_size);
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> void NumPtrC<T>::heapSort(){
  int n=this->getNumber();
  if (n<=1) return;
  for (int i=n/2-1;i>=0;i--) {
    heapify(i,n);
  }
  for (int i=n-1;i>=1;i--) {
    swap(0,i);
    heapify(0,i);
  }
#ifdef DEBUG
//printOn(cout);
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> int NumPtrC<T>::binarySearch(T v) const {
  unsigned int l=0,r=this->getNumber()-1;
  while (r>=l) {
    unsigned int i=(l+r)/2;
    T value=this->operator[](i);
    if (v==value) return i;
    if (v<value) r=i-1;
    else l=i+1;
  }
  return -1;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> void NumPtrC<T>::printOn(ostream &os) const {
  os << "NumPtrC: " << endl;
  NumPtr<T>::printOn(os);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> void NumPtrWithValue<T>::heapify(int index,
int heap_size,bool descending) {
  int l=2*index+1;
  int r=l+1;
  int s=index;
  if (descending) {
    if (l<heap_size && value_array[s]>value_array[l]) s=l;
    if (r<heap_size && value_array[s]>value_array[r]) s=r;
  } else {
    if (l<heap_size && value_array[s]<value_array[l]) s=l;
    if (r<heap_size && value_array[s]<value_array[r]) s=r;
  }
  if (s!=index) {
    swap(s,index);
    heapify(s,heap_size,descending);
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> void NumPtrWithValue<T>::heapSort(bool descending){
  if (getNumber()<=1) return;
  for (int i=getNumber()/2-1;i>=0;i--) {
    heapify(i,getNumber(),descending);
  }
  for (int i=getNumber()-1;i>=1;i--) {
    swap(0,i);
    heapify(0,i,descending);
  }
#ifdef DEBUG
//printOn(cout);
  for (int i=0;i<getNumber()-1;i++) {
    if (descending) {
      CHECK_TEST(!(value_array[i]<value_array[i+1]));
    } else {
      CHECK_TEST(!(value_array[i]>value_array[i+1]));
    }
  }
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> int NumPtrWithValue<T>::binarySearch(int v) const {
  int l=0,r=value_array.getNumber()-1;
  while (r>=l) {
    int i=(l+r)/2;
    int value=value_array[i];
    if (v==value) return i;
    if (v<value) r=i-1;
    else l=i+1;
  }
  return -1;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> void NumPtrWithValue<T>::printOn(ostream &os) const {
  os << "NumPtrWithValue: value = " << value << endl;
  for (int i=0;i<getNumber();i++) {
    os << i << "\t: " << NumPtr<T>::operator[](i) << " " 
       << value_array[i] << endl;
  }
  NumPtr<T>::printOn(os);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T,class V> void NumPtrWithCompare<T,V>::heapify(
int index,int heap_size,bool descending) {
  int l=2*index+1;
  int r=l+1;
  int s=index;
  if (descending) {
    if (l<heap_size && *value_array[s] > *value_array[l]) s=l;
    if (r<heap_size && *value_array[s] > *value_array[r]) s=r;
  } else {
    if (l<heap_size && *value_array[s] < *value_array[l]) s=l;
    if (r<heap_size && *value_array[s] < *value_array[r]) s=r;
  }
  if (s!=index) {
    swap(s,index);
    heapify(s,heap_size,descending);
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T,class V> void NumPtrWithCompare<T,V>::heapSort(
bool descending) {
  if (getNumber()<=1) return;
  for (int i=getNumber()/2-1;i>=0;i--) {
    heapify(i,getNumber(),descending);
  }
  for (int i=getNumber()-1;i>=1;i--) {
    swap(0,i);
    heapify(0,i,descending);
  }
#ifdef DEBUG
//printOn(cout);
  for (int i=0;i<getNumber()-1;i++) {
    if (descending) {
      CHECK_TEST(!(*value_array[i] < *value_array[i+1]));
    } else {
      CHECK_TEST(!(*value_array[i] > *value_array[i+1]));
    }
  }
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T,class V> void NumPtrWithCompare<T,V>::printOn(
ostream &os) const{
  os << "NumPtrWithCompare" << endl;
  for (int i=0;i<getNumber();i++) {
    os << i << "\t: " << NumPtr<T>::operator[](i) << " " 
       << value_array[i] << endl;
  }
  NumPtr<T>::printOn(os);
}
