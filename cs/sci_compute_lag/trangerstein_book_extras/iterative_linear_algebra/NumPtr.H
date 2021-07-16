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
#ifndef NUMPTR_H
#define NUMPTR_H "$Header: /home/faculty/johnt/cvs/deal_new/memdebug/NumPtr.H,v 1.1 2009/08/20 17:33:32 johnt Exp $"
#include <fstream>
#include <iostream>
#include <limits.h>
#include "Errors.H"
#include "MemoryDebugger.H"

template<class T,class N> class NumPtrBase {
  protected:
    T *ptr_data;
    N number;

    inline void setData(T *p) { ptr_data=p; }
    inline void copyDataFrom(const NumPtrBase &p) { 
      ptr_data=p.ptr_data; 
    }
    inline void zeroData() { ptr_data=0; }
    inline void setNumber(N n) { number=n; }
    inline void copyNumberFrom(const NumPtrBase &p) { number=p.number; }
    inline void decrementNumber() { number--; }
    inline void decrementNumberBy(N i) { number-=i; }
    inline void incrementNumber() { number++; }
    inline void incrementNumberBy(N i) { number+=i; }
  public:
    inline NumPtrBase() : ptr_data(0),number(0) {;}
    inline NumPtrBase(N s,T *p) : ptr_data(p),number(s) {
      CHECK_TEST((s==0) || (p!=0))
    }
    inline NumPtrBase(const NumPtrBase &p) : ptr_data(p.ptr_data),
      number(p.number) {;}
    inline virtual ~NumPtrBase() {;}

    inline void cleanup() { ptr_data=0; number=0; }
//  inline T operator*() const { 
//    CHECK_TEST(ptr_data);
//    return *ptr_data; 
//  }
//  inline T& operator*() { 
//    CHECK_TEST(ptr_data);
//    return *ptr_data; 
//  }
    inline NumPtrBase& operator=(const NumPtrBase &p) {
      CHECK_TEST(number==0)
      CHECK_TEST(ptr_data==0)
      ptr_data=p.ptr_data; number=p.number;
      return *this;
    }
    inline bool operator!() const { return !ptr_data; }
    inline bool operator<(T *p) const { return ptr_data < p; }
    inline bool operator<(const NumPtrBase &p) const { 
      return ptr_data < p.ptr_data;
    }
    inline bool operator>(T *p) const { return ptr_data > p; }
    inline bool operator>(const NumPtrBase &p) const { 
      return ptr_data > p.ptr_data;
    }
    inline bool operator<=(T *p) const { return ptr_data <= p; }
    inline bool operator<=(const NumPtrBase &p) const { 
      return ptr_data <= p.ptr_data;
    }
    inline bool operator>=(T *p) const { return ptr_data >= p; }
    inline bool operator>=(const NumPtrBase &p) const { 
      return ptr_data >= p.ptr_data;
    }
    inline bool operator==(const NumPtrBase &p) const { 
      return ptr_data==p.ptr_data 
      ;
    } 
    inline bool operator!=(const NumPtrBase &p) const { 
      return ptr_data!=p.ptr_data;
    } 

    inline T* getData() const { return ptr_data; }
    inline N getNumber() const { return number; }
    inline void undefine() {
      ptr_data=static_cast<T*>(0); 
      number=0;
    }

    virtual void printOn(ostream &os = cerr) const {
      os << "NumPtrBase: this = " << this
         << ", ptr_data = " << reinterpret_cast<void*>(ptr_data)
         << ", number = " << number
         << endl;
    }
};

template<class T> class NumPtr : public NumPtrBase<T,int> {
  private:
    bool owns_data;
  protected:
    inline void deleteData() { delete [] NumPtrBase<T,int>::ptr_data; }
    inline void decrementData() { NumPtrBase<T,int>::ptr_data--; }
    inline void decrementDataBy(int i) {NumPtrBase<T,int>::ptr_data-=i;}
    inline void incrementData() { NumPtrBase<T,int>::ptr_data++; }
    inline void incrementDataBy(int i) {NumPtrBase<T,int>::ptr_data+=i;}

    inline void zeroData() { NumPtrBase<T,int>::zeroData(); }
    inline void setNumber(int n) { NumPtrBase<T,int>::setNumber(n); }
    inline void decrementNumber() { 
      NumPtrBase<T,int>::decrementNumber(); 
    }
    inline void decrementNumberBy(int i) { 
      NumPtrBase<T,int>::decrementNumberBy(i);
    }
    inline void incrementNumber() { 
      NumPtrBase<T,int>::incrementNumber();
    }
    inline void incrementNumberBy(int i) { 
      NumPtrBase<T,int>::incrementNumberBy(i);
    }
  public:
    inline NumPtr() : owns_data(false) {;}
    inline explicit NumPtr(unsigned int s) : owns_data(false) {
      if (s>0) { 
        setData(OPERATOR_NEW_BRACKET(T,s)); 
        setNumber(s);
        owns_data=true; 
      } 
    }
    inline NumPtr(int s,T *p) : NumPtrBase<T,int>(s,p),owns_data(false)
      {}
    inline NumPtr(const NumPtr &p) : NumPtrBase<T,int>(p),
      owns_data(false) {;}
    inline NumPtr(const NumPtr &p,int s) : owns_data(false) {
      CHECK_TEST(s<=p.getNumber())
      setData(p.getData()+s);
      setNumber(p.getNumber()-s);
    }
    inline NumPtr(const NumPtr &p,unsigned int s,unsigned int sz) : owns_data(false) {
      CHECK_TEST(s+sz<=p.getNumber())
      setData(p.getData()+s);
      setNumber(sz);
    }
    inline bool ownsData() const { return owns_data; }
    inline T* getData() const { return NumPtrBase<T,int>::getData(); }
    inline int getNumber() const { 
      return NumPtrBase<T,int>::getNumber(); 
    }
    inline void cleanup() {
      if (owns_data && getData()!=0) {
        deleteData();
      }
      NumPtrBase<T,int>::cleanup();
      owns_data=false;
    }
//  calling NumPtr::~NumPtr can mess up the virtual function table
    virtual ~NumPtr() { cleanup(); }
//  the following is used by PatchArr::makePatchesFrom
    void copy(NumPtr &c) {
      CHECK_POSITIVE(c.getNumber())
      if (owns_data && getNumber()>0) {
#ifdef DEBUG
        if (getData()==0) {
          cerr << "assertion getData()!=0 failed in file " << __FILE__
               << " at line " << (__LINE__-2) << endl;
          abort();
        }
#endif
        deleteData();
      }
      copyNumberFrom(c);
      copyDataFrom(c);
      owns_data=c.owns_data; c.owns_data=false;
    }
    inline NumPtr& operator=(const NumPtr &p) {
      NumPtrBase<T,int>::operator=(p);
      owns_data=false;
      return *this;
    }
    inline T& operator[](int i) { 
      CHECK_BOUNDS(i,0,getNumber())
      return *(getData()+i);
    }
    inline const T& operator[](int i) const {
      CHECK_BOUNDS(i,0,getNumber())
      return *(getData()+i);
    }
    inline const NumPtr<T> operator+(int i) const {
      return NumPtr<T>(getNumber()-i,getData()+i);
    }
    inline NumPtr<T>& operator++() { // prefix: ++p
      CHECK_TEST(!owns_data)
      incrementData(); decrementNumber();
      CHECK_TEST(getNumber()>=0)
      if (getNumber()==0) zeroData();
      return *this;
    }
//  inline NumPtr<T>& operator++(int) { // postfix: p++
//    CHECK_TEST(!owns_data)
//    NumPtr<T> save(*this);
//    incrementData(); decrementNumber();
//    CHECK_TEST(getNumber()>=0)
//    if (getNumber()==0) zeroData();
//    return save;
//  }
    inline NumPtr<T>& operator+=(int i) {
      CHECK_TEST(!owns_data)
      incrementDataBy(i); decrementNumberBy(i);
      CHECK_TEST(getNumber()>=0)
      if (getNumber()==0) zeroData();
      return *this;
    }
    inline bool contains(const NumPtr<T> &t) const {
      return t>=getData() && t<getData()+getNumber();
    }
    inline bool contains(const T *t) const {
      return t>=getData() && t<getData()+getNumber();
    }
    inline void allocate(int n) {
      CHECK_TEST(!owns_data)
      setNumber(n);
      if (n>0) {
        setData(OPERATOR_NEW_BRACKET(T,n));
        owns_data=true;
      } else zeroData();
    }
    void allocate(int n,T *t) {
      CHECK_TEST(!owns_data)
      setNumber(n);
      if (n>0) setData(t);
      else zeroData();
      owns_data=false;
    }
    void undefine() {
      CHECK_TEST(!owns_data)
      NumPtrBase<T,int>::undefine();
    }
    void initialize(T) const ;
    virtual void copyFrom(const NumPtr&) const ;
    virtual void write(ofstream &os) const {
      CHECK_NONNEGATIVE(getNumber())
      CHECK_TEST((getNumber()==0) == (getData()==0))
      int n=getNumber();
      os.write(reinterpret_cast<const char*>(&n),sizeof(int));
      if (getNumber()>0) {
        os.write(reinterpret_cast<char*>(getData()),
                 getNumber()*sizeof(T));
      }
    }
    virtual void read(ifstream &is) {
      int len;
      is.read(reinterpret_cast<char*>(&len),sizeof(len));
      CHECK_SAME(len,getNumber())
      if (getNumber()>0) {
        is.read(reinterpret_cast<char*>(getData()),
                getNumber()*sizeof(T));
      }
    }

    virtual void printOn(ostream &os = cerr) const {
      os << "NumPtr: owns_data = " << owns_data << endl;
      NumPtrBase<T,int>::printOn(os);
    }
};
#define INSTANTIATE_NUMPTR(c) template class NumPtrBase<c,int>; template class NumPtr<c>;

//bool operator<(const T&,const T&) and
//bool operator==(const T&,const T&) must be defined
template<class T> class NumPtrC : public NumPtr<T> {
  private:
    void swap(int i,int j) {
      T tmp=this->operator[](i);
      this->operator[](i)=this->operator[](j);
      this->operator[](j)=tmp;
    }
    void heapify(int,int);
  public:
    inline NumPtrC() {;}
    inline NumPtrC(unsigned int s) : NumPtr<T>(s) {;}
    inline NumPtrC(int s,T *p) : NumPtr<T>(s,p) {;}
    inline NumPtrC(const NumPtrC &p) : NumPtr<T>(p) {;}
    inline NumPtrC(const NumPtrC &p,int s) : NumPtr<T>(p,s) {;}
    inline NumPtrC(const NumPtrC &p,unsigned int s,unsigned int sz) : NumPtr<T>(p,s,sz)
      {;}

//  heapSort puts array in increasing order only
    void heapSort();
    int binarySearch(T) const;
    virtual void printOn(ostream &os = cerr) const;
};
#define INSTANTIATE_NUMPTRC(c) template class NumPtrBase<c,int>; template class NumPtr<c>; template class NumPtrC<c>;

template<class T> class NumPtrWithValue : public NumPtr<T> {
  private:
    NumPtr<int> value_array;
    int value;

//  T& operator[](int i) { return NumPtr<T>::operator[](i); }
//  const T& operator[](int i) const { return NumPtr<T>::operator[](i);}
    void swap(int i,int j) {
      T tmp=NumPtr<T>::operator[](i);
      int tmpv=value_array[i];
      NumPtr<T>::operator[](i)=NumPtr<T>::operator[](j);
      value_array[i]=value_array[j];
      NumPtr<T>::operator[](j)=tmp;
      value_array[j]=tmpv;
    }
    virtual void copyFrom(const NumPtr<T> &p) const {
      OBSOLETE("need NumPtrWithValue arg")
    }
  protected:
    inline void setNumber(int n) { NumPtr<T>::setNumber(n); }
    inline void decrementNumber() { NumPtr<T>::decrementNumber(); }
    inline void decrementNumberBy(int i) { 
      NumPtr<T>::decrementNumberBy(i);
    }
    inline void incrementNumber() { NumPtr<T>::incrementNumber(); }
    inline void incrementNumberBy(int i) { 
      NumPtr<T>::incrementNumberBy(i);
    }
    void heapify(int,int,bool);
  public:
    NumPtrWithValue() : value(0) {;}
    NumPtrWithValue(unsigned int s) : NumPtr<T>(s),value_array(s),value(0)
    {
      value_array.initialize(0);
    }
    NumPtrWithValue(int s,T *p,int *v) : NumPtr<T>(s,p),
    value_array(s,v),value(0) {
      for (int i=0;i<s;i++) value += value_array[i];
    }
    NumPtrWithValue(const NumPtrWithValue &p) : NumPtr<T>(p),
      value_array(p.value_array),value(p.value) {;}
    NumPtrWithValue(const NumPtrWithValue &p,int s) : NumPtr<T>(p,s),
    value_array(p.value_array,s),value(0) {
      for (int i=0;i<getNumber();i++) value += value_array[i];
    }
    NumPtrWithValue(const NumPtrWithValue &p,unsigned int s,unsigned int sz) : 
    NumPtr<T>(p,s,sz),value_array(p.value_array,s,sz),value(0) {
      for (int i=0;i<getNumber();i++) value += value_array[i];
    }
    inline int getNumber() const { return NumPtr<T>::getNumber(); }
    inline T* getData() const { return NumPtr<T>::getData(); }
    void cleanup() {
      NumPtr<T>::cleanup();
      value_array.cleanup();
      value=0;
    }
//  calling NumPtr::~NumPtr can mess up the virtual function table
    virtual ~NumPtrWithValue() { cleanup(); }
//  the following is used by PatchArr::makePatchesFrom
    void copy(NumPtrWithValue &c) {
      NumPtr<T>::copy(c);
      value_array.copy(c.value_array);
      value=c.value;
    }
    NumPtrWithValue& operator=(const NumPtrWithValue &p) {
      NumPtr<T>::operator=(p);
      value_array=p.value_array;
      value=p.value;
      return *this;
    }
    int getTotalValue() const { return value; } 
    int getValue(int i) const { return value_array[i]; }
    void getValue(int i,T &t,int &v) const { 
      t=NumPtr<T>::operator[](i); v=value_array[i];
    }
    void setValue(int i,const int &v) { 
      value += v-value_array[i];
      value_array[i]=v;
    }
    void setValue(int i,const T &t,const int &v) { 
      value += v-value_array[i];
      NumPtr<T>::operator[](i)=t; value_array[i]=v;
    }
    const NumPtrWithValue<T> operator+(int i) const {
      return NumPtrWithValue<T>(getNumber()-i,getData()+i,
        value_array.getData()+i);
    }
    NumPtrWithValue<T>& operator++() {
      value-=value_array[0];
      NumPtr<T>::operator++();
      value_array.operator++();
      return *this;
    }
    NumPtrWithValue<T>& operator+=(int i) {
      for (int j=0;j<i;j++) value-=value_array[j];
      NumPtr<T>::operator+=(i);
      value_array+=i;
      return *this;
    }
    bool contains(const NumPtrWithValue<T> &t) const {
      return NumPtr<T>::contains(t) 
          && value_array.contains(t.value_array);
    }
    bool contains(const T *t) const {
      return NumPtr<T>::contains(t);
    }
    void allocate(int n) {
      NumPtr<T>::allocate(n);
      value_array.allocate(n);
      value=0;
    }
    void allocate(int n,T *t,int *v) {
      NumPtr<T>::allocate(n,t);
      value_array.allocate(n,v);
      value=0;
    }
    void undefine() {
      NumPtr<T>::undefine();
      value_array.undefine();
      value=INT_MAX;
    }
    void initialize(T t,int v) {
      NumPtr<T>::initialize(t);
      value_array.initialize(v);
      value=getNumber()*v;
    }
    virtual void copyFrom(const NumPtrWithValue &p) {
      NumPtr<T>::copyFrom(p);
      value_array.copyFrom(p.value_array);
      value=0;
      for (int i=0;i<getNumber();i++) value += value_array[i];
    }
    void heapSort(bool=true);
    int binarySearch(int) const;
    virtual void write(ofstream &os) const {
      NumPtr<T>::write(os);
      value_array.write(os);
      os.write(reinterpret_cast<const char*>(&value),sizeof(value));
    }
    virtual void read(ifstream &is) {
      NumPtr<T>::read(is);
      value_array.read(is);
      is.read(reinterpret_cast<char*>(&value),sizeof(value));
    }

    virtual void printOn(ostream &os = cerr) const;
};

//work with pointer type V to avoid copying *V in swaps
//class *V must have operators <, >
template<class T,class V> class NumPtrWithCompare : public NumPtr<T> {
  private:
    NumPtr<V> value_array;

    T& operator[](int i) { return NumPtr<T>::operator[](i); }
    const T& operator[](int i) const { return NumPtr<T>::operator[](i);}
    void swap(int i,int j) {
      T tmp=NumPtr<T>::operator[](i);
      V tmpv=value_array[i];
      NumPtr<T>::operator[](i)=NumPtr<T>::operator[](j);
      value_array[i]=value_array[j];
      NumPtr<T>::operator[](j)=tmp;
      value_array[j]=tmpv;
    }
    virtual void copyFrom(const NumPtr<T> &p) const {
      OBSOLETE("need NumPtrWithCompare arg")
    }
  protected:
    inline void setNumber(int n) { NumPtr<T>::setNumber(n); }
    inline void decrementNumber() { NumPtr<T>::decrementNumber(); }
    inline void decrementNumberBy(int i) { 
      NumPtr<T>::decrementNumberBy(i);
    }
    inline void incrementNumber() { NumPtr<T>::incrementNumber(); }
    inline void incrementNumberBy(int i) { 
      NumPtr<T>::incrementNumberBy(i);
    }
    void heapify(int,int,bool);
  public:
    NumPtrWithCompare() {;}
    NumPtrWithCompare(unsigned int s) : NumPtr<T>(s),value_array(s) {
      value_array.initialize(0);
    }
    NumPtrWithCompare(int s,T *p,V *v) : NumPtr<T>(s,p),
      value_array(s,v) { }
    NumPtrWithCompare(const NumPtrWithCompare &p) : NumPtr<T>(p),
      value_array(p.value_array) {;}
    NumPtrWithCompare(const NumPtrWithCompare &p,int s) : 
      NumPtr<T>(p,s),value_array(p.value_array,s) { }
    NumPtrWithCompare(const NumPtrWithCompare &p,unsigned int s,unsigned int sz) : 
      NumPtr<T>(p,s,sz),value_array(p.value_array,s,sz) { }
    void cleanup() {
      NumPtr<T>::cleanup();
      value_array.cleanup();
    }
//  calling NumPtr::~NumPtr can mess up the virtual function table
    virtual ~NumPtrWithCompare() { cleanup(); }
    inline int getNumber() const { return NumPtr<T>::getNumber(); }
    inline T* getData() const { return NumPtr<T>::getData(); }
    void copy(NumPtrWithCompare &c) {
      NumPtr<T>::copy(c);
      value_array.copy(c.value_array);
    }
    NumPtrWithCompare& operator=(const NumPtrWithCompare &p) {
      NumPtr<T>::operator=(p);
      value_array=p.value_array;
      return *this;
    }
    void getValue(int i,T &t,V &v) const { 
      t=NumPtr<T>::operator[](i); v=value_array[i];
    }
    void setValue(int i,const T &t,const V &v) { 
      NumPtr<T>::operator[](i)=t; value_array[i]=v;
    }
    const NumPtrWithCompare operator+(int i) const {
      return NumPtrWithCompare(getNumber()-i,getData()+i,
        value_array.getData()+i);
    }
    NumPtrWithCompare& operator++() {
      NumPtr<T>::operator++();
      value_array.operator++();
      return *this;
    }
    NumPtrWithCompare& operator+=(int i) {
      NumPtr<T>::operator+=(i);
      value_array+=i;
      return *this;
    }
    bool contains(const NumPtrWithCompare &t) const {
      return NumPtr<T>::contains(t) 
          && value_array.contains(t.value_array);
    }
    bool contains(const T *t) const {
      return NumPtr<T>::contains(t);
    }
    void allocate(int n) {
      NumPtr<T>::allocate(n);
      value_array.allocate(n);
    }
    void allocate(int n,T *t,V *v) {
      NumPtr<T>::allocate(n,t);
      value_array.allocate(n,v);
    }
    void undefine() {
      NumPtr<T>::undefine();
      value_array.undefine();
    }
    void initialize(T t,V v) {
      NumPtr<T>::initialize(t);
      value_array.initialize(v);
    }
    virtual void copyFrom(const NumPtrWithCompare &p) {
      NumPtr<T>::copyFrom(p);
      value_array.copyFrom(p.value_array);
    }
    void heapSort(bool=true);
    virtual void write(ofstream &os) const {
      NumPtr<T>::write(os);
      value_array.write(os);
    }
    virtual void read(ifstream &is) {
      NumPtr<T>::read(is);
      value_array.read(is);
    }

    virtual void printOn(ostream &os = cerr) const;
};
#endif
