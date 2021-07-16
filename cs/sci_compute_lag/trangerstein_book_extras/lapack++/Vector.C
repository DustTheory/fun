#include "Vector.H"

template<typename F,typename Z> void Vector<F,Z>::fillWith(Z scalar) {
  for (int i=0; i<sz; i++) data[i] = scalar;
}

template<typename F,typename Z> void Vector<F,Z>::printOn(ostream& s)
const {
  s << "Vector(" << sz << ")\n" ;
  for (int i=0; i<sz; i++) s << data[i] << "  ";
  s << "\n";
}

template<typename F,typename Z> void testVector(const F &abscalar,
const Z &scalar) {
  Vector<F,Z> a1;
  cout << "after Vector(), a1 = " << endl;
  a1.printOn(cout);

  Vector<F,Z> *t;
  {
    int n=3;
    Vector<F,Z> a2(n);
    cout << "\nafter Vector(int), a2 = " << endl;
    a2.printOn(cout);

    {
      Vector<F,Z> a3(n,scalar);
      cout << "\nafter a3=Vector(n,scalar), a3 = " << endl;
      a3.printOn(cout);
    }

    a1.resize(n);
    cout << "\nafter resize, a1 = " << endl;
    a1.printOn(cout);

    a2=static_cast<F>(2.)*scalar;
    cout << "\nafter a2=Z, a2 = " << endl;
    a2.printOn(cout);

    a1=scalar;
    t=a1+a2;
    cout << "\na1+a2 = " << endl;
    t->printOn(cout);
    delete t; t=0;

    t=a1-a2;
    cout << "\na1-a2 = " << endl;
    t->printOn(cout);
    delete t; t=0;

    t=a1*scalar;
    cout << "\na1*scalar = " << endl;
    t->printOn(cout);
    delete t; t=0;

    t=a1/scalar;
    cout << "\na1/scalar = " << endl;
    t->printOn(cout);

    cout << "\na1 = " << endl;
    a1.printOn(cout);
    cout << "a1.amax = " << a1.amax() << endl;
    cout << "a1.asum = " << a1.asum() << endl;
    cout << "a1.nrm2 = " << a1.nrm2() << endl;

    cout << "\na2 = " << endl;
    a2.printOn(cout);
    a1.copy(a2);
    cout << "after a1.copy(a2), a1 = " << endl;
    a1.printOn(cout);

    {
      Vector<F,Z> a3(5);
      a3.copyFrom(n,a2);
      cout << "\nafter a3.copyFrom(a2), a3 = " << endl;
      a3.printOn(cout);
    }

    if (sizeof(F)==sizeof(Z)) {
      cout << "\na1.dot(a2) = " << a1.dot(a2) << endl;
    } else {
      cout << "\na1.dotc(a2) = " << a1.dotc(a2) << endl;
      cout << "a1.dotu(a2) = " << a1.dotu(a2) << endl;
    }
  }

  a1.scal(scalar);
  cout << "\nafter a1.scal, a1 = " << endl;
  a1.printOn(cout);

  a1.swap(*t);
  cout << "\nafter a1.swap(t)" << endl;
  cout << "a1 = " << endl;
  a1.printOn(cout);
  cout << "t = " << endl;
  t->printOn(cout);

  a1.rot(*t,abscalar,abscalar*2.);
  cout << "\nafter a1.rot(t,abscalar,abscalar*2.)" << endl;
  cout << "a1 = " << endl;
  a1.printOn(cout);
  cout << "t = " << endl;
  t->printOn(cout);

  a1.axpy(scalar,*t);
  cout << "\nafter a1.axpy(scalar,t)" << endl;
  cout << "a1 = " << endl;
  a1.printOn(cout);

  {
    Vector<F,F> c(a1.size());
    Vector<F,Z> s(a1.size());
    s.copy(*t);
    a1.largv(s,c);
    cout << "\nafter largv, a1 = " << endl;
    a1.printOn(cout);
    cout << "c = " << endl;
    c.printOn(cout);
    cout << "s = " << endl;
    s.printOn(cout);

    a1.lartv(*t,c,s);
    cout << "\nafter lartv, a1 = " << endl;
    a1.printOn(cout);
    cout << "t = " << endl;
    t->printOn(cout);
  }

  delete t; t=0;
}
