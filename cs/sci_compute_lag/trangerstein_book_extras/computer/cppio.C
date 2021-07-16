#include <cstring>
#include <fstream>
#include <iostream>
using namespace std;

int main(int /*argc*/,char** /*argv*/) {
  cout << boolalpha;
  const char *str="abcdefghijklmnopqrstuvwxyz";
//formatted to stdout
  cout << "in cppio" << endl;
  for (char c='A';c<='Z';c++) cout << c;
  cout << endl;
  cout << str << endl;
  bool b=true;
  cout << "bool = " << b << endl;
  short s=1;
  cout << "short = " << s << endl;
  int i=2;
  cout << "int = " << i << endl;
  long l=3;
  cout << "long = " << l << endl;
  float f=4.;
  cout << "float = " << f << endl;
  double d=5.;
  cout << "double = " << d << endl;

//formated to file
  cout << "\nformatted write to file" << endl;
  ofstream outs;
  outs.open("cppio_formatted_output");
  for (char c='A';c<='Z';c++) outs << c;
  outs << endl;
  outs << str << endl;
  outs << "bool = " << b << endl;
  outs << "short = " << s << endl;
  outs << "int = " << i << endl;
  outs << "long = " << l << endl;
  outs << "float = " << f << endl;
  outs << "double = " << d << endl;
  outs.close();

//formated from file
  cout << "\nformatted read from file" << endl;
  ifstream ins;
  ins.open("cppio_formatted_output");
  for (char c='A';c<='Z';c++) {
    ins >> c;
    cout << c;
  }
  char *str2=new char[128];
  char *str3=new char[128];
  ins.getline(str3,128);
  cout << endl;
  ins >> str2;
  cout << str2 << endl;
  ins.getline(str3,128);
  char c2;
  bool b2;
  ins >> str3 >> c2 >> b2;
  cout << "bool = " << b2 << endl;
  ins.getline(str3,128);
  short s2;
  ins >> str3 >> c2 >> s2;
  cout << "short = " << s2 << endl;
  ins.getline(str3,128);
  int i2;
  ins >> str3 >> c2 >> i2;
  cout << "int = " << i2 << endl;
  ins.getline(str3,128);
  long l2;
  ins >> str3 >> c2 >> l2;
  cout << "long = " << l2 << endl;
  ins.getline(str3,128);
  float f2;
  ins >> str3 >> c2 >> f2;
  cout << "float = " << f2 << endl;
  ins.getline(str3,128);
  double d2;
  ins >> str3 >> c2 >> d2;
  cout << "double = " << d2 << endl;
  ins.getline(str3,128);
  ins.close();

//unformated to file
  cout << "\nunformatted write to file" << endl;
  outs.open("cppio_unformatted_output");
  for (char c='A';c<='Z';c++) outs.put(c);
  outs.put('\n');
//cout << "after upcase, tellp = " << outs.tellp() << endl;
  int n=strlen(str);
  outs.write(reinterpret_cast<const char*>(&n),sizeof(n));
  outs.write(str,n);
  outs.write(reinterpret_cast<const char*>(&b),sizeof(b));
  outs.write(reinterpret_cast<const char*>(&s),sizeof(s));
  outs.write(reinterpret_cast<const char*>(&i),sizeof(i));
  outs.write(reinterpret_cast<const char*>(&l),sizeof(l));
  outs.write(reinterpret_cast<const char*>(&f),sizeof(f));
  outs.write(reinterpret_cast<const char*>(&d),sizeof(d));
  outs.close();

//unformated from file
  cout << "\nunformatted read from file" << endl;
  ins.open("cppio_unformatted_output");
  while (true) {
    char cc;
    ins.read(&cc,1);
    cout << cc;
    if (cc=='\n') break;
  }
//cout << "after upcase, tellg = " << ins.tellg() << endl;
  int nn;
  ins.read(reinterpret_cast<char*>(&nn),sizeof(int));
  char *strstr=new char[nn+1];
  ins.read(strstr,nn);
  strstr[nn]='0';
  cout << strstr << endl;
  bool bb;
  ins.read(reinterpret_cast<char*>(&bb),sizeof(bb));
  cout << "bool = " << bb << endl;
  short ss;
  ins.read(reinterpret_cast<char*>(&ss),sizeof(ss));
  cout << "short = " << ss << endl;
  int ii;
  ins.read(reinterpret_cast<char*>(&ii),sizeof(ii));
  cout << "int = " << ii << endl;
  long ll;
  ins.read(reinterpret_cast<char*>(&ll),sizeof(ll));
  cout << "long = " << ll << endl;
  float ff;
  ins.read(reinterpret_cast<char*>(&ff),sizeof(ff));
  cout << "float = " << ff << endl;
  double dd;
  ins.read(reinterpret_cast<char*>(&dd),sizeof(dd));
  cout << "double = " << dd << endl;
  ins.close();
}
