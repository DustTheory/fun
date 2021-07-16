X='compute f_n = 0.5*f_{n-1} + 1., f_1 = 1 until rounding error'; disp(X);
format hex;
fold=single(1.);
fnew=single(1.5); 
while (fnew<2.)
  X=num2str(fold,24);disp(X);disp(fold);fold=fnew; fnew=1.+0.5*fnew; 
end; 
X=num2str(fnew,24);disp(X);disp(fnew);

format short;
X='compute f_n = 2.*f_{n-1} + 1., f_1 = 1 until overflow'; disp(X);
format hex;
f=single(1.);
while isfinite(f)
  X=num2str(f,24);disp(X);disp(f);f=2.*f+1.;
end;
X=num2str(f,24);disp(X);disp(f);f=2.*f+1.;

format short;
X='compute f_n = 0.5*f_{n-1}, f_1 = 2-epsilon until underflow'; disp(X);
format hex;
while (fold>0.)
  X=num2str(fold,24);disp(X);disp(fold);fold=0.5*fold;
end;
X=num2str(fold,24);disp(X);disp(fold)
