X='find smallest number eps such that 1 + eps > 1:';disp(X);
small=eps('double');
for i=1:2
  lambda=1.+small;
  rho=(lambda-1.)/small;
  pow=round(log2(small));
  X=sprintf('eps = 2**{%d} , ((1.+eps)-1.)/eps = %1.0f\n',pow,rho);disp(X);
  small=0.5*small;
end
X='find largest number eps such that 1 + eps = 1:';disp(X);
small=eps('double');
for i=1:2
  lambda=1.+small;
  rho=(lambda-1.)/small;
  pow=round(log2(small));
  X=sprintf('eps = 2**{%d} , ((1.+eps)-1.)/eps = %1.0f\n',pow,rho);disp(X);
  small=0.5*small;
end
X='find largest number eps such that ( 1 / eps - 1 ) * eps = 1:';disp(X);
small=eps('double');
for i=1:3
  lambda=1./small;
  rho=(lambda-1.)*small;
  pow=round(log2(small));
  X=sprintf('eps = 2**{%d} , (1./eps-1.)*eps = %18.16f\n',pow,rho);disp(X);
  small=0.5*small;
end
X='find largest number eps such that ( 1 - 1 / eps ) + 1 / eps = 0:';disp(X);
small=eps('double');
for i=1:3
  lambda=1./small;
  rho=(1.-lambda)+lambda;
  pow=round(log2(small));
  X=sprintf('eps = 2**{%d} , (1.-1./eps)+1./eps = %1.0f\n',pow,rho);disp(X);
  small=0.5*small;
end
