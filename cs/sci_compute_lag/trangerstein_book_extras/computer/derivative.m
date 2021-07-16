exact=cos(1.);
h=1.;  
logh=zeros(64);
logerror=zeros(64);
for n=1:64
  diff=(sin(1.+h)-sin(1.))/h;
  logh(n)=-log10(h);
  logerror(n)=log10(abs(diff-exact));
  h=h*0.5; 
end
plot(logh,logerror);
