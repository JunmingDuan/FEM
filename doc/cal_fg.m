syms x y lambda mu
u=x^2*y^2;
v=x^2*sin(y);
%f1=-mu*(diff(u,x,2)+diff(u,y,2))-(lambda+mu)*(diff(diff(u,x)+diff(v,y),x));
%f2=-mu*(diff(v,x,2)+diff(v,y,2))-(lambda+mu)*(diff(diff(u,x)+diff(v,y),y));
f1=simplify(-diff(s11,x)-diff(s12,y));
f2=simplify(-diff(s21,x)-diff(s22,y));

e11=diff(u,x);
e22=diff(v,y);
e12=(diff(u,y)+diff(v,x))/2;
e21=e12;

s11=lambda*(e11+e22)+2*mu*e11;
s22=lambda*(e11+e22)+2*mu*e22;
s12=2*mu*e12;
s21=2*mu*e21;

g1l = simplify(-s11);
g2l = simplify(-s12);
g1r = simplify(s11);
g2r = simplify(s12);
g1u = simplify(s12);
g2u = simplify(s22);

%diary fg.dat
%diary on
%fprintf('%s\n',latex(f1));
%fprintf('%s\n',latex(f2));
%fprintf('%s\n',latex(g1l));
%fprintf('%s\n',latex(g2l));
%fprintf('%s\n',latex(g1r));
%fprintf('%s\n',latex(g2r));
%fprintf('%s\n',latex(g1u));
%fprintf('%s\n',latex(g2u));
%diary off

