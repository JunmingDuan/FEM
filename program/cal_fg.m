syms x y lambda mu
u = x*y;
v = -x*y;

e11=diff(u,x);
e22=diff(v,y);
e12=(diff(u,y)+diff(v,x))/2;
e21=e12;

s11=lambda*(e11+e22)+2*mu*e11;
s22=lambda*(e11+e22)+2*mu*e22;
s12=2*mu*e12;
s21=2*mu*e21;

f1=simplify(-diff(s11,x)-diff(s12,y))
f2=simplify(-diff(s21,x)-diff(s22,y))

gl1 = simplify(-s11)
gl2 = simplify(-s12)
gr1 = simplify(s11)
gr2 = simplify(s12)
gu1 = simplify(s12)
gu2 = simplify(s22)

