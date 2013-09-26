function [ferr,Jerr] = testnine
z=rand(1,8); p=rand(1,2);
eps=1e-7;
f=ninepointf(z,p);
df=ninepointdf(z,p);
E=eye(8);
J=zeros(1,8);
for ii=1:8
    zz=z+eps*E(ii,:);
    ff=ninepointf(zz,p);
    J(ii) = (ff-f)/eps;
end;
Jerr = df-J;
C=rand(1,1)+1i*rand(1,1);
D=rand(1,1)+1i*rand(1,1);
P=rand(1,1)+1i*rand(1,1);
t=rand(1,1)+1i*rand(1,1);
th=exp(1i*rand(1,1));
CC=t+th*C;
DD=t+th*D;
PP=t+th*P;
A=(DD+D)/2+1i*rand(1,1)*(DD-D);
B=(CC+C)/2+1i*rand(1,1)*(CC-C);
x=D-P; y=C-P; a=A-P; b=B-P;
z=[x x' y y' a a' b b'];
p=[PP-P (PP-P)'];
ferr=ninepointf(z,p);


function [f]=ninepointf(z,p)
% function value
% ordering in z is [ x xx y yy a aa b bb ]
% p is 2-vector [ \delta, \bar\delta ] = precision point
x=z(1); xx=z(2); y=z(3); yy=z(4);
a=z(5); aa=z(6); b=z(7); bb=z(8);
d=p(1); dd=p(2);
c1 = (aa-dd)*x;
c2 = ( a- d)*xx;
c3 = d*(aa-xx) + dd*( a- x) - d*dd;
d1 = (bb-dd)*y;
d2 = ( b- d)*yy;
d3 = d*(bb-yy) + dd*( b- y) - d*dd;
g1 = c2*d3-d2*c3;
g2 = c3*d1-d3*c1;
g3 = c1*d2-d1*c2;
f = g1*g2+g2*g3+g3*g1;

function [df]=ninepointdf(z,p)
% Jacobian
% ordering in z is [ x xx y yy a aa b bb ]
% p is 2-vector [ \delta, \bar\delta ] = precision point
x=z(1); xx=z(2); y=z(3); yy=z(4);
a=z(5); aa=z(6); b=z(7); bb=z(8);
d=p(1); dd=p(2);
c1 = (aa-dd)*x;
c2 = ( a- d)*xx;
c3 = d*(aa-xx) + dd*( a- x) - d*dd;
d1 = (bb-dd)*y;
d2 = ( b- d)*yy;
d3 = d*(bb-yy) + dd*( b- y) - d*dd;
g1 = c2*d3-d2*c3;
g2 = c3*d1-d3*c1;
g3 = c1*d2-d1*c2;
%f = g1*g2+g2*g3+g3*g1;
dg1 = [ dd*d2             (a-d)*d3+d*d2 -dd*c2           
        -d*c2-(b-d)*c3 xx*d3-dd*d2 -d*d2     dd*c2-yy*c3 d*c2];
dg2 = [ -dd*d1-(aa-dd)*d3 -d*d1         (bb-dd)*c3+dd*c1 d*c1          
         dd*d1       d*d1-x*d3 -dd*c1      y*c3-d*c1]; 
dg3 = [ (aa-dd)*d2        -(a-d)*d1     -(bb-dd)*c2      (b-d)*c1   
         -xx*d1      x*d2      yy*c1       -y*c2];
df = dg1*(g2+g3) + dg2*(g1+g3) + dg3*(g2+g1);
