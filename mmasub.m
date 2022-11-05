function [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2, ...
f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a,c,d);
epsimin = sqrt(m+n)*10^(-9);
feps = 0.000001;
asyinit = 0.5;  
asyincr = 1.2;
asydecr = 0.7;  
albefa = 0.1;
een = ones(n,1);
zeron = zeros(n,1);
asyinit = 2.e-3;            %step size

% Calculation of the asymptotes low and upp :
if iter < 2.5
    low = xval - asyinit*(xmax-xmin);
    upp = xval + asyinit*(xmax-xmin);
else
    zzz = (xval-xold1).*(xold1-xold2);
    factor = een;
    factor(find(zzz > 0)) = asyincr;
    factor(find(zzz < 0)) = asydecr;
    low = xval - factor.*(xold1 - low);
    upp = xval + factor.*(upp - xold1);
end

% Calculation of the bounds alfa and beta :
zzz = low + albefa*(xval-low);
alfa = max(zzz,xmin);
zzz = upp - albefa*(upp-xval);
beta = min(zzz,xmax);

% Calculations of p0, q0, P, Q and b.

ux1 = upp-xval;         ux2 = ux1.*ux1;           ux3 = ux2.*ux1;
xl1 = xval-low;         xl2 = xl1.*xl1;           xl3 = xl2.*xl1;
ul1 = upp-low;          ulinv1 = een./ul1;       uxinv1 = een./ux1;
xlinv1 = een./xl1;      uxinv3 = een./ux3;      xlinv3 = een./xl3;
diap = (ux3.*xl1)./(2*ul1);         diaq = (ux1.*xl3)./(2*ul1);

p0 = zeron;
p0(find(df0dx > 0)) = df0dx(find(df0dx > 0));
p0 = p0 + 0.001*abs(df0dx) + feps*ulinv1;
p0 = p0.*ux2;

q0 = zeron;
q0(find(df0dx < 0)) = -df0dx(find(df0dx < 0));
q0 = q0 + 0.001*abs(df0dx) + feps*ulinv1;
q0 = q0.*xl2;

dg0dx2 = 2*(p0./ux3 + q0./xl3);
del0 = df0dx2 - dg0dx2;
delpos0 = zeron;
delpos0(find(del0 > 0)) = del0(find(del0 > 0));

p0 = p0 + delpos0.*diap;
q0 = q0 + delpos0.*diaq;

P = zeros(m,n);
P(find(dfdx > 0)) = dfdx(find(dfdx > 0));
P = P * diag(ux2);

Q = zeros(m,n);
Q(find(dfdx < 0)) = -dfdx(find(dfdx < 0));
Q = Q * diag(xl2);

dgdx2 = 2*(P*diag(uxinv3) + Q*diag(xlinv3));

del = dfdx2 - dgdx2;
delpos = zeros(m,n);
delpos(find(del > 0)) = del(find(del > 0));

P = P + delpos*diag(diap);
Q = Q + delpos*diag(diaq);
b = P*uxinv1 + Q*xlinv1 - fval ;

%%% Solving the subproblem by a primal-dual Newton method
[xmma,ymma,zmma,lam,xsi,eta,mu,zet,s] = ...
subsolv(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d);