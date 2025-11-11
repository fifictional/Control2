syms theta(t) alpha(t) u t 
dt = diff(theta);
d2t = diff(theta, 2);
da = diff(alpha);
d2a = diff(alpha,2);

ts = 0.002;

Rm = 7.5;
kt = 0.042;
km = 0.042;

mr = 0.095;
r = 0.085;
Jr = (mr*r^2)/3;
br = 1e-3;

mp = 0.024;
Lp = 0.129;
I = Lp/2;
Jp = (mp*Lp^2)/3;
bp = 5e-5;

g = 9.81;


d = Jp*(Jr+mp*r^2)-(mp*I*r)^2;
im = (u-km*dt)/Rm;

% eq = (Jr+mp*r^2)*d2t - mp*I*r*d2a + br*dt == kt*im;
% eq2 = mp*I*r*d2t -Jp*d2a + mp*g*I*alpha - bp*da == 0;
% [SS,Sbs] = odeToVectorField(eq, eq2);
% SS(4)


A = zeros(4,4); B = zeros(4,1); C = [1, 0, 0, 0; 0, 1, 0, 0]; D = 0;
A(1, 1) = 1;
A(2, 2) = 1;

A(3, 2) = (1/d)*((mp*I*r)*(mp*g*I));
A(3, 3) = -(1/d)*(Jp*(br + kt*(km/Rm)));
A(3, 4) = (1/d)*(mp*I*r*(-bp));

A(4, 2) = (1/d)*(Jr+mp*r^2)*(mp*g*I);
A(4, 3) = (1/d)*(mp*I*r*(-br-kt*(km/Rm)));
A(4, 4) = (1/d)*(Jr+mp*r^2)*(-bp);

B(3,1) = (1/d)*Jp*(kt/Rm);
B(4,1) = (1/d)*(mp*I*r*(kt*(1/Rm)));

Ad = eye(4,4)+ts*A
Bd = B*ts
