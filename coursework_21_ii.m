syms theta(t) alpha(t) u t Rm kt km mr r Jr br mp Lp I Jp bp g
dt = diff(theta);
d2t = diff(theta, 2);
da = diff(alpha);
d2a = diff(alpha,2);

ts = 0.002;

% Rm = 7.5;
% kt = 0.042;
% km = 0.042;
% 
% mr = 0.095;
% r = 0.085;
% Jr = 2.41e-4;
% br = 1e-3;
% 
% mp = 0.024;
% Lp = 0.129;
% I = 0.0645;
% Jp = 1.33e-4;
% bp = 5e-5;
% 
% g = 9.81;


d = Rm*(Jp*(Jr+mp*r^2)-(mp*I*r)^2);

eq = (Jr+mp*r^2)*d2t - mp*I*r*d2a + br*dt == kt*im;
eq2 = mp*I*r*d2t -Jp*d2a + mp*g*I*alpha - bp*da == 0;
eq3 = u - Rm*im - km*dt == 0;
[SS,Sbs] = odeToVectorField(eq, eq2,eq3)

%%

A = zeros(4,4); B = zeros(4,1); C = [1, 0, 0, 0; 0, 1, 0, 0]; D = 0;
A(1, 1) = 1;
A(2, 2) = 1;

A(3, 2) = (1/d)*((mp*I*r)*(mp*g*I));
A(3, 3) = -(1/d)*(Jp*(br + kt*(km/Rm)));
A(3, 4) = (1/d)*(Jp*(kt*1/Rm));

A(4, 2) = (1/d)*(Jr+mp*r^2)*(mp*g*I);
A(4, 3) = (1/d)*(mp*I*r*(-br-kt*(km/Rm)));
A(4, 4) = (1/d)*(Jr+mp*r^2)*(-bp);

B(3,1) = (1/d)*Jp*(kt/Rm);
B(4,1) = (1/d)*(mp*I*r*(kt*(1/Rm)));

Ad = eye(4,4)+A;
Bd = B*ts;