clc; clear all;

%parameters
jm=1.62e-6;
kt=0.044;
res=17;
n=6.3;
r=0.011;
mc=0.86;
jo=1.45;
m=3.7;
h1=0.025;
h2=0.075;
grav=9.81;

damping=1.6e-6;

%coefficients
a=jm*(n/r)^2+mc;
b=-(kt^2/res+damping)*(n/r)^2/a;  %xp
c=-mc*grav/a;              %theta
d=kt/res*n/r/a;            %v

e=((m*grav*h1+mc*grav*h2)+mc*h2*c)/jo;  %theta
f=-mc*grav/jo;   %x
g=(mc*h2*b)/jo;     %xp
h=(mc*h2*d)/jo;     %v

%INPUTS%
A = [0 1 0 0; 0 b c 0; 0 0 0 1; f g e 0]
B = [0; d; 0; h]
C = [0 0 1 0]
D = 0


%DESIRED ROOTS%
%POLES = -2+2.2i, -2-2.2i, -5, -5.5
desDiscrete = c2d(tf(1, [1 176 9143.14 131702 1172080]),...
    0.001, 'zoh');
Z = pole(desDiscrete)

% %State Space%
sys_Discrete = c2d(ss(A, B, C, D), 0.001, 'zoh')
K = place(sys_Discrete.A,sys_Discrete.B,Z)
N = inv(sys_Discrete.C*inv(eye(4) - (sys_Discrete.A-...
     sys_Discrete.B*K))*sys_Discrete.B)