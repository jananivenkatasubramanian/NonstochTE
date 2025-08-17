%% Initialization - Sequential learning and control
% clear all;
% clc;

% A, B and frequencies
% freqs=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9];
freqs=[0:floor(T/fsel)*(1/T):1-(1/T)];
% freqs=[0,0.1,0.2,0.3,0.4,0.5,0.6];
L=size(freqs,2);

% choose A, B
% 2x2
% A=[0.49 0.49; 0 0.49];
% B=[0;0.49];

% 3x3
% A=[0.49 0.49 0;0 0.49 0.49;0 0 0.49];
% B=[0; 0; 0.49];

% 4x4
% A=[0.49 0.49 0 0;0 0.49 0.49 0; 0 0 0.49 0.49;0 0 0 0.49];
% B=[0; 0; 0; 0.49];

% Chain of 2 mass-spring-dampers: 4x4
Ts=0.1; % sampling period
m1=20; % mass-1
m2=10; % mass-2
k1=3; % spring constant 1
k2=1; % spring constant 2
d1=0.01; % damping coefficient 1
d2=0.03; % damping coefficient 2
A=[1,Ts,0,0;...
  (-Ts*(k1+k2)/m2),(1-(Ts*(d1+d2)/m1)), (Ts*k2/m1), (Ts*d2/m1);...
  0,0,1,Ts;...
  (Ts*k1/m2), (Ts*d2/m2), (-Ts*k2/m2), (1-(Ts*d2/m2)) ];
B=[0;0;0;(Ts/m2)];

nx=size(A,2);
nu=size(B,2);
nphi=nx+nu;
thetatr=[A(:);B];

x0=zeros(nx,1);

% Omega_T set of possible frequencies
Omega_T=[];
for i=1:T
    Omega_T(i)=(i-1)/T;
end

eps=0.5;

