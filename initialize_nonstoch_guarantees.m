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
A=[0.49 0.49 0;0 0.49 0.49;0 0 0.49];
B=[0; 0; 0.49];

% 4x4
% A=[0.49 0.49 0 0;0 0.49 0.49 0; 0 0 0.49 0.49;0 0 0 0.49];
% B=[0; 0; 0; 0.49];

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
