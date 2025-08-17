 %% Main --- Non-stoch targeted exploration
clear all;
clc;
simpars=[];
%%
% Exploration time
T=100;
fsel=20;

run initialize_nonstoch_guarantees.m

% Set initial bound on estimate D0 - with structure so that  bounds on \hat{A}_0 
% and \hat{B}_0 can be computed 
D0tilde=1e6*eye(nphi);
D0=kron(D0tilde,eye(nx));
D0inv=inv(D0);

% initial estimate \hat{A}_0 and \hat{B}_0
theta0=thetatr + 4e-4;
theta0_check=(theta0-thetatr)'*D0*(theta0-thetatr);
if theta0_check>1
    disp("Initial uncertainty bound not satisfied!");
    disp(theta0_check);
else
    disp(theta0_check);
end

p0=reshape(theta0,[nx,nphi]);
A0=p0(:,1:nx);
B0=p0(:,nx+1);

% Generate \theta samples to derive scenario bounds
thetas=[];
params=[];
n=0;
Ns=20000; % set number of samples based on required confidence
cbar=chi2inv(0.9,(nx*(nx+nu))); %to scale D0 - to get samples from normrnd
while n<Ns
    t1=(mvnrnd(theta0,D0inv/cbar))';
    t2=(theta0-t1)'*D0*(theta0-t1);
    if t2<=1
        thetas=[thetas,t1];
        paramt=reshape(t1,[nx,nphi]);
        params=[params;paramt];
        n=n+1;
    end
end

% Construct transfer matrices
run transfermatrices.m

% Scenario bounds gamma_x, gamma_phi, Gamma_x, Gamma_phi, l1, l2
run thetasamples.m

disp(eig(Gamma_phi-Vphi0*Vphi0'));


%% -----------------------------------------------------------------------
% exploration time
% mul=1;
T=100;
gw=[1e-4,1e-3,1e-2,1e-1,1e0,1e1]; %size 6

% Noise energy bound
gamma_w=gw(6);

%Set initial U_e and \tilde{L} and D_{des}

% Desired bound D_des
D_des_tilde=1e4*eye(nphi);
D_des=kron(D_des_tilde,eye(nx));
D_des_half=chol(D_des);
D_des_tilde_half=chol(D_des_tilde);

% Bounds on noise spectral terms \bar{W}_\phi and \bar{W}_Z

barWphi=(gamma_w/T)*Gyphi;

l2_1=norm(D_des)*(gamma_x^2)*(gamma_w)/T;
l2_2=(gamma_w/T)*(gamma_phi^2);
l2=l2_1+l2_2;
barWZ=l2*eye((nphi)+(nx*nphi*nphi));

% -----------------------------------------------------------------------
% Set initial candidate \tilde{U}
Utilde=[];
ut=1.5e2*ones(nu,L);
for i=1:L
    Utilde=blkdiag(Utilde,ut(:,i));
end


Z_1=(D_des_tilde_half')*kron(ones(1,L)*(Utilde')*(Vx0'),eye(nphi));
Z_2=kron((Vphi0*Utilde),eye(nx*nphi));
Zhat=[Z_1;Z_2];
Lpinv=pinv(Zhat);

% Non-stochastic exploration: run iterations to reduce suboptimality
% from convex-relaxation
gtemp=Inf;
for i=1:20
    run exploration.m
    if(abs((gtemp-gammae)/gtemp))<1e-2
        break;
    end
    gtemp=gammae;
end

% Compute simulation parameters after exploration

run compareexp.m

%%
% simpars=[]; % set to interested values
% temp=[D0(1,1),D_des_tilde(1,1),gamma_w,T,err_ns,werr_ns,G_ns,normGP,nDdesinv];
% simpars=[simpars;temp];



