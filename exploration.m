%% Exploration LMIs

cvx_begin SDP quiet
cvx_solver SDPT3
cvx_precision high


%% Variables

dtemp=zeros(nx*nphi*nphi,nx*nphi*nphi);

variable D1((nphi)+(nx*nphi*nphi),(nphi)+(nx*nphi*nphi)) symmetric
variable D2((nphi)+(nx*nphi*nphi),(nphi)+(nx*nphi*nphi)) symmetric
variable D3((nphi)+(nx*nphi*nphi),(nphi)+(nx*nphi*nphi)) symmetric

variable tau1(1)
variable tau2(1)
variable tau3(1)

variable u(nu,L)

variable gammae(1)

Ue=[];
for i=1:L
    Ue=blkdiag(Ue,u(:,i));
end

%% Exploration energy and D>=0

Se_energy=[gammae, ones(1,L)*(Ue'); Ue*ones(L,1), gammae*eye(nu*L)];

Se_D=D1+D2+D3;

%% Exploration LMI-1: S_exploration_1
% Ue=zeros(nu*L,L);
% D1=zeros((nx*nphi)+(nx*nphi)^2);
% tau1=0;

% S_exploration_1 = S1_1 - tau1 * S1_2
% S1_1
s11_11=zeros(L*nu*nphi);
s11_12=[(1-eps)*(D_des_tilde_half')*kron(ones(1,L)*Ue',eye(nphi));zeros((nx*nphi*nphi),(L*nu*nphi))]';
s11_21=s11_12';
s11_22=-D1;

S1_1=[s11_11, s11_12; s11_21, s11_22];

% S1_2
s12_11=-eye(L*nu*nphi);
s12_12=(kron(Vx0',eye(nphi)))*(Zhat');
s12_21=s12_12';
s12_22=(Zhat)*(kron((Gamma_x-(Vx0*(Vx0'))),eye(nphi)))*(Zhat');
s12_22=(s12_22+s12_22')/2;

S1_2=[s12_11, s12_12; s12_21, s12_22];

S_exploration_1=S1_1-(tau1*S1_2);


%% Exploration LMI-2: S_exploration_2
% D2=zeros((nx*nphi)+(nx*nphi)^2);
% tau2=0;

% S_exploration_2 = S2_1 - tau1 * S2_2
% S2_1

s21_11=zeros(L*nx*nu*nphi);
s21_12=(1-eps)*kron(Ue,eye(nx*nphi))*(Zhat');
s21_21=s21_12';
s21_22= -((1-eps)/2)*(Zhat*Zhat') - ((1-eps)/eps)*barWZ - D2;

S2_1=[s21_11, s21_12; s21_21, s21_22];

% S2_2
s22_11=-eye(L*nx*nu*nphi);
s22_12=[zeros((nphi),(L*nx*nu*nphi)); kron(Vphi0,eye(nx*nphi))]';
s22_21=s22_12';
s22_22=blkdiag( zeros((nphi)) , kron( (Gamma_phi-(Vphi0*Vphi0')) , eye(nx*nphi) ) );
s22_22=(s22_22+s22_22')/2;

S2_2=[s22_11, s22_12; s22_21, s22_22];

S_exploration_2=S2_1-(tau2*S2_2);

%% Exploration LMI-3: S_exploration_3
% D311=zeros(nx*nphi);
% tau3=0;
% dtemp=zeros((nx*nphi)^2,(nx*nphi)^2);

% S_exploration_3 = S3_1 - tau3 * S3_2
% S3_1
s31_11=(1-eps)*(Ue*Utilde' + Utilde*Ue' - Utilde*Utilde');
s31_12=zeros(L*nu,nphi+(nx*nphi*nphi));
s31_21=s31_12';
s31_22_temp=-((1-eps)/eps)*barWphi - (1/T)*gamma_w*D_des_tilde;
s31_22=blkdiag(s31_22_temp,dtemp) - D3;

S3_1=[s31_11, s31_12; s31_21, s31_22];

% S3_2
s32_11= -eye(nu*L);
s32_12= [Vphi0;zeros(nx*nphi*nphi,L*nu)]';
s32_21=s32_12';
s32_22=blkdiag(Gamma_phi-(Vphi0*Vphi0'),dtemp);
s32_22=(s32_22+s32_22')/2;

S3_2=[s32_11, s32_12; s32_21, s32_22];

S_exploration_3= S3_1-(tau3*S3_2);

%%

minimize gammae

Se_energy>=0;
Se_D>=0;
% D1+D2>=0;
tau1>=0;
tau2>=0;
tau3>=0;
S_exploration_1>=0;
S_exploration_2>=0;
S_exploration_3>=0;

cvx_end

% Update for the next iteration
Utilde=Ue;
L_1=(D_des_tilde_half')*kron(ones(1,L)*(Utilde')*(Vx0'),eye(nphi));
L_2=kron((Vphi0*Utilde),eye(nx*nphi));
Zhat=[L_1;L_2];




