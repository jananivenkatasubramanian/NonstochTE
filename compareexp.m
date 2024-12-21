% Simulate with targeted inputs

U_ns=zeros(L,T);
tset=1:T;
for i=1:L
    cosi=cos(2*pi*freqs(i)*tset);
    U_ns(i,:)=Ue(i,i)*cosi;
end
U_ns=sum(U_ns);


% Compute excitation
X_ns=zeros(nx,T+1);
Xns=[];
Phi_ns=zeros(nphi,T);

x0=zeros(nx,1);
X_ns(:,1)=x0;
for m=1:T
    Phi_ns(:,m)=[X_ns(:,m);U_ns(:,m)];
%     w_ns=(sqrt(gamma_w/T)*[zeros(nx-1,1);-cos(X_ns(1,m))]);
    w_ns=(sqrt(gamma_w/T)*[-cos(X_ns(1,m));zeros(nx-1,1)]);
%     w_ns=sqrt(gamma_w/(T*L))*(-cos(X_ns(:,m)));
%     w_ns=[sqrt(gamma_w)/T;zeros(nx-1,1)];
%     w_ns=w(m);
    xns=A*X_ns(:,m)+B*U_ns(:,m)+w_ns;
    X_ns(:,m+1)=xns;
    Xns=[Xns;xns];
end

D_ns=Phi_ns*Phi_ns';

P_ns=kron(inv(D_ns),eye(nx));

theta_ns=P_ns*kron(Phi_ns,eye(nx))*Xns;

G_ns=gamma_w+(theta_ns'*inv(P_ns)*theta_ns)-(Xns'*Xns);

Dbar_ns=G_ns*P_ns;

normGP=norm(Dbar_ns);

nDdesinv=norm(inv(D_des));

err_ns=sqrt((thetatr-theta_ns)'*(thetatr-theta_ns));

werr_ns=(thetatr-theta_ns)'*(D_des)*(thetatr-theta_ns);

