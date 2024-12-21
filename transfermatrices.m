%% Construct \hat{V}_x and \hat{V}_\phi with prior estimates \hat{A}_0 and \hat{B}_0

Vx0=[];
Vphi0=[];

Vx_tr=[];
Vphi_tr=[];

for i=1:L
    vxi=inv(exp(1i*2*pi*freqs(i))*eye(nx)-A0)*B0;
    vphii=[vxi;eye(nu)];
    Vx0=blkdiag(Vx0,vxi);
    Vphi0=[Vphi0,vphii];
end

for i=1:L
    vxi=inv(exp(1i*2*pi*freqs(i))*eye(nx)-A)*B;
    vphii=[vxi;eye(nu)];
    Vx_tr=blkdiag(Vx_tr,vxi);
    Vphi_tr=[Vphi_tr,vphii];
end

%% Construct \hat{Y}_x and \hat{Y}_\phi with prior estimates \hat{A}_0 and \hat{B}_0

Yx0=[];
Yphi0=[];

Yx_tr=[];
Yphi_tr=[];

for i=1:L
    yxi=inv(exp(1i*2*pi*freqs(i))*eye(nx)-A0);
    yphii=[yxi;zeros(nu,nx)];
    Yx0=blkdiag(Yx0,yxi);
    Yphi0=[Yphi0,yphii];
end


for i=1:L
    yxi=inv(exp(1i*2*pi*freqs(i))*eye(nx)-A0);
    yphii=[yxi;zeros(nu,nx)];
    Yx_tr=blkdiag(Yx_tr,yxi);
    Yphi_tr=[Yphi_tr,yphii];
end
g_phi_tr=norm(Yphi_tr);
g_x_tr=norm(Yx_tr);