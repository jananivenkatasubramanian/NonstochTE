%% Transfer functions Y_x and Y_\phi bounds:
% Determine bound: ||Y_{x,tr}|| <= \gamma_x
% Determine bound: ||Y_{\phi,tr}|| <= \gamma_\phi

Yx_s=[];
Yphi_s=[];
Yx_norm=[];
Yphi_norm=[];
for i=1:nx:nx*Ns
    pa=params(i:i+nx-1,1:nx);
    pb=params(i:i+nx-1,nx+1:nphi);
    for j=1:L
        ytemp=inv(exp(1i*2*pi*freqs(j))*eye(nx)-pa);
        Yx_s=blkdiag(Yx_s,ytemp);
        Yphi_s=[Yphi_s,[ytemp;zeros(nu,nx)]];
    end
    Yx_norm=[Yx_norm,norm(Yx_s)];
    Yphi_norm=[Yphi_norm,norm(Yphi_s)];
    Yx_s=[];
    Yphi_s=[];
end

gamma_x=max(Yx_norm);
gamma_phi=max(Yphi_norm);

%% Transfer functions Y_\phi bounds: Robust control bound
% % \tilde{Y}_\phi \tilde{Y}_\phi^H <= Gyphi

% Yphi_s=[];
% cvx_begin SDP quiet
% cvx_precision high
% 
% variable Gyphitilde(nphi,nphi) semidefinite
% 
% minimize trace(Gyphitilde)
% subject to
% for i=1:nx:nx*Ns
%     pa=params(i:i+nx-1,1:nx);
%     pb=params(i:i+nx-1,nx+1:nphi);
%     for j=1:L
%         ytemp=inv(exp(1i*2*pi*freqs(j))*eye(nx)-pa);
%         Yphi_s=[Yphi_s,[ytemp;zeros(nu,nx)]];
%     end
%     Gyphitilde-((Yphi_s-Yphi0)*(Yphi_s-Yphi0)')>=0;
%     Yphi_s=[];
% end
% 
% cvx_end

%% Transfer functions Y_\x bounds: Robust control bound

% Yx_s=[];
% cvx_begin SDP quiet
% cvx_precision high
% 
% variable Gyx(nx*L,nx*L) semidefinite
% 
% minimize trace(Gyx)
% subject to
% for i=1:nx:nx*Ns
%     pa=params(i:i+nx-1,1:nx);
%     pb=params(i:i+nx-1,nx+1:nphi);
%     for j=1:L
%         ytemp=inv(exp(1i*2*pi*freqs(j))*eye(nx)-pa);
%         Yx_s=blkdiag(Yx_s,ytemp);
%     end
%     Gyx-((Yx_s)*(Yx_s)')>=0;
%     Yx_s=[];
% end
% 
% cvx_end

%% Transfer functions Y_\phi bounds: Robust control bound
% Y_\phi Y_\phi^H <= Gyphi

Yphi_s=[];
cvx_begin SDP quiet
cvx_precision high

variable Gyphi(nphi,nphi) semidefinite

minimize trace(Gyphi)
subject to
for i=1:nx:nx*Ns
    pa=params(i:i+nx-1,1:nx);
    pb=params(i:i+nx-1,nx+1:nphi);
    for j=1:L
        ytemp=inv(exp(1i*2*pi*freqs(j))*eye(nx)-pa);
        Yphi_s=[Yphi_s,[ytemp;zeros(nu,nx)]];
    end
    Gyphi-((Yphi_s)*(Yphi_s)')>=0;
    Yphi_s=[];
end

cvx_end

%% Transfer functions V_x bounds: Robust control bound
% \tilde{V}_x \tilde{V}_x^H <= \Gamma_x  

Vx_s=[];
cvx_begin SDP quiet
% cvx_solver mosek
cvx_precision high

variable Gamma_x((nx*L),(nx*L)) semidefinite

minimize trace(Gamma_x)
subject to
for i=1:nx:nx*Ns
    pa=params(i:i+nx-1,1:nx);
    pb=params(i:i+nx-1,nx+1:nphi);
    for j=1:L
        vtemp=inv(exp(1i*2*pi*freqs(j))*eye(nx)-pa)*pb;
        Vx_s=blkdiag(Vx_s,vtemp);
    end
    Gamma_x-((Vx_s-Vx0)*(Vx_s-Vx0)')>=0;
    Vx_s=[];
end

cvx_end

%% Transfer functions V_\phi bounds: Robust control bound
% \tilde{V}_\phi \tilde{V}_\phi^H <= \Gamma_\phi 

Vphi_s=[];
cvx_begin SDP quiet
% cvx_solver mosek
cvx_precision high

variable Gamma_phi((nx+nu),(nx+nu)) semidefinite

minimize trace(Gamma_phi)
subject to
for i=1:nx:nx*Ns
    pa=params(i:i+nx-1,1:nx);
    pb=params(i:i+nx-1,nx+1:nphi);
    for j=1:L
        vtemp=inv(exp(1i*2*pi*freqs(j))*eye(nx)-pa)*pb;
        Vphi_s=[Vphi_s,[vtemp;eye(nu)]];
    end
    Gamma_phi-((Vphi_s-Vphi0)*(Vphi_s-Vphi0)')>=0;
    Vphi_s=[];
end
cvx_end








