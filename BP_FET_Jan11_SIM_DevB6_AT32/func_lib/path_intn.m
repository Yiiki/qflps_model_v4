function out=path_intn(par,vgs0,vds0,nid)
% test parameters
Vgs=vgs0; % V
Vds=vds0; % V
acc_a=par.acc_a;
acc_r=par.acc_r;
%% defined parameters

% geometrical parameters
% material parameters
eg=par.eg; % eV
phi_n=par.pn; %eV
phi_p=par.pp; % eV

%% derived parameters, independent with bias

tempmd=max(1e-3,abs(par.TT*nid));

par.TT=tempmd;

% par.TT=par.TT*nid;
par=parLic(par);
VT=par.VT;
kapa_e=par.kae;
kapa_h=par.kah;
% s to In conversion factor
sI=par.sI;

%% boundary conditions, denpendent with bias
%      |   x=0  |   x=1   |
xi_mat=([0+phi_n, Vds+phi_n;...          % xi_n
         0-phi_p, Vds-phi_p]-Vgs)./VT;   % xi_p
         
% xi_vec=[sol_bnd(kapa_e,kapa_h,xi_mat(1,1),xi_mat(2,1)),...% x=0
%         sol_bnd(kapa_e,kapa_h,xi_mat(1,2),xi_mat(2,2))];  % x=1
% uv_mat=[log_exp_plus(xi_vec(1)-xi_mat(1,1)),log_exp_plus(xi_vec(2)-xi_mat(1,2));... % u
%         log_exp_plus(xi_mat(2,1)-xi_vec(1)),log_exp_plus(xi_mat(2,2)-xi_vec(2))];   % v
f=@(x)udiag(x);
% g=@(x)vdiag(x);

i_n=sI.*integral(f,xi_mat(1,1),xi_mat(1,2),'RelTol',acc_r,'AbsTol',acc_a);
% i_p=tI.*integral(g,xi_mat(2,1),xi_mat(2,2),'RelTol',acc_r,'AbsTol',acc_a);
% out=i_n+i_p;
out=i_n;
% xi find in serial
function u=udiag(xi_n0)
    u=xi_n0;
    for i=1:length(xi_n0)
        xi_n=xi_n0(i);
        xi=fzero(@(x)bdeq_n(kapa_e,kapa_h,xi_n,xi_n-eg/VT,x),ap_sol(kapa_e,kapa_h,xi_n,xi_n-eg/VT));
        u(i)=log_exp_plus(xi-xi_n);
    end
end
% function v=vdiag(xi_p0)
%     v=xi_p0;
%     for i=1:length(xi_p0)
%         xi_p=xi_p0(i);
%         xi=fzero(@(x)bdeq_p(kapa_e,kapa_h,xi_p+eg/VT,xi_p,x),ap_sol(kapa_e,kapa_h,xi_p+eg/VT,xi_p));
%         v(i)=log_exp_plus(xi_p-xi);
%     end
% end

end