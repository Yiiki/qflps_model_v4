function out=path_intnp_bpt_v2_fsd(par,a,rsd,phi,phi0,Vgs,Vds,nid)

% Vgs=Vgs-par.threshold;
%
% Author : Zhao-Yi Yan
%
% Data : 2021-07-22
%
% Internal Learning Only, Copy Right Reserved
%
% test parameters
% Vgs=vgs0; % V
% Vds=vds0; % V
%% derived parameters, independent with bias

tempmd=max(1e-3,abs(par.TT*nid));

par.TT=tempmd;

% par.TT=par.TT*nid;
par=parLic(par);

%% defined parameters
acc_a=par.acc_a;
acc_r=par.acc_r;
pn=par.pn;
pp=par.pp;
VT=par.VT;
kae=par.kae;
kah=par.kah;
egVT=par.eg/VT;
% s to In conversion factor
sI=par.sI;
tI=par.tI;
An=par.An;
%% boundary conditions, denpendent with bias
%      |   x=0  |   x=1   |
xi_mat=([0+pn, Vds+pn;...          % xi_n
         0-pp, Vds-pp]-Vgs)./VT;   % xi_p

xi0=fzero(@(x)bdeq(kae,kah,xi_mat(1,1),xi_mat(2,1),x),ap_sol(kae,kah,xi_mat(1,1),xi_mat(2,1)));
xi1=fzero(@(x)bdeq(kae,kah,xi_mat(1,2),xi_mat(2,2),x),ap_sol(kae,kah,xi_mat(1,2),xi_mat(2,2)));

% xi0=1;
% xi1=1;

x0=xi0-xi_mat(1,1);
x1=xi1-xi_mat(1,2);

y0=xi_mat(2,1)-xi0;
y1=xi_mat(2,2)-xi1;

u0=log_exp_plus(x0);
u1=log_exp_plus(x1);

v0=log_exp_plus(y0);
v1=log_exp_plus(y1);

Ie_squ=sI.*kae.*0.5.*(u0.^2-u1.^2);
Ih_squ=tI.*kah.*0.5.*(v1.^2-v0.^2);

gg=@(x) log_exp_plus(x);

ff=@(x) log_exp_plus(x).*(1+exp(egVT+x)).^-1;
% 
% Ie_exp=sI.*integral(gg,x1,x0,'RelTol',acc_r,'AbsTol',acc_a);
% Ih_exp=tI.*integral(gg,y0,y1,'RelTol',acc_r,'AbsTol',acc_a);


Ie_cor=sI.*integral(@(x) ff(x).*kah + gg(x),x1,x0,'RelTol',acc_r,'AbsTol',acc_a);
Ih_cor=tI.*integral(@(x) ff(x).*kae + gg(x),y0,y1,'RelTol',acc_r,'AbsTol',acc_a);

out0=Ie_squ+Ih_squ+Ie_cor+Ih_cor;
% out0=(1 + rsd.*An.*turonf(a,Vgs-pn)).^-1.*(Ie_squ+Ie_cor)+Ih_squ+Ih_cor;

out = out0.*(1+exp((phi-0).*phi0.^-1)).*(1+exp((phi-Vds).*phi0.^-1)).^-1;

% xi find in serial

function y=turonf(a,x)
    y = a.*log_exp_plus(x.*a.^-1);
end

function res=bdeq(a,b,xi_n,xi_p,xi)
   res=a.*log_exp_plus(xi-xi_n)-b.*log_exp_plus(xi_p-xi)+xi;
end
function xi=ap_sol(a,b,xi_n,xi_p)
    xi=((max(0,xi_p)>xi_n).*a.*xi_n+(min(0,xi_n)<xi_p).*b.*xi_p)...
        ./(1+(max(0,xi_p)>xi_n).*a+(min(0,xi_n)<xi_p).*b);
end
function y=log_exp_plus(x)
% numerical limitation is set by the exp(-inf)
bdn=-6;
y(x<bdn)=log_exp_plus_expand(x(x<bdn));
y(x>=bdn)=log(exp(x(x>=bdn))+1);

    function y=log_exp_plus_expand(x)
        y=exp(x)-exp(2.*x)./2+exp(3.*x)./3-exp(4.*x)./4+...
            exp(5.*x)./5-exp(6.*x)./6+exp(7.*x)./7-exp(8.*x)./8;
    end
end

% function y=dilog_exp_plus(x)
% 
% % numerical limitation is set by the exp(-inf)
% bdn=-6;
% y(x<bdn)=dilog_exp_plus_expand(x(x<bdn));
% y(x>=bdn)=-dilog(1+exp(x(x>=bdn)));
% % y(isnan(y))=x(isnan(y));
% 
%     function y=dilog_exp_plus_expand(x)
%         y=exp(x)-exp(2.*x)./4+exp(3.*x)./9-exp(4.*x)./16+exp(5.*x)./25-exp(6.*x)./36+exp(7.*x)./49-exp(8.*x)./64;
%     end
% end

end