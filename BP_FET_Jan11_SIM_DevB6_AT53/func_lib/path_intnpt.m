function out=path_intnpt(par,Vgs,Vds,nid1,nid2)
% Vgs=Vgs-par.threshold;
%
% Author : Zhao-Yi Yan
%
% Data : 2022-07-31
%
% Copy Right Reserved
%
% Version 2
% Take-off kT division
%% derived parameters, independent with bias

%% defined parameters

acc_a=par.acc_a;
acc_r=par.acc_r;

pn=par.pn;
pp=par.pp;

eg=par.eg;
kae=par.kae;
kah=par.kah;

tempmd1=max(1e-3,abs(par.TT*nid1));
tempmd2=max(1e-3,abs(par.TT*nid2));

par.TT=tempmd1;
par1=parLic(par);

par.TT=tempmd2;
par2=parLic(par);

vt1=par1.VT;
vt2=par2.VT;

kaet=kae.*vt1;
kaht=kah.*vt2;
% s to In conversion factor
sI=par1.sI.*vt1.^-1;
tI=par2.tI.*vt2.^-1;

%% boundary conditions, denpendent with bias
%      |   x=0  |   x=1   |
xi_mat=([0+pn, Vds+pn;...          % xi_n
         0-pp, Vds-pp]-Vgs);   % xi_p
         
ff=@(x)udiag(x);
gg=@(x)vdiag(x);

i_n=sI.*integral(ff,xi_mat(1,1),xi_mat(1,2),'RelTol',acc_r,'AbsTol',acc_a);
i_p=tI.*integral(gg,xi_mat(2,1),xi_mat(2,2),'RelTol',acc_r,'AbsTol',acc_a);
out=i_n+i_p;
% xi find in serial
function u=udiag(xi_n0)
    u=xi_n0;
    for i=1:length(xi_n0)
        xi_n=xi_n0(i);
        xi=fzero(@(x)bdeq(kaet,kaht,xi_n,xi_n-eg,x),ap_sol(kae,kah,xi_n,xi_n-eg));
        u(i)=log_exp_plus((xi-xi_n)./vt1);
    end
end
function v=vdiag(xi_p0)
    v=xi_p0;
    for i=1:length(xi_p0)
        xi_p=xi_p0(i);
        xi=fzero(@(x)bdeq(kaet,kaht,xi_p+eg,xi_p,x),ap_sol(kae,kah,xi_p+eg,xi_p));
        v(i)=log_exp_plus((xi_p-xi)./vt2);
    end
end
function res=bdeq(a,b,xi_n,xi_p,xi)
   res=a.*log_exp_plus((xi-xi_n)./vt1)-b.*log_exp_plus((xi_p-xi)./vt2)+xi;
end
function xi=ap_sol(a,b,xi_n,xi_p)
    xi=((max(0,xi_p)>xi_n).*a.*xi_n+(min(0,xi_n)<xi_p).*b.*xi_p)...
        ./(1+(max(0,xi_p)>xi_n).*a+(min(0,xi_n)<xi_p).*b);
end
function y=log_exp_plus(x)
% numerical limitation is set by the exp(-inf)
bdn=-6;
y=(x<bdn).*log_exp_plus_expand(x)+(x>=bdn).*log(exp(x)+1);
y(isnan(y))=x(isnan(y));

    function y=log_exp_plus_expand(x)
        y=exp(x)-exp(2.*x)./2+exp(3.*x)./3-exp(4.*x)./4+...
            exp(5.*x)./5-exp(6.*x)./6+exp(7.*x)./7-exp(8.*x)./8;
    end
end
end