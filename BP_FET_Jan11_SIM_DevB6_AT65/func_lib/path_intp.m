function out=path_intp(par,Vgs,Vds,nid)
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

tempmd=max(1e-3,abs(par.TT*nid));

par.TT=tempmd;

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
tI=par.tI;

%% boundary conditions, denpendent with bias
%      |   x=0  |   x=1   |
xi_mat=([0+pn, Vds+pn;...          % xi_n
         0-pp, Vds-pp]-Vgs)./VT;   % xi_p
         

gg=@(x)vdiag(x);
% gg =@(x) arrayfun(@vdiag,gpuArray(x));
% gg =@(x) arrayfun(@vdiag,x);

i_p=tI.*integral(gg,xi_mat(2,1),xi_mat(2,2),'RelTol',acc_r,'AbsTol',acc_a);

out=i_p;

% xi find in serial

function v=vdiag(xi_p0)
    v=xi_p0;
    for i=1:length(xi_p0)
        xi_p=xi_p0(i);
        xi=fzero(@(x)bdeq(kae,kah,xi_p+egVT,xi_p,x),ap_sol(0,kah,1e4,xi_p));
        v(i)=log_exp_plus(xi_p-xi);
    end
end

% function v=gg(x)
% y=gpuArray(x);
% v=arrayfun(vdiag,y);
% v=gather(v);
% end
% 
% function v=vdiag(xi_p)
% 
%         xi=fzero(@(x)bdeq(kae,kah,xi_p+egVT,xi_p,x),ap_sol(0,kah,1e4,xi_p));
%         v=log_exp_plus(xi_p-xi);
% 
% end

function res=bdeq(~,b,~,xi_p,xi)
   res=-b.*log_exp_plus(xi_p-xi)+xi;
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