function out=path_intnp_opt_v3_tnd_pub(lab,optvec,args_pts,sigval,V1,V2,V3)
%                                                                  G, D, S
% a stand-alone function file
% Author : Zhao-Yi Yan
% 
% dev log ----------------------------------- Data : 2022-12-19
% pab-acceptable two-methods (svec+cardinal) enable ver.
%
% dev log ----------------------------------- Data : 2022-12-28
% optvec enable ver.
%_____________________________________________________________________
vth2Ntrp=2.5145e+12;% V --> cm^-2
nid2phit=0.026;% 1--> V
Vpn=0.2-2;
Vpp=0.19+2;

Vs=min(V2,V3);
Vd=max(V2,V3);
Vgs=V1-Vs;
Vds=Vd-Vs;

% check

if length(args_pts)~=7
    error('invalid args_pts length')
end

if sum(args_pts)~=length(optvec)
    error('invalid optvec met')
end

% parm gen
parm=parmf(optvec,Vgs,args_pts,sigval);
ue=parm(1).*1e-4;
pn=parm(2).*vth2Ntrp.^-1+Vpn;
uh=parm(3).*1e-4;
pp=parm(4).*vth2Ntrp.^-1+Vpp;
nid=parm(5).*nid2phit.^-1;
phi=parm(6);
phi0=parm(7);

%% device specified parameters
% first rank ---------------------------------------------------------
%     mue=3.6366e-04;% mu_e [m^2 V^-1 s^-1]
%     muh=0.0016;% mu_h [m^2 V^-1 s^-1]
%     pn=2.4155;  % Vth,n [V]
%     pp=-0.2696; % Vth,p [V]
kae=24.9195;% Cq,e/Cox
kah=23.2582;% Cq,h/Cox
% second rank --------------------------------------------------------
VT0 = 0.0259;% kT/q [V] at 300 K
elchg=1.6022e-19;% elementary charge [C]
lwv=[2.8153464./1.6405017,2.8708355./1.5427191,2.9038664./1.4809202,...
    3.3380508./1.144557,3.22944./1.5943659].^-1;
W2L=sum(lwv(lab(1):lab(2)-1)).^-1;% W/L
NcDOS0=1.6199e+16;% [m^-2] NcDOS at 300 K
NvDOS0=1.5119e+16;% [m^-2] NcDOS at 300 K
% third rank ---------------------------------------------------------
kbm_mat =[...% si,ti = kb1*egVT + kb2
%    kb1        kb2
0.4659    5.6248    % s1
0.1477    1.6150    % t1
0.4541    6.3953    % s2
0.1530    1.7150];  % t2
% forth rank (Vgs-dependent constant) --------------------------------
VT=VT0.*nid; 
egVT=(pn+pp)./VT;
NcDOS=NcDOS0.*nid;
NvDOS=NvDOS0.*nid;
sI=ue.*VT.*NcDOS.*W2L.*elchg;
tI=uh.*VT.*NvDOS.*W2L.*elchg;
svec = kbm_mat*[egVT;1];
% fifth rank ---------------------------------------------------------
N=10;h=1;

%% boundary conditions, denpendent with bias
%      |   x=0  |   x=1   |
xi_mat=([0+pn, Vds+pn;...          % xi_n
0-pp, Vds-pp]-Vgs)./VT;   % xi_p

sa = svec(1);ta = svec(2);sb = svec(3);tb = svec(4);
xi0=ap_sol_sft_NIOS(kae,kah,xi_mat(1,1),xi_mat(2,1),sa,ta,sb,tb);
xi1=ap_sol_sft_NIOS(kae,kah,xi_mat(1,2),xi_mat(2,2),sa,ta,sb,tb);

x0=xi0-xi_mat(1,1);
x1=xi1-xi_mat(1,2);

y0=xi_mat(2,1)-xi0;
y1=xi_mat(2,2)-xi1;

u0=log_exp_plus(x0);
u1=log_exp_plus(x1);

v0=log_exp_plus(y0);
v1=log_exp_plus(y1);

Ie_squ=sI.*kae.*0.5.*(u0-u1).*(u0+u1);
Ih_squ=tI.*kah.*0.5.*(v1-v0).*(v0+v1);

w=(-N:N).*h;
Ie1 = sI.*h.*sum(fseh(w,x0,x1,kae,egVT));
Ih1 = tI.*h.*sum(fseh(w,y1,y0,kah,egVT));

out=Ie_squ+Ih_squ+Ie1+Ih1;

dirc = 2.*(V2>V3) -1;

out = dirc.*out.*(1+exp((phi-Vds).*phi0.^-1)).^-1;

function y=fseh(w,x0,x1,ka,egVT)
y = 0.5.*(x0-x1).*zfeh(xf(x1,x0,zf(w)),ka,egVT).*dzf(w);
end

function y=zfeh(x,ka,egVT)
y = log_exp_plus(x).*(1 + ka.*(1+exp(egVT+x)).^-1);
end

function x=xf(a,b,z)
x=0.5.*(b-a).*(z+1) + a;
end

function z=zf(w)
z=(exp(w)-1).*(exp(w)+1).^-1;
end

function dzdw=dzf(w)
dzdw = 2.*exp(w).*(exp(w)+1).^-2;
end

%% Section 2: parameter-scaling method with NIOS technique
function xi=ap_sol_sft_NIOS(a,b,xi_n,xi_p,sa,ta,sb,tb)
xi=(fq(xi_n,xi_p,sa,ta,b).*a.*xi_n+fq(-xi_p,-xi_n,sb,tb,a).*b.*xi_p)...
./(1+fq(xi_n,xi_p,sa,ta,b).*a+fq(-xi_p,-xi_n,sb,tb,a).*b);
% One-Shot Newton-Iteration
xi = xi - F(a,b,xi_n,xi_p,xi).*dFdx(a,b,xi_n,xi_p,xi).^-1;
end

function z = F(a,b,xi_n,xi_p,xi)
z = a.*log_exp_plus(xi-xi_n) - b.*log_exp_plus(xi_p-xi) + xi;
end

function z = dFdx(a,b,xi_n,xi_p,xi)
z = a.*(1+exp(-xi+xi_n)).^-1 + b.*(1+exp(-xi_p+xi)).^-1 + 1;
end

function  z = fq(x,y,s,t,a)
z = (exp((x-s.*log_exp_plus(a.*(a+1).^-1.*y.*s.^-1)).*t.^-1) + 1 ).^-1;
end

%% y = ln(1+exp(x))
function y=log_exp_plus(x)

% 'bdn' is set according to exp(x) == inf
bdn=-6;
y(x<bdn)=log_exp_plus_expand(x(x<bdn));
y(x>=bdn)=log(exp(x(x>=bdn))+1);

function y=log_exp_plus_expand(x)
y=exp(x)-exp(2.*x)./2+exp(3.*x)./3-exp(4.*x)./4+...
exp(5.*x)./5-exp(6.*x)./6+exp(7.*x)./7-exp(8.*x)./8;
end
end

%% version 3
function parm=parmf(optvec,vgsbch,args_pts,sigval)

% sigval=0.1;

shape_tag=0;

if size(vgsbch,1)>1% row vgsbch preferred
    vgsbch=vgsbch';
    shape_tag=1;
end


for i=1:7
    uevec_parm=optvec(1:args_pts(1));
    vthnvec_parm=optvec(args_pts(1)+1:sum(args_pts(1:2)));
    uhvec_parm=optvec(sum(args_pts(1:2))+1:sum(args_pts(1:3)));
    vthpvec_parm=optvec(sum(args_pts(1:3))+1:sum(args_pts(1:4)));
    nidvec_parm=optvec(sum(args_pts(1:4))+1:sum(args_pts(1:5)));
    phivec_parm=optvec(sum(args_pts(1:5))+1:sum(args_pts(1:6)));
    phi0vec_parm=optvec(sum(args_pts(1:6))+1:sum(args_pts(1:7)));
end

vgsvec_parm_ue=linspace(-3,3,args_pts(1));
vgsvec_parm_vthn=linspace(-3,3,args_pts(2));
vgsvec_parm_uh=linspace(-3,3,args_pts(3));
vgsvec_parm_vthp=linspace(-3,3,args_pts(4));
vgsvec_parm_nid=linspace(-3,3,args_pts(5));
vgsvec_parm_phi=linspace(-3,3,args_pts(6));
vgsvec_parm_phi0=linspace(-3,3,args_pts(7));


uevec=wtbnd_gaus(vgsvec_parm_ue,uevec_parm,vgsbch,sigval);
vthnvec=wtbnd_gaus(vgsvec_parm_vthn,vthnvec_parm,vgsbch,sigval);
uhvec=wtbnd_gaus(vgsvec_parm_uh,uhvec_parm,vgsbch,sigval);
vthpvec=wtbnd_gaus(vgsvec_parm_vthp,vthpvec_parm,vgsbch,sigval);
nidvec=wtbnd_gaus(vgsvec_parm_nid,nidvec_parm,vgsbch,sigval);
phivec=wtbnd_gaus(vgsvec_parm_phi,phivec_parm,vgsbch,sigval);
phi0vec=wtbnd_gaus(vgsvec_parm_phi0,phi0vec_parm,vgsbch,sigval);


parm=[
    uevec
    vthnvec
    uhvec
    vthpvec
    nidvec
    phivec
    phi0vec
];

if shape_tag==1% recover shape
    parm=parm';
end

end

function out=wtbnd_gaus(faket,pvec,tinq,sig)

if size(pvec,1)==1% unified as collumn vector
pvec=pvec';
end

if size(faket,1)==1% unified as collumn vector
faket=faket';
end

pvec=inflect_extrapolate_func(pvec);
faket=inflect_extrapolate_func(faket);

pvect=pvec-pvec(1);

out=tinq;

for i=1:length(tinq)
    out(i)=inptdat(sig,faket,pvect,tinq(i));
end


out=out+pvec(1);


end

function yinv=inflect_extrapolate_func(y)

yinv=[2.*y(1)-flip(y(2:end));y;2.*y(end)-flip(y(1:end-1))];

end

function y=inptdat(sig,faket,pvec,tinq)
% y=integral(@(t) gaus(sig,tinq-t).*interp1(faket,pvec,t),min(faket),max(faket));
uvec=(faket-tinq).*(sqrt(2).*sig).^-1;
amu=diff(pvec).*diff(faket).^-1;
bmu=pvec(1:end-1)-amu.*uvec(1:end-1).*(sqrt(2).*sig);
erfvec=erf111s(uvec);
exp2vc=exp(-uvec.^2);
y=0.5.*bmu'*diff(erfvec)-sig.*(2.*pi).^-0.5.*amu'*diff(exp2vc);

% rectfiying
y=sign(faket(end)-faket(1)).*y;
end

function y=erf111s(x)
y=x;

xa=0.46875;
xb=4.0;

xan=-xa;
xbn=-xb;

for i=1:length(x)
    xi=x(i);
    if abs(xi)<= xa
        yi=AR11(xi);
    elseif xi>xa && xi<xb
        yi=1-BR11(xi);
    elseif xi>=xb
        yi=1-CR11(xi);
    elseif xi<xan && xi>xbn
        yi=BR11(-xi)-1;
    else
        yi=CR11(-xi)-1;
    end
    y(i)=yi;
end

function y=AR11(x)
% p 
p0=3.6767877e00;
p1=-9.7970465e-02;
% q
q0=3.2584593e00;
q1=1.0e00;
% numerator
nu=p0.*x.^0+p1.*x.^2;
% denominator
de=q0.*x.^0+q1.*x.^2;
% y
y=x.*nu.*de.^-1;
end

function y=BR11(x)
% p 
p0=7.3033e-01;
p1=-2.3877e-02;
% q
q0=6.6211e-01;
q1=1.0e00;
% numerator
nu=p0.*x.^0+p1.*x.^1;
% denominator
de=q0.*x.^0+q1.*x.^1;
% y
y=exp(-x.^2).*nu.*de.^-1;
end

function y=CR11(x)
% p 
p0=-1.24368544e-01;
p1=-9.68210364e-02;
% q
q0=4.40917061e-01;
q1=1.0e00;
% numerator
nu=p0.*x.^-0+p1.*x.^-2;
% denominator
de=q0.*x.^-0+q1.*x.^-2;
% y
y=exp(-x.^2).*x.^-1.*(pi.^-0.5+x.^-2.*nu.*de.^-1);
end

end
% function y=AR44(x)
% p0=3.209377589138469472562;
% p1=3.209377589138469472562;



end