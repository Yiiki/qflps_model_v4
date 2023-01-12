function out=path_intnp_opt_v3_tmd(pab,V1,V2,V3)
    %                                  G, D, S
    % a stand-alone function file
    % Author : Zhao-Yi Yan
    % Data : 2022-12-19
    % dev log ----------------------------------- -
    % pab-acceptable two-methods (svec+cardinal) enable ver.
    %_____________________________________________________________________

    Vs=min(V2,V3);
    Vd=max(V2,V3);
    Vgs=V1-Vs;
    Vds=Vd-Vs;
    
    %% defined parameters
%     acc_a=pab(1);
%     acc_r=pab(2);
    pn=pab(3);
    pp=pab(4);
    VT=pab(5);
    kae=pab(6);
    kah=pab(7);
    egVT=pab(8);
    sI=pab(9);
    tI=pab(10);
    phi=pab(11);
    phi0=pab(12);


    %% device specified parameters
    % first rank ---------------------------------------------------------
%     mue=3.6366e-04;% mu_e [m^2 V^-1 s^-1]
%     muh=0.0016;% mu_h [m^2 V^-1 s^-1]
%     pn=2.4155;  % Vth,n [V]
%     pp=-0.2696; % Vth,p [V]
%     kae=24.9195;% Cq,e/Cox
%     kah=23.2582;% Cq,h/Cox
    % second rank --------------------------------------------------------
%     VT0 = 0.0259;% kT/q [V] at 300 K
%     elchg=1.6022e-19;% elementary charge [C]
%     W2L=1.0587;% W/L
%     NcDOS0=1.6199e+16;% [m^-2] NcDOS at 300 K
%     NvDOS0=1.5119e+16;% [m^-2] NcDOS at 300 K
    % third rank ---------------------------------------------------------
%     beta = [...
%         7.3534
%         5.1961
%         1.2351
%         2.5424];
    kbm_mat =[...% si,ti = kb1*egVT + kb2
    %    kb1        kb2
        0.4659    5.6248    % s1
        0.1477    1.6150    % t1
        0.4541    6.3953    % s2
        0.1530    1.7150];  % t2
    % forth rank (Vgs-dependent constant) --------------------------------
%     nid = modelfun1(beta,Vgs);
%     VT=VT0.*nid; 
%     egVT=(pn+pp)./VT;
%     NcDOS=NcDOS0.*nid;
%     NvDOS=NvDOS0.*nid;
%     sI=mue.*VT.*NcDOS.*W2L.*elchg;
%     tI=muh.*VT.*NvDOS.*W2L.*elchg;
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
    
%     function y=modelfun1(b,x)
%     y = b(4)+ b(1).*exp(-((x-b(2))./b(3)).^2);
%     end
    
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
    
end