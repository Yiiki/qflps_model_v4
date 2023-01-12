function out=path_intnp_opt_v5_tmd(pab,V1,V2,V3)
    %                                  G, D, S
    % a stand-alone function file
    % Author : Zhao-Yi Yan
    % Data : 2023-01-12
    % dev log ----------------------------------- -
    % polylog alogrithm
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
    X=3;
    kbm_mat =[...% si,ti = kb1*egVT + kb2
    %    kb1        kb2
        0.5931    2.3844    % s1
        0.1169    2.3933    % t1
        0.5893    2.9467    % s2
        0.1201    2.5495];  % t2
    % forth rank (Vgs-dependent constant) --------------------------------
%     nid = modelfun1(beta,Vgs);
%     VT=VT0.*nid; 
%     egVT=(pn+pp)./VT;
%     NcDOS=NcDOS0.*nid;
%     NvDOS=NvDOS0.*nid;
%     sI=mue.*VT.*NcDOS.*W2L.*elchg;
%     tI=muh.*VT.*NvDOS.*W2L.*elchg;
    svec = kbm_mat*[egVT;1];

    %% boundary conditions, denpendent with bias
    %      |   x=0  |   x=1   |
    xi_mat=([0+pn, Vds+pn;...          % xi_n
             0-pp, Vds-pp]-Vgs)./VT;   % xi_p
    
    sa = svec(1);ta = svec(2);sb = svec(3);tb = svec(4);
    xi0=ap_sol_sft_NIXS(X,kae,kah,xi_mat(1,1),xi_mat(2,1),sa,ta,sb,tb);
    xi1=ap_sol_sft_NIXS(X,kae,kah,xi_mat(1,2),xi_mat(2,2),sa,ta,sb,tb);
    
    x0=xi0-xi_mat(1,1);
    x1=xi1-xi_mat(1,2);
    
    y0=xi_mat(2,1)-xi0;
    y1=xi_mat(2,2)-xi1;
    
    u0=log_exp_plus(x0);
    u1=log_exp_plus(x1);
    
    v0=log_exp_plus(y0);
    v1=log_exp_plus(y1);

    eu0=exp(x0);
    eu1=exp(x1);
    
    ev0=exp(y0);
    ev1=exp(y1);
    
    Ie_squ=sI.*kae.*0.5.*(u0-u1).*(u0+u1);
    Ih_squ=tI.*kah.*0.5.*(v1-v0).*(v0+v1);
    
    Ie_exp=sI.*(polylog2_neg(-eu1)-polylog2_neg(-eu0));
    Ih_exp=tI.*(polylog2_neg(-ev0)-polylog2_neg(-ev1));
    
    Ie_cor=sI.*kah.*(corfunc(egVT,eu0)-corfunc(egVT,eu1));
    Ih_cor=tI.*kae.*(corfunc(egVT,ev1)-corfunc(egVT,ev0));

    out=Ie_squ+Ih_squ+Ie_exp+Ih_exp+Ie_cor+Ih_cor;

    dirc = 2.*(V2>V3) -1;

    out = dirc.*out.*(1+exp((phi-Vds).*phi0.^-1)).^-1;
        
    %% Section 2: parameter-scaling method with NIOS technique
    function y=corfunc(a,x)
        y=polylog2_neg(-(x+exp(-a)).*(1-exp(-a)).^-1)-polylog2_neg(-x)-log(1-exp(-a)).*log((x+exp(-a)).*(1-exp(-a)).^-1);
    end

    function y=polylog2_neg(x)
    
    % x<=0 [required]
    
    y=x;
    if ~isempty(find(x>0,1))
        error('pos x met')
    end
    xt=x(x>=-1);
    y(x>=-1)=-0.5.*log(1-xt).^2-polyLog2_func012(xt.*(xt-1).^-1);
    x2=x(x<-1);
    y(x<-1)=-pi.^2.*6.^-1-0.5.*log(1-x2).*(2.*log(-x2)-log(1-x2))+polyLog2_func012((1-x2).^-1);
    end

    function ytest=polyLog2_func012(xtest)
    P8=[
    -.424860941577076578425815513413e5
    +.183244531779306992722635513651e6
    -.327526669297550162522223586638e6
    +.313122440690220044970422956820e6
    -.172218697819825427580876931586e6
    +.544281929261565482937728842213e5
    -.927317883956725077226221332049e4
    +.723979649523300321037522923595e3
    -.168151615514046449580975569676e2
    ];
    
    Q9=[
    -.424860941577076578425815660209e5
    +.193866055318733907183294751241e6
    -.371272505998599455115517656852e6
    +.387055275261533981150821382566e6
    -.238147200659810439902217446162e6
    +.875889100897049125433126845284e5
    -.185677658098068966095883646794e5
    +.205632929098181836999766662708e4
    -.954141875991525514141234686148e2
    +1
    ];
    
    P8=flipud(P8);
    Q9=flipud(Q9);
    
    % xtest=linspace(-1,0,1e3);
    % tic
    ytest=xtest.*polyval(P8,xtest).*polyval(Q9,xtest).^-1;
    % toc
    % tic
    % ytest0=polylog(2,xtest);
    % toc
    % figure
    % plot(xtest,ytest0,'o',xtest,ytest,'-')
    
    
    end


    function xi=ap_sol_sft_NIXS(X,a,b,xi_n,xi_p,sa,ta,sb,tb)
    %     xi=((max(0,xi_p)>xi_n).*a.*xi_n+(min(0,xi_n)<xi_p).*b.*xi_p)...
    %         ./(1+(max(0,xi_p)>xi_n).*a+(min(0,xi_n)<xi_p).*b);
        xi = (fq(xi_n,xi_p,sa,ta,b).*a.*xi_n+fq(-xi_p,-xi_n,sb,tb,a).*b.*xi_p)...
            ./(1+fq(xi_n,xi_p,sa,ta,b).*a+fq(-xi_p,-xi_n,sb,tb,a).*b);
    
        for i=1:X
        % X-Shot Newton-Iteration
        xi = xi - F(a,b,xi_n,xi_p,xi).*dFdx(a,b,xi_n,xi_p,xi).^-1;
        end
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