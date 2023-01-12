function pab=parLiz_FET(ue,vthn,uh,vthp,nid,phi,phi0)
%% overwrite
pa0=parLib;

pa0.ue=ue;
pa0.uh=uh;

pa0.pn=vthn;
pa0.pp=vthp;

temp0=pa0.TT;
tempmd=max(1e-3,abs(temp0*nid));
pa0.TT = tempmd;

pa1=parLic(pa0);
VT=pa1.VT;

pab=[
    pa1.acc_a;                      % 01
    pa1.acc_r;                      % 02
    pa1.pn;                         % 03
    pa1.pp;                         % 04
    VT;                             % 05
    pa1.kae;                        % 06
    pa1.kah;                        % 07
    pa1.eg/VT;                      % 08
    pa1.sI;                         % 09
    pa1.tI;                         % 10
    phi;                            % 11
    phi0                            % 12
    ];

end