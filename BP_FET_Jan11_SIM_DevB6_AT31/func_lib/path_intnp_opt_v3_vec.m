function idsvec=path_intnp_opt_v3_vec(optparm,Vgs,vdsmin,vdsmax,vdspts,chnl)

ue=optparm(1);
vthn=optparm(2);
uh=optparm(3);
vthp=optparm(4);
nid=optparm(5);
phi=optparm(6);
phi0=optparm(7);

pab=parLiz_FET(ue,vthn,uh,vthp,nid,phi,phi0);

vdslis=linspace(vdsmin,vdsmax,vdspts)';

idsvec=vdslis;

if chnl==0

for i=1:length(idsvec)
    Vds=vdslis(i);
    ids = path_intnp_opt_v3(pab,Vgs,Vds);
    idsvec(i) = ids;
end

elseif chnl==1
    for i=1:length(idsvec)
    Vds=vdslis(i);
    ids = path_intnp_opt_v3_tmd(pab,Vgs,Vds,0);
    idsvec(i) = ids;
    end
else
    error('non-valid chnl signal met.')
end

end