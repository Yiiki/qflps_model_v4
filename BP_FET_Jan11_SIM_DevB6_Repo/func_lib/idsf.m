function idsm=idsf(optvec,vgsbch,vdslis0,args_pts,sigval)

chnl=1;

vdspts=length(vdslis0);
vgspts=length(vgsbch);
parmat=parmf(optvec,vgsbch,args_pts,sigval);
idsm=zeros(vdspts,vgspts);

for i=1:length(vgsbch)
ue=parmat(1,i);
vthn=parmat(2,i);
uh=parmat(3,i);
vthp=parmat(4,i);
nid=parmat(5,i);
phi=parmat(6,i);
phi0=parmat(7,i);

pab=parLiz_FET(ue,vthn,uh,vthp,nid,phi,phi0);

Vgs=vgsbch(i);

if chnl==0

    for ii=1:vdspts
        Vds=vdslis0(ii);
        ids = path_intnp_opt_v3(pab,Vgs,Vds);
        idsm(ii,i) = ids;
    end

elseif chnl==1

    for ii=1:vdspts
        Vds=vdslis0(ii);
        % ids = path_intnp_opt_v3_tmd(pab,Vgs,Vds,0);
        ids = path_intnp_opt_v5_tmd(pab,Vgs,Vds,0);
        idsm(ii,i) = ids;
    end

else
    error('non-valid chnl sig met.')
end


end

end
