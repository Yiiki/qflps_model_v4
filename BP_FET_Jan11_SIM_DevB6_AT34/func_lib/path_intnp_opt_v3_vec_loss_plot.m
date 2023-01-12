function loss=path_intnp_opt_v3_vec_loss_plot(optparm,vgsid,vds_idx_start,vds_spac,vds_idx_end, ...
    VgsBgn,VgsSpc,Vgs_1stidx, ...
    vdsmin,vdsmax,vdspts, ...
    chnl)

ue=optparm(1);
vthn=optparm(2);
uh=optparm(3);
vthp=optparm(4);
nid=optparm(5);
phi=optparm(6);
phi0=optparm(7);

pab=parLiz_FET(ue,vthn,uh,vthp,nid,phi,phi0);

Vgs=VgsBgn+VgsSpc.*(vgsid-Vgs_1stidx);

vdslis0=linspace(vdsmin,vdsmax,vdspts)';
idsvec0=vdslis0;

if chnl==0

    for i=vds_idx_start:vds_spac:vds_idx_end
        Vds=vdslis0(i);
        ids = path_intnp_opt_v3(pab,Vgs,Vds);
        idsvec0(i) = ids;
    end

elseif chnl==1

    for i=vds_idx_start:vds_spac:vds_idx_end
    Vds=vdslis0(i);
    ids = path_intnp_opt_v3_tmd(pab,Vgs,Vds,0);
    idsvec0(i) = ids;
    end

else
    error('non-valid chnl sig met.')
end

idsvec = idsvec0(vds_idx_start:vds_spac:vds_idx_end);

idsbench0=fet_data(vgsid);

idsben=idsbench0(vds_idx_start:vds_spac:vds_idx_end);

res=(idsben-idsvec).*idsben.^-1;

loss = vdspts.^-0.5.*(res'*res).^0.5;

figure
yyaxis left
semilogy(vdslis0(vds_idx_start:vds_spac:vds_idx_end),idsben,'o',...
    vdslis0(vds_idx_start:vds_spac:vds_idx_end),idsvec,'--')
yyaxis right
plot(vdslis0(vds_idx_start:vds_spac:vds_idx_end),idsben,'o',...
    vdslis0(vds_idx_start:vds_spac:vds_idx_end),idsvec,'--')

end