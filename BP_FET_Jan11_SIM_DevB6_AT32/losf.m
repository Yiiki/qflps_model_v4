function loss_tot=losf(optvec,args_pts,sigval)
chnl=1;
% sigval=0.1;
vgslis=linspace(-3,3,21);
vgsbch=vgslis(1:1:21);
parmat=parmf(optvec,vgsbch,args_pts,sigval);
loss_vec=(vgsbch.^0-1).*inf;

vdsmin=0;vdsmax=3;vdspts=101;
vds_idx_start=2;vds_spac=9;vds_idx_end=101;

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

vdslis0=linspace(vdsmin,vdsmax,vdspts)';
idsvec0=vdslis0;


if chnl==0

    for ii=vds_idx_start:vds_spac:vds_idx_end
        Vds=vdslis0(ii);
        ids = path_intnp_opt_v3(pab,Vgs,Vds);
        idsvec0(ii) = ids;
    end

elseif chnl==1

    for ii=vds_idx_start:vds_spac:vds_idx_end
        Vds=vdslis0(ii);
        % ids = path_intnp_opt_v3_tmd(pab,Vgs,Vds,0);
        ids = path_intnp_opt_v5_tmd(pab,Vgs,Vds,0);
        idsvec0(ii) = ids;
    end

else
    error('non-valid chnl sig met.')
end

idsvec = idsvec0(vds_idx_start:vds_spac:vds_idx_end);

vgsid=find(vgslis==Vgs);
idsbench0=fet_data(vgsid);

idsben=idsbench0(vds_idx_start:vds_spac:vds_idx_end);

res=(idsben-idsvec).*idsben.^-1;

% loss = length(res).^-0.5.*(res'*res).^0.5;
loss = length(res).^-1.*sum(abs(res));

loss_vec(i)=loss;

end

% loss_tot=sqrt(length(loss_vec).^-1.*sum(loss_vec.^2));
loss_tot=length(loss_vec).^-1.*sum(loss_vec);