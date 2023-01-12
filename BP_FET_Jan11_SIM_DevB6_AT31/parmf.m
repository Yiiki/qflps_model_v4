function parm=parmf(optvec,vgsbch,args_pts,sigval)

% sigval=0.1;

shape_tag=0;

if size(vgsbch,1)>1% row vgsbch preferred
    vgsbch=vgsbch';
    shape_tag=1;
end

if length(args_pts)~=7
    error('invalid args_pts length')
end

if sum(args_pts)~=length(optvec)
    error('invalid optvec met')
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
