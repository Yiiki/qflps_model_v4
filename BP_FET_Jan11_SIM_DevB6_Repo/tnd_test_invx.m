% 

Vinlis=linspace(-0.4,3,31)';
vinlen=length(Vinlis);
% lab=[1,6];
% lablis=[(lab(1):lab(2)-1);lab(1)+1:lab(2)]';


idxSet0=[
6,3
5,1
5,2
5,4
6,1
6,2
6,4
];

idxVdd=[
1
2
3
];

idxSet=[kron(ones(3,1),idxSet0),kron(idxVdd,ones(7,1))];

% inv_idx=1;
for inv_idx=1:size(idxSet,1)
k1=idxSet(inv_idx,1);
k2=idxSet(inv_idx,2);
k3=idxSet(inv_idx,3);

invtag=sprintf('%d-%d',k1,k2);
Vdd=k3;
lablis=invtag2lab(invtag);
fetnum=size(lablis,1);
prtnum=fetnum-1;

iddmat=zeros(vinlen,fetnum);
voumat=zeros(vinlen,prtnum);

for i=1:vinlen
    fprintf('i=%d/%d\n',i,vinlen)
    if i==1
        x0=linspace(0,Vdd,prtnum+2);
        x0=x0(2:end-1);
    else
        x0=voumat(i-1,:);
    end
    vin=Vinlis(i);
    fsfun=@(vo) fsfun0(vin,vo,Vdd,lablis,set_vec,set_cfg,set_sig);
    vo=fsolve(fsfun,x0);
    voumat(i,:)=vo;

    idm=fsfun1(vin,vo,Vdd,lablis,set_vec,set_cfg,set_sig);
    iddmat(i,:)=idm;
end

vovgpac_intp{k1,k2,k3}=[Vinlis,voumat];
idvgpac_intp{k1,k2,k3}=[Vinlis,iddmat];

% figure
% plot(Vinlis,voumat,'LineWidth',2)
% xlabel('{\it V}_{in} (V)','FontSize',12)
% ylabel('{\it V}_{ou} (V)','FontSize',12)
% legend({'port-2','port-3','port-4','port-5'},'FontSize',12)
% 
% figure
% plot(Vinlis,iddmat,'LineWidth',2)
% xlabel('{\it V}_{in} (V)','FontSize',12)
% ylabel('{\it I}_{dd} (A)','FontSize',12)
end

figure
for i=1:size(idxSet,1)

    k1=idxSet(i,1);
    k2=idxSet(i,2);
    k3=idxSet(i,3);
 
    vovgmat = vovgpac_intp{k1,k2,k3};
    vgslis = vovgmat(:,1);
    voumat = vovgmat(:,2:end);
    stag = sprintf('INV#%d-%d@%dV',k1,k2,k3);

    subplot(3,7,i)
    plot(vgslis,voumat)
    ylim([0 k3])
    title(stag)
end


function res=fsfun0(vin,vo,Vdd,lablis,set_vec,set_cfg,set_sig)

idm=fsfun1(vin,vo,Vdd,lablis,set_vec,set_cfg,set_sig);% row vo expected

res=idm(2:end)-idm(1);

res=res.*1e9;

end


function idm=fsfun1(vin,vo,Vdd,lablis,set_vec,set_cfg,set_sig)

fetnum=size(lablis,1);

vsvec=[0,vo];% row vo expected
vdvec=[vo,Vdd];

idm=vsvec.^0-1;

for i=1:fetnum

Vg=vin;
Vd=vdvec(i);
Vs=vsvec(i);

k1=lablis(i,1);
k2=lablis(i,2);

lab=[k1,k2];

optvec=set_vec{k1,k2};
args_pts=set_cfg{k1,k2};
sigval=set_sig{k1,k2};

ids=path_intnp_opt_v5_tnd(lab,optvec,args_pts,sigval,Vg,Vd,Vs);

idm(i)=ids;
end

end
