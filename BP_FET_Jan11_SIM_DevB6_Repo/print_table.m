% set_vec [data required]
% set_sig [data required]
vth2Ntrp=2.5145e+12;% V --> cm^-2
nid2phit=0.026;% 1--> V
Vpn=0.2-2;
Vpp=0.19+2;

idxSet=[
    3,1
    3,2
    3,4
    5,3
    6,5
    ];

for idx_FET=1:size(idxSet,1)
% for idx_FET=1

% fprintf('idx_FET=%d\n',idx_FET)

lab=idxSet(idx_FET,:);

k1=idxSet(idx_FET,1);
k2=idxSet(idx_FET,2);

optvec=set_vec{k1,k2};
args_pts=set_cfg{k1,k2};
sigval=set_sig{k1,k2};

cellve=vec2cell(optvec',args_pts);

data2print=zeros(3,8).*NaN;

for i=1:7
    data2print(1:length(cellve{i}),i)=cellve{i};
end

data2print(1,8)=sigval;


dev_name=sprintf('Dev#%d-%d',lab);
tab_tag = {dev_name,'ue','Ntrp_e','uh','Ntrp_h','phi_t','phi','phi0','sig'};

idx=[1 2 3]';
ue=data2print(:,1).*1e4;% cm^2V^-1s^-1
Ntrp_e=(-Vpn+data2print(:,2)).*vth2Ntrp;% cm^-2
uh=data2print(:,3).*1e4;% cm^2V^-1s^-1
Ntrp_h=(-Vpp+data2print(:,4)).*vth2Ntrp;% cm^-2
phi_t=data2print(:,5).*nid2phit;% V
vphi_d=data2print(:,6);% V
vphi_d0=data2print(:,7);% V
sig=data2print(:,8);% V
format_tab=[idx,ue,Ntrp_e,uh,Ntrp_h,phi_t,vphi_d,vphi_d0,sig];

fprintf('%8.12s%12.12s%12.12s%12.12s%12.12s%12.12s%12.12s%12.12s%12.12s\n', ...
    dev_name,'ue','Ntrp_e','uh','Ntrp_h','phi_t','phi','phi0','sig')

fprintf('%s\n','-------------------------------------------------------------------------------------------------------------')

for j=1:3
    fprintf('%8.1d%12.3e%12.3e%12.3e%12.3e%12.3e%12.3e%12.3e%12.3e\n',format_tab(j,:))
end

fprintf('\n')

filled_vec=[ue;Ntrp_e;uh;Ntrp_h;phi_t;vphi_d;vphi_d0]';

filted_vec=filled_vec;

filted_vec(isnan(filled_vec))=[];

set_vec_published{k1,k2}=filted_vec;

end

