idxSet=[
3,1
3,2
3,4
5,3
6,5
];

% for i=1:size(idxSet,1)
addpath(genpath(pwd))
% i=1;
i=3;
    k1=idxSet(i,1);
    k2=idxSet(i,2);
    filename=sprintf('1-10-6-FET-idvd-vbg-3-3(21points)-%d-%d',k1,k2);
    idvdmat=readmatrix(['exp_data',filesep,'fets',filesep,filename,'.xlsx']);

%     idvdpac{idxSet(i,1),idxSet(i,2)}=idvdmat;
% end


vdsmat0 = idvdmat(:,1);
idsmat0 = idvdmat(:,2:end);

vdsfwd = vdsmat0(1:101,:);
idsfwd = idsmat0(1:101,:);

vdsbwd = flipud(vdsmat0(102:end,:));
idsbwd = flipud(idsmat0(102:end,:));

vdsmid = 0.5.*(vdsfwd+vdsbwd);
idsmid = 0.5.*(idsfwd+idsbwd);

% smooth
idsmidx=idsmid;
for i=1:size(idsmid,2)
    if i<15
        sig=0.03;
    else
        sig=0.03;
    end
    idsmidx(:,i)=wtbnd_gaus(vdsmid,idsmid(:,i),vdsmid,sig);
end

%% formal-plot
vspac=1;
% newcolor21
% 
% set(groot,'defaultAxesColorOrder',color21,...
%       'defaultAxesLineStyleOrder','-|--|:')

figure
plot(vdsmat0(1:vspac:202),idsmat0(1:vspac:202,:),'o',vdsmid,idsmid,'-','LineWidth',2)

figure
plot(vdsmid,idsmid,'-',vdsmid,idsmidx,'--','LineWidth',2)

figure
semilogy(vdsmid(2:end,:),idsmid(2:end,:),'-', ...
    vdsmid(2:end,:),idsmidx(2:end,:),'--','LineWidth',2)

figure
semilogy(vdsfwd,idsfwd,'-','LineWidth',2)

figure
semilogy(vdsbwd,idsbwd,'-','LineWidth',2)

figure
semilogy(vdsmid(2:end,1:2:end),idsmidx(2:end,1:2:end),'-','LineWidth',2)

% set(groot,'defaultAxesLineStyleOrder','remove')
% 
% set(groot,'defaultAxesColorOrder','remove')
