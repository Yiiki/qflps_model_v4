%      6     3     1
%      5     1     1
%      5     2     1
%      5     4     1
%      6     1     1
%      6     2     1
%      6     4     1
%      6     3     2
%      5     1     2
%      5     2     2
%      5     4     2
%      6     1     2
%      6     2     2
%      6     4     2
%      6     3     3
%      5     1     3
%      5     2     3
%      5     4     3
%      6     1     3
%      6     2     3
%      6     4     3
%    Vdd   GND   Vdd (V)
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


for i=1:size(idxSet,1)
    k1=idxSet(i,1);
    k2=idxSet(i,2);
    k3=idxSet(i,3);
    filename=sprintf('1-10-6-mix-inv%d%d_%dV',k1,k2,k3);
    vovgmat=readmatrix(['exp_data',filesep,'INV',filesep,filename,'.xlsx'],'Sheet',1);
    idvgmat=readmatrix(['exp_data',filesep,'INV',filesep,filename,'.xlsx'],'Sheet',2);

    vovgpac{k1,k2,k3}=vovgmat;
    idvgpac{k1,k2,k3}=idvgmat;
end

% newcolor21
% 
% set(groot,'defaultAxesColorOrder',color21,...
%       'defaultAxesLineStyleOrder','-|--|:')

figure
for i=1:size(idxSet,1)

    k1=idxSet(i,1);
    k2=idxSet(i,2);
    k3=idxSet(i,3);
 
    vovgmat = vovgpac{k1,k2,k3};
    vgslis = vovgmat(:,1);
    voumat = vovgmat(:,2:end);
    stag = sprintf('INV#%d-%d@%dV',k1,k2,k3);

    subplot(3,7,i)
    plot(vgslis,voumat)
    ylim([0 k3])
    title(stag)
end

% set(groot,'defaultAxesLineStyleOrder','remove')
% 
% set(groot,'defaultAxesColorOrder','remove')


figure
for i=1:size(idxSet,1)

    k1=idxSet(i,1);
    k2=idxSet(i,2);
    k3=idxSet(i,3);
 
    idvgmat = idvgpac{k1,k2,k3};
    vgslis = idvgmat(:,1);
    idsmat = idvgmat(:,2:end);
    stag = sprintf('INV#%d-%d@%dV',k1,k2,k3);

    subplot(3,7,i)
    plot(vgslis,idsmat)
%     ylim([0 1e-4])
    title(stag)
end
