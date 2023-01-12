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

% newcolor21
% 
% set(groot,'defaultAxesColorOrder',color21,...
%       'defaultAxesLineStyleOrder','-|--|:')

figure
for i=1:size(idxSet,1)

    k1=idxSet(i,1);
    k2=idxSet(i,2);
    k3=idxSet(i,3);
 
    vovgmat = vovgpac_intp{k1,k2,k3};
    vgslis = vovgmat(:,1);
    voumat = vovgmat(:,2:end);
    stag = sprintf('INV#%d-%d@%dV',k1,k2,k3);

    % data
    vovgmat0 = vovgpac{k1,k2,k3};
    vgslis0 = vovgmat0(:,1);
    voumat0 = vovgmat0(:,2:end);

    subplot(3,7,i)
    plot(vgslis+0.4,voumat,'-',vgslis0,voumat0,'o')
    xlim([-0.5 3.5])
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
 
    idvgmat = idvgpac_intp{k1,k2,k3};
    vgslis = idvgmat(:,1);
    idsmat = idvgmat(:,2:end);
    stag = sprintf('INV#%d-%d@%dV',k1,k2,k3);

    subplot(3,7,i)
    plot(vgslis,idsmat)
%     ylim([0 1e-4])
    title(stag)
end
