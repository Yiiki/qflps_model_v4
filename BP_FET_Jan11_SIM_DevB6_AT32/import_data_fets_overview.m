idxSet=[1,2
        1,3
        1,4
        1,5
        1,6
        2,3
        2,4
        2,5
        2,6
        3,4
        3,5
        3,6
        4,5
        4,6
        5,6];
for i=1:size(idxSet,1)
    k1=idxSet(i,1);
    k2=idxSet(i,2);
    filename=sprintf('U-1-7-8-idvd(0-3)-vbg-3-3(21points)-port%d-%d(%d-gnd-%d-vdd)',k1,k2,k1,k2);
    idvdmat=readmatrix(['exp_data',filesep,'fets',filesep,filename,'.xlsx']);

    idvdpac{idxSet(i,1),idxSet(i,2)}=idvdmat;
end

z1=size(idvdpac,1);
z2=size(idvdpac,2)-1;

newcolor32

set(groot,'defaultAxesColorOrder',color32,...
      'defaultAxesLineStyleOrder','-|--|:')

figure
for i=1:size(idxSet,1)

    k1=idxSet(i,1);
    k2=idxSet(i,2);
 
    k3=(k1-1).*z2 + k2-1;

    idvdmat = idvdpac{k1,k2};
    vgslis = idvdmat(:,1);
    idsmat = idvdmat(:,2:end);
    stag = sprintf('%d-%d',k1,k2);

    subplot(z1,z2,k3)
    plot(vgslis,idsmat)
    ylim([0 8e-6])
    title(stag)
end

set(groot,'defaultAxesLineStyleOrder','remove')

set(groot,'defaultAxesColorOrder','remove')
