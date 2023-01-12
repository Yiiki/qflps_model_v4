% harvest_data.m [required]
% or external_data_source.mat [required]

z1=5;z2=5;


%% output curves

newcolor32

set(groot,'defaultAxesColorOrder',color32(1:3:32,:),...
      'defaultAxesLineStyleOrder','-|--|:')

figure
for i=1:size(idxSet,1)

    k1=idxSet(i,1);
    k2=idxSet(i,2);
 
    k3=(k1-1).*z2 + k2-1;

    stag = sprintf('%d-%d',k1,k2);

    idsbco=ids_output_dat{k1,k2};
    output_ids=ids_output_sim{k1,k2};

    subplot(z1,z2,k3)
    plot(output_vds,idsbco,'o',output_vds,output_ids,'-')
    ylim([0 1].*1e-5)
    title(stag)
end

set(groot,'defaultAxesLineStyleOrder','remove')

set(groot,'defaultAxesColorOrder','remove')


%% transfer curves

newcolor32

set(groot,'defaultAxesColorOrder',color32(1:3:32,:),...
      'defaultAxesLineStyleOrder','-|--|:')

figure
for i=1:size(idxSet,1)

    k1=idxSet(i,1);
    k2=idxSet(i,2);
 
    k3=(k1-1).*z2 + k2-1;

    stag = sprintf('%d-%d',k1,k2);

    idsbct=ids_transf_dat{k1,k2};
    transf_ids=ids_transf_sim{k1,k2};

    subplot(z1,z2,k3)
    semilogy(vgslis',idsbct,'o',transf_vgs,transf_ids,'-')
    title(stag)
end

set(groot,'defaultAxesLineStyleOrder','remove')

set(groot,'defaultAxesColorOrder','remove')

