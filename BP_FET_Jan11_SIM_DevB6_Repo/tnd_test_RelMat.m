idxSet=[
    3,1
    3,2
    3,4
    5,3
    6,5
    ];

fet_num=size(idxSet,1);
MatchLis=(1:fet_num).^0-1;

figure

for idx_FET = 1:size(idxSet,1)

lab=idxSet(idx_FET,:);

k1=idxSet(idx_FET,1);
k2=idxSet(idx_FET,2);

optvec=set_vec{k1,k2};
args_pts=set_cfg{k1,k2};
sigval=set_sig{k1,k2};

vgslen=21;
vdslen=101;

Vglis=linspace(-3,3,vgslen);
Vdlis=linspace(0,3,vdslen);
Vdlis=Vdlis(11:9:101);
Vdlis(Vdlis==0)=[];% filterd 0
vdslen=length(Vdlis);
idstf=kron(Vdlis',Vglis).^0-1;

for jj=1:vgslen
    Vg=Vglis(jj);
    for ii=1:vdslen
        Vd=Vdlis(ii);
        Vs=0;
        ids=path_intnp_opt_v5_tnd(lab,optvec,args_pts,sigval,Vg,Vd,Vs);
        idstf(ii,jj)=ids;
    end
end

% output
idsbco=ids_output_date{k1,k2};
vgsbch=linspace(-3,3,21);
vdsbch=linspace(0,3,101);


% check idx list
vgs_ck_lis=1:vgslen;
vds_ck_lis=1:vdslen;

for kkk=1:vgslen
    vgs_ck_lis(kkk)=find(vgsbch==Vglis(kkk));
end

for kkk=1:vdslen
    vds_ck_lis(kkk)=find(vdsbch==Vdlis(kkk));
end

idsbco_ck=idsbco(vds_ck_lis,vgs_ck_lis);

RelMat=(idstf-idsbco_ck).*idsbco_ck.^-1;

MeanRelEr=(vgslen.*vdslen).^-1.*sum(sum(abs(RelMat)));
% MeanRelEr=((vgslen.*vdslen).^-1.*sum(sum(abs(RelMat).^2))).^0.5;

% row_er=vdslen.^-0.5.*sum(RelMat.^2,1).^0.5;
% MeanRelEr=sqrt(vgslen.^-1.*sum(row_er.^2));
MatchLis(idx_FET)=100.*(1-MeanRelEr);
MaxMat=1-min(min(abs(RelMat)));
MinMat=1-max(max(abs(RelMat)));
MidMat=1-median(median(abs(RelMat)));

fprintf('Dev#%d-%d MeanRelEr=%.2f%%, Best = %.2f%%, Worst = %.2f%%, Median = %.2f%%\n', ...
    lab,100.*(1-MeanRelEr),1e2*MaxMat,1e2*MinMat,1e2*MidMat)
% 
% 
% 
subplot(2,3,idx_FET)
plot(Vdlis',idsbco_ck,'o',Vdlis',idstf,'-')
xlim([0 3])
ylim([0 1].*1e-4)
yticks(linspace(0,1e-4,6))
fet_tag=sprintf('Dev#%d-%d',lab);
title(fet_tag)

pause(0.1)

dev_ticks{idx_FET}=sprintf('%d-%d',lab);
end


figure
b=bar(1:fet_num,MatchLis);
xlim([0.25 0.75+fet_num])
ylim([0 100])
ylabel('FoM (%)','FontSize',12)
xlabel('Dev#','FontSize',12)
xticks(1:fet_num)
xticklabels(dev_ticks)
xtickangle(0)

a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Arial','fontsize',12)

xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
for j=1:length(labels1)
    labels2{j}=sprintf('%.0f%%',str2num(labels1{j}));
end
text(xtips1,ytips1,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',12)
