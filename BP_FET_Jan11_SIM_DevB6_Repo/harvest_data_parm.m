% collect the variable span the folders
% idxSet=[1,2
%         1,3
%         1,4
%         1,5
%         1,6
%         2,3
%         2,4
%         2,5
%         2,6
%         3,4
%         3,5
%         3,6
%         4,5
%         4,6
%         5,6];

idxSet=[
    3,1
    3,2
    3,4
    5,3
    6,5
    ];

args_tag = {'ue','vthn','uh','vthp','nid','phi','phi0','vd0','vs0'};

fet_num=size(idxSet,1);

vgslis_t=linspace(-3,3,201);

parm_eval_mat=zeros(9,201,fet_num);

for i=1:size(idxSet,1)
% for i=1
% for i=[11 6]
fprintf('%d\n',i)

foldname=sprintf('../BP_FET_Jan11_SIM_DevB6_AT%d%d/',idxSet(i,:));
datapath=[foldname,'export_data/cube_parm_pvt.mat'];
funcpath_args=[foldname,'args_config_func.m'];
funcpath_fets=[foldname,'fet_data.m'];
funcpath_plib=[foldname,'parLib.m'];

% copy
copyfile(funcpath_args)
copyfile(funcpath_fets)
copyfile(funcpath_plib)

% bios parameters setting [remote folder data recieved required]
args_cfg=args_config_func(1:4);
args_pts=[2 2 2 2 args_cfg(1:3)];
args_len=sum(args_pts);
sigs_num=1;
iter_num=8;

% data import
optctr0=load(datapath).cube_parm_pvt;
parm_opt=optctr0(:,iter_num,1);% pick up the best
parm_opts=parm_opt(1:args_len);% first-args_len rows for control points
sigvalm=parm_opt(args_len+1:args_len+sigs_num);% last row for sigma

% package setting
parm_optn=parmf(parm_opts,vgslis_t,args_pts,sigvalm);
phid=parm_optn(6,:);
phid0=parm_optn(7,:);
uemat=parm_optn(1,:);
uhmat=parm_optn(3,:);
Lfet=parLib().ch;
vd0=Lfet.^-1.*phid0.*uemat.*exp(-phid.*phid0.^-1);
vs0=Lfet.^-1.*phid0.*uhmat.*exp(-phid.*phid0.^-1);
parm_eval_mat(:,:,i)=[parm_optn;vd0;vs0];

% remove
delete('args_config_func.m')
delete('fet_data.m')
delete('parLib.m')
end

parm_eval_bkc=permute(parm_eval_mat,[2,3,1]);

for j=1:fet_num
label_tag{j}=sprintf('Dev#%d-%d',idxSet(j,:));
end

newcolor32

set(groot,'defaultAxesColorOrder',color32(floor(linspace(1,31,fet_num)),:),...
      'defaultAxesLineStyleOrder','-|--|:')

figure
t = tiledlayout(3,3);
%,'TileSpacing','compact');
for i=1:9
    if i==1 || i==8 || i==9
    nexttile
    semilogy(vgslis_t',parm_eval_bkc(:,:,i),'LineWidth',2)
    title(args_tag{i})
    else
    nexttile
    plot(vgslis_t',parm_eval_bkc(:,:,i),'LineWidth',2)
    title(args_tag{i})
    end
end
% lgd = legend;
% lgd.Layout.Tile = 8;
% lgd.String = label_tag;
% lgd.NumColumns=2;
% lgd.FontName='Arial';
% lgd.FontSize=12;
% set(groot,'defaultAxesLineStyleOrder','remove')
% set(groot,'defaultAxesColorOrder','remove')