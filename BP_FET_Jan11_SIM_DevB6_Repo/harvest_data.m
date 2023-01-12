idxSet=[
    3,1
    3,2
    3,4
    5,3
    6,5
    ];

fet_num=size(idxSet,1);

for i=1:size(idxSet,1)
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
parm_optm=vec2cell(parm_opts',args_pts);% converted to cell array

% package setting

set_vec{idxSet(i,1),idxSet(i,2)}=parm_opts';
set_cfg{idxSet(i,1),idxSet(i,2)}=args_pts;
set_sig{idxSet(i,1),idxSet(i,2)}=sigvalm;
set_pts{idxSet(i,1),idxSet(i,2)}=parm_optm;


% calculate current

% test cond
vgslis=linspace(-3,3,21);
vdslis=linspace(0,3,101);

vdslis_len=length(vdslis);
vgslis_len=length(vgslis);

fet_data_mat=fet_data(1:vgslis_len);

idsbco=fet_data_mat(1:vdslis_len,1:2:vgslis_len);
idsbct=fet_data_mat(11:9:vdslis_len,1:1:vgslis_len)';

% test data package

ids_output_dat{idxSet(i,1),idxSet(i,2)}=idsbco;
ids_transf_dat{idxSet(i,1),idxSet(i,2)}=idsbct;

% test data package (no deleted)

idsbcoe=fet_data_mat(1:vdslis_len,1:vgslis_len);
ids_output_date{idxSet(i,1),idxSet(i,2)}=idsbcoe;

% output curves setting
output_vds=linspace(0,3,101)';
output_vgs=linspace(-3,3,11);

output_vds_len=length(output_vds);
output_vgs_len=length(output_vgs);

% simulated output data and package

output_ids=idsf(parm_opts,output_vgs,output_vds,args_pts,sigvalm);

ids_output_sim{idxSet(i,1),idxSet(i,2)}=output_ids;

% transf curves setting
transf_vds=linspace(0.3,3,11);
transf_vgs=linspace(-3,3,201)';

transf_vds_len=length(transf_vds);
transf_vgs_len=length(transf_vgs);

% simulated transf data and package

transf_ids=idsf(parm_opts,transf_vgs',transf_vds,args_pts,sigvalm)';

ids_transf_sim{idxSet(i,1),idxSet(i,2)}=transf_ids;

figure
subplot(2,1,1)
plot(output_vds,idsbco,'o',output_vds,output_ids,'-')
ylim([0 1].*1e-4)

subplot(2,1,2)
semilogy(vgslis',idsbct,'o',transf_vgs,transf_ids,'-')

tag_AT=sprintf('%d-%d',idxSet(i,:));

sgtitle(tag_AT)

pause(0.1)
% remove
delete('args_config_func.m')
delete('fet_data.m')
delete('parLib.m')
end