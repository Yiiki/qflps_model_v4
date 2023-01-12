parpool("local",50)

addpath(genpath(pwd))
args_tag = {'ue','vthn','uh','vthp','nid','phi','phi0'};
samp_num = 5e1;
iter_num = 8;    
args_num = length(args_tag);
args_pts = [2,2,2,2,3,3,3];
sigval = 0.1;
sigs_num = 1;

args_len = sum(args_pts);

file_fold='export_fopen';

if exist(file_fold,'dir')
    rmdir(file_fold,'s')
end

mkdir(file_fold)

filn = 'rand_data_func';
randfh=str2func(filn);


parfor i=1:samp_num

initial_here=feval(randfh,i);

initial_here=round(initial_here,3,"significant");

func1=sprintf('%s%d','func_parm',i);
func2=sprintf('%s%d','func_loss',i);

foldertag=sprintf('./%s/',file_fold);

fprintf('i=%d\n',i)

fstag=sprintf('./%s/log%d.txt',file_fold,i);
fileID=fopen(fstag,'w');
fprintf(fileID,'say hellow\n');

% solution container
parm_opt_mat=zeros(args_len+sigs_num,iter_num);% args
loss_opt_mat=zeros(1,iter_num);% loss

% start work -------------------

for j=1:iter_num

%% rand initial generation
loup = [...
1e-5,   3e-2;  % ue
-1,      6;      % vthn
1e-5,   3e-2; % uh
-1,      6;    % vthp
1e-1,    20;     % nid
-1,   1e1;     % phi
1e-3,   1e1;    % phi0
];

lov = loup(:,1);
upv = loup(:,2);
 
fprintf(fileID,'i=%d/%d, j=%d/%d\n',i,samp_num,j,iter_num);

optf=@(x) losf(x(1:args_len),args_pts,x(args_len+1:args_len+sigs_num));

if j==1
rand_initv_here=[initial_here;sigval];
else
rand_initv_here=parm_opt_mat(:,j-1);
end

lov_here=zeros(args_len,1);
upv_here=zeros(args_len,1);
for ii=1:args_num
if ii==1
    lov_here(1:args_pts(1))=lov(ii);
    upv_here(1:args_pts(1))=upv(ii);
else
    lov_here(sum(args_pts(1:ii-1))+1:sum(args_pts(1:ii)))=lov(ii);
    upv_here(sum(args_pts(1:ii-1))+1:sum(args_pts(1:ii)))=upv(ii);
end
end

lovs=1e-2;
upvs=1e1;

try 
    [parm_opt,fval] = fminsearchcon(optf,rand_initv_here,[lov_here;lovs],[upv_here;upvs]);
catch 
    parm_opt=rand_initv_here;
    fval=inf;
end

parm_opt_mat(:,j)=parm_opt;
loss_opt_mat(:,j)=fval;


% summary of j-th opt
fprintf(fileID,'loss=%.6e\n',fval);

end

write_mat_dat(foldertag,func1,parm_opt_mat);
write_mat_dat(foldertag,func2,loss_opt_mat);

% close work -------------------
fclose(fileID);

end

%% collect data
addpath(genpath(pwd))

cube_parm=zeros(args_len+sigs_num,iter_num,samp_num);
cube_loss=zeros(1,iter_num,samp_num);


out_fold='export_collect';

if exist(out_fold,'dir')
    rmdir(out_fold,'s')
end

mkdir(out_fold)

addpath(genpath(pwd))

for epoch=1:samp_num

func1=sprintf('%s%d','func_parm',epoch);
func2=sprintf('%s%d','func_loss',epoch);

cube_parm(:,:,epoch)=feval(str2func(func1),1:(args_len+sigs_num),1:iter_num);
cube_loss(:,:,epoch)=feval(str2func(func2),1,1:iter_num);

end

% pivot
[min_loss,~]=min(cube_loss(1,:,:),[],2);
[~,sort_idx]=sort(min_loss);% ascent
cube_parm_pvt=cube_parm(:,:,sort_idx);
cube_loss_pvt=cube_loss(:,:,sort_idx);

parm_opt=cube_parm_pvt(:,iter_num,1);% best in this folder
parm_optm=vec2cell(parm_opt(1:args_len)',args_pts);
sigvalm=parm_opt(args_len+1:args_len+sigs_num);
parm_opts=parm_opt(1:args_len);
%% plot
vgslis_t=linspace(-3,3,201);
parm_optn=parmf(parm_opts,vgslis_t,args_pts,sigvalm);

colorPck=newcolor12();
color3=colorPck.c3;
color11=colorPck.c11;
color12=colorPck.c12;

figure
for i=1:7
    subplot(3,4,i)
    plot(linspace(-3,3,args_pts(i)),parm_optm{i},'o', ...
        vgslis_t,parm_optn(i,:),'-')
    title(args_tag{i})
end

% loss
cube_loss_pmt=permute(cube_loss_pvt,[2 3 1]);
loss_best=cube_loss_pmt(:,3:-1:1,1);
subplot(3,4,8)
set(groot,'defaultAxesColorOrder',color3,...
      'defaultAxesLineStyleOrder','-|--|:')
plot(1:iter_num,loss_best,'o-')
title('loss')
set(groot,'defaultAxesLineStyleOrder','remove')
set(groot,'defaultAxesColorOrder','remove')  

% ids

vgslis=linspace(-3,3,21);
vdslis=linspace(0,3,101);
idsmto=idsf(parm_opts,vgslis,vdslis,args_pts,sigvalm);

vgslis_t=linspace(-3,3,201);
vdslis_t=linspace(0.03,3,12);
idsmtt=idsf(parm_opts,vgslis_t,vdslis_t,args_pts,sigvalm);

idsmtb=zeros(101,21);
for ii=1:21
    idsmtb(:,ii)=fet_data(ii);
end

vdslis_idx=1:2:101;
vdslis_idx_log=2:2:101;

set(groot,'defaultAxesColorOrder',color11,...
      'defaultAxesLineStyleOrder','-|--|:')
subplot(3,4,9)
plot(vdslis(vdslis_idx)',idsmto(vdslis_idx,1:2:21),'-', ...
    vdslis(vdslis_idx)',idsmtb(vdslis_idx,1:2:21),'o')


subplot(3,4,10)
plot(vdslis(vdslis_idx)',idsmto(vdslis_idx,1:2:21),'-', ...
    vdslis(vdslis_idx)',idsmtb(vdslis_idx,1:2:21),'o')
ylim([0 1].*1e-5)

subplot(3,4,11)
semilogy(vdslis(vdslis_idx_log)',idsmto(vdslis_idx_log,1:2:21),'-', ...
    vdslis(vdslis_idx_log)',idsmtb(vdslis_idx_log,1:2:21),'o')

set(groot,'defaultAxesLineStyleOrder','remove')
set(groot,'defaultAxesColorOrder','remove')  

set(groot,'defaultAxesColorOrder',color12,...
      'defaultAxesLineStyleOrder','-|o|:')
subplot(3,4,12)
semilogy(vgslis_t,idsmtt,'-', ...
    vgslis,idsmtb(2:9:101,:),'o')
set(groot,'defaultAxesLineStyleOrder','remove')
set(groot,'defaultAxesColorOrder','remove')  

fval=cube_loss_pvt(:,iter_num,1);
sgtitle(['loss = ',num2str(fval),', sigval = ',num2str(sigvalm)])


set(groot,'defaultAxesLineStyleOrder','remove')
set(groot,'defaultAxesColorOrder','remove')  

saveas(gcf,['.',filesep,'export_fig',filesep,'parm_panel_at_',datestr(now,'HH-MM-SS-FFF'),'.fig'])
save('./export_data/cube_parm_pvt.mat','cube_parm_pvt')
save('./export_data/cube_loss_pvt.mat','cube_loss_pvt')
save('./export_data/cube_idvd_pvt.mat','cube_idvd_pvt')
% delete(gcp('nocreate'))

% exit