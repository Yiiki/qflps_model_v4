batch_size=5e2;
args_num = 7;
args_pts = [2,2,2,2,3,3,3];
sigval = 0.1;
args_len = sum(args_pts);

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

lov=loup(:,1);
upv=loup(:,2);

% % neat algorithm when args_pts is a scalar
% rand_intm_set=rand(args_num,args_pts,batch_size).*(upv-lov)+lov;
% rand_intv_set=reshape(rand_intm_set,[],batch_size);

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

rand_intv_set=rand(args_len,batch_size).*(upv_here-lov_here)+lov_here;

% partition to 10 fold
fold_num=1e1;
idx_mat=reshape(1:batch_size,[],fold_num);

func_name='rand_data_func';
func_fold='export_rand_func';
func_sub_fold_head='rand_data';

if exist(func_fold,'dir')
    rmdir(func_fold,"s")
end

mkdir(func_fold)

for i=1:fold_num
    func_path=sprintf('.%s%s%s%s%d%s',filesep,func_fold,filesep,func_sub_fold_head,i,filesep);
    mkdir(func_path)
    fileID=fopen([func_path,func_name,'.m'],'w');
    fprintf(fileID,'function y=%s(npidx)\n',func_name);
    for j=1:size(idx_mat,1)
    fprintf(fileID,'rand_data_vec(:,%d)=[\n',j);
    fprintf(fileID,'%.6e\n',rand_intv_set(:,idx_mat(j,i))');
    fprintf(fileID,'];\n');
    end
    fprintf(fileID,'y=rand_data_vec(:,npidx);\n');
    fprintf(fileID,'end');
    fclose(fileID);
end