% harvest_data.m [required]
% harvest_data_parm_extend.m [required]
% harvest_data_parm.m [not-required]

folderpath='export_txt';
if exist(folderpath,'dir')
    rmdir(folderpath,'s')
end
mkdir(folderpath)

generapath='export_txt_va';
if exist(generapath,'dir')
    rmdir(generapath,'s')
end
mkdir(generapath)

idxSet=[1,2
        2,3
        3,4
        4,5
        5,6
        1,3
        2,4
        3,5
        4,6
        1,4
        2,5
        3,6
        1,5
        2,6
        1,6];

args_tag = {'ue','vthn','uh','vthp','nid','phi','phi0'};
args_len = length(args_tag);

for epoch=1:size(idxSet)

idxnow=idxSet(epoch,:);
k1=idxnow(1);
k2=idxnow(2);
parm_optmZ=set_ptZ{k1,k2};

filepath=sprintf('%s/%d-%d',folderpath,k1,k2);

mkdir(filepath)

%% write section A
filenameA=sprintf('patch_SecA_Dev#%d-%d.txt',idxnow);
fA=['./',filepath,'/',filenameA];

fileID=fopen(fA,'w');

fprintf(fileID,'%s%d%d%s\n','module ambpfet',idxnow,' (V1,V2,V3,gnd);');

fclose(fileID);

%% write section B

filenameB=sprintf('patch_SecB_Dev#%d-%d.txt',idxnow);
fB=['./',filepath,'/',filenameB];

fileID=fopen(fB,'w');

fprintf(fileID,'%s%.2f%s\n','parameter W2L=',wl_set{k1,k2},'*W/L;//W/L');
fprintf(fileID,'%s\n','// fourth rank ----DEVICE-----DESIGN------PARAMETERS---[MODEL-DEPENDENT]-');
fprintf(fileID,'%s%.6e%s\n','parameter sigmaval = ',set_sig{k1,k2},';// sigmaval');

for j=1:args_len
fprintf(fileID,'%s%s%s%d-%d\n\n','// ',args_tag{j},' parameter Dev#',k1,k2);
parm_data=parm_optmZ{j};
parm_cols=size(parm_data,2);
for jj=1:parm_cols
fprintf(fileID,'%s%s%d%s%.6e%s\n','parameter A',args_tag{j},jj,' = ',parm_data(2,jj),';');
end
fprintf(fileID,'\n');
for jj=1:parm_cols
fprintf(fileID,'%s%s%d%s%.6e%s\n','parameter B',args_tag{j},jj,' = ',parm_data(1,jj),';');
end
fprintf(fileID,'\n');
end

%% write section C
filenameC=sprintf('patch_SecC_Dev#%d-%d.txt',idxnow);
fC=['./',filepath,'/',filenameC];

fileID=fopen(fC,'w');
args_pts=set_cfg{k1,k2};

for kk=1:args_len

fprintf(fileID,'\t%s%s%d%s',args_tag{kk},'=modelpt',args_pts(kk),'(');

parm_data=parm_optmZ{kk};
parm_cols=size(parm_data,2);

for pku=1:parm_cols
fprintf(fileID,'%s%s%d%s','A',args_tag{kk},pku,',');
end

for pku=1:parm_cols
fprintf(fileID,'%s%s%d%s','B',args_tag{kk},pku,',');
end

fprintf(fileID,'sigmaval,Vgs);\n');

end

fclose(fileID);

% fprintf(fileID,);

%% concatenate them
source1='./va_sources/fixed_part_A.txt';
source2='./va_sources/fixed_part_B.txt';
source3='./va_sources/fixed_part_C.txt';

filenameD=sprintf('amos_qflps3.0Dev%d%d.va',idxnow);
fD=['./',generapath,'/',filenameD];

copyfile(fA,'f1.txt');
copyfile(source1,'f2.txt');
copyfile(fB,'f3.txt');
copyfile(source2,'f4.txt');
copyfile(fC,'f5.txt');
copyfile(source3,'f6.txt');

system('copy /b f1.txt+f2.txt+f3.txt+f4.txt+f5.txt+f6.txt f7.txt');

copyfile('f7.txt',fD);

delete('f1.txt')
delete('f2.txt')
delete('f3.txt')
delete('f4.txt')
delete('f5.txt')
delete('f6.txt')
delete('f7.txt')

end

fclose('all');