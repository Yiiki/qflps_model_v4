% import_data_fets_aimed [required]
vdsmat = vdsmid;
idsmat = zeros(101,21);
for i=2:size(idsmidx,1)
    idsmat(i,:)=10.^interp1(1:size(idsmidx,2),log10(idsmidx(i,:)),linspace(1,size(idsmidx,2),21));
end


filn = 'fet_data';
fileID = fopen([filn,'.m'],'w');
% function idsdata=fet_dat(i)
fprintf(fileID,'%s\n',['function idsdata=',filn,'(i)']);

% idsmat(:,1 -> ...) = [...];
for j=1:size(idsmat,2)
    fprintf(fileID,'%s\r\n',['idsmat(:,',num2str(j),')=[...']);
    fprintf(fileID,'%.16e\r\n',idsmat(:,j));
    fprintf(fileID,'%6s\n','];');
end

fprintf(fileID,'%s\r\n','idsdata=idsmat(:,i);');

% end
fprintf(fileID,'%s\n','end');
fclose(fileID);