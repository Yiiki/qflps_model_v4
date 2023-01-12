function write_mat_dat(func_path,func_name,mat_data)

fileID=fopen([func_path,func_name,'.m'],'w');
fprintf(fileID,'function y=%s(rows,cols)\n',func_name);
for i=1:size(mat_data,2)
fprintf(fileID,'mat_data(:,%d)=[\n',i);
fprintf(fileID,'%.16e\n',mat_data(:,i));
fprintf(fileID,'];\n');
end
fprintf(fileID,'y=mat_data(rows,cols);\n');
fprintf(fileID,'end');
fclose(fileID);

end