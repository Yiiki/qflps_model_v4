function write_vec_dat(func_path,func_name,mat_data)

fileID=fopen([func_path,func_name,'.m'],'w');
fprintf(fileID,'function y=%s(idx)\n',func_name);

if size(mat_data,2)>1
    mat_data=mat_data';
end

fprintf(fileID,'mat_data=[\n');
fprintf(fileID,'%.16e\n',mat_data);
fprintf(fileID,'];\n');

fprintf(fileID,'if nargin~=0\n');
fprintf(fileID,'y=mat_data(idx);\n');
fprintf(fileID,'else\n');
fprintf(fileID,'y=mat_data;\n');
fprintf(fileID,'end\n');

fprintf(fileID,'end');
fclose(fileID);

end