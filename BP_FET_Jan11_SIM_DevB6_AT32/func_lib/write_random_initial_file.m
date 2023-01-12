function write_random_initial_file(randmat,filn)
fileID = fopen([filn,'.m'],'w');
fprintf(fileID,'%6s %12s\n','function',['ic = ',filn,'(i)']);

% write x-cor vec
    
for j=1:size(randmat,2)
    fprintf(fileID,['randmat(:,',num2str(j),')=[']);
fprintf(fileID,'%.16e\r\n',randmat(:,j));
fprintf(fileID,'%6s\n','];');
end

% write ux vx wx
fprintf(fileID,'ic=randmat(:,i);\n');

fprintf(fileID,'%6s\n','end');

fclose(fileID);
end