function write_solution(filename,solution)
%
% write_solution(file,s)
%    writes the solution s to a file in PHCpack format
%
% INPUTs:
% filename   file name
% s          solution in structure format
%
id = fopen(filename,'wt');

tempsolution = rmfield(solution,'time');
tempsolution = rmfield(tempsolution,'multiplicity');
tempsolution = rmfield(tempsolution,'err');
tempsolution = rmfield(tempsolution,'rco');
tempsolution = rmfield(tempsolution,'res');

names = fieldnames(tempsolution);
n_var = size(names,1);
tmpsolstr = ['==========================================================================\n'];
format long e;
for k=1:size(solution,2)  
    str_1 = ['solution ' num2str(k) ': \n'];
    str_2 = ['t : ' num2str(real(solution(k).time), '%17.15e')  '  ' num2str(imag(solution(k).time), '%17.15e') '\n'];
    str_3 = ['m : ' num2str(solution(k).multiplicity) '\n'];
    varstr = ['the solution for t : \n'];
    for nv=1:n_var
        temp_2 = [' ' (names{nv}) ' : ' num2str(real(solution(k).(names{nv})), '%17.15e') '  ' num2str(imag(solution(k).(names{nv})), '%17.15e') '\n'];
        varstr = strcat(varstr,temp_2);
    end
    str_4 = ['== err : ' num2str(solution(k).err) ' = rco : ' num2str(solution(k).rco) ' = res : ' num2str(solution(k).res) ' == \n'];
    tmpsolstr = strcat(tmpsolstr,str_1,str_2,str_3,varstr, str_4);
end
solsstr = [num2str(size(solution,2)) ' ' num2str(n_var) '\n' tmpsolstr]; 
fprintf(id, solsstr);
fclose(id);

