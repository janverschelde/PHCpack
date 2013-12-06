function write_system(filename,system)
%
% write_system(filename,system)
%   write a polynomial system to a file.
%
% INPUT:
%   filename   file name
%   system     polynomial system in tableau or cell array format
%              see help solve_system for details
%
formatflag = 0; % default tableau format
if(size(system,2)>1) % system in tableau format
   [sysf,sysf_npoly,sysf_nvar]=make_system(system);
else 
   formatflag = 1; % need to clean file if the system is in cell array format
   sysf_nvar = size(system,1);
   sysf_npoly = sysf_nvar; % for now, assume the system is square
   sysf = system;
end
id = fopen(filename,'wt');
temps = [num2str(sysf_npoly) ' ' num2str(sysf_nvar) '\n']; 
fprintf(id, temps);
clear temps;
for k=1:size(sysf,1)
    fprintf(id, ' ');
    fprintf(id, sysf{k,1});
    fprintf(id, ';\n');
end
fclose(id);
