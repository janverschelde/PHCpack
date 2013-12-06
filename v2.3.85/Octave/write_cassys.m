function write_cassys(filename,system)
%
% write_cassys(filename,system)
% A utility function special for cascade algorithm.
% The 1st row of written polynomial system only contains number of poly.
% Assuming the system is square.
% INPUT:
% filename =  file name
% system   =  polynomial system in tableau or cell array format
%             see help solve_system for details
%
sysf_npoly = size(system,1);
sysf = system;

id = fopen(filename,'wt');
temps = [num2str(sysf_npoly) '\n']; 
fprintf(id, temps);
clear temps;
for k=1:size(sysf,1)
    fprintf(id, ' ');
    fprintf(id, sysf{k,1});
    fprintf(id, ';\n');
end
fclose(id);
