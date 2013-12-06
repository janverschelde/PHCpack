function [system]=read_system(filename)
%
% [system]=read_system(filename)
%    Reads polynomial system from a file.
%
% INPUT:
%   filename =  file name
%   system   =  polynomial system in cell array format
%               see help solve_system for details
% EXAMPLE: [T] = read_system('noon3');
%   noon3 is the name of a file which contains the polynomial system.
%
% Note: Please make sure there are no extra blank lines between polynomials.
%

% read system 
verbose = 0;
id = fopen(filename,'rt');
temps = char(fread(id))';
fclose(id);

newline = findstr(temps, char(10));
dim = str2num(temps(1:newline(1)));
if(size(dim,2)>1) % # of poly. and # of var. are read in
    n_poly = dim(1);
else 
    n_poly = dim;
end
semicolon = findstr(temps, ';');
system = cell(n_poly,1);
system{1,1} = temps(newline(1)+1:semicolon(1)-1);
for j=2:n_poly
    system{j,1} = temps(semicolon(j-1)+1:semicolon(j)-1);
end
for j=1:n_poly
    nl = findstr(system{j,1}, char(10));
    system{j,1}(nl) = '';
end
if(verbose)
    fprintf('%d polynomials in the given system.\n', n_poly);
    for k=1:n_poly
        disp(system{k,1});
    end
end
