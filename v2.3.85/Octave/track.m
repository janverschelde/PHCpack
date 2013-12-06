function [L]=track(t,s,startsolution)

%
% [L]=track(t,s,startsolutions)
%    tracks the specified solutions to a given start system 
%    returning the corresponding solutions of a given target system;
%
% INPUTs:
%   t = target system in tableau or cell array format
%   s = start system in tableau or cell array format
%   startsolution = solutions to the start system
%
% OUTPUT:
%   L = structure of solution
%       fields: time,multiplicity,err,rco,res,x1,x2,x3,...
% 
% EXAMPLE: 
% target system: x^2+y^2-1;
%                x^3+y^3-1;
% t = [ 1  2 0;
%       1  0 2;
%       -1 0 0;
%       0  0 0;
%       1  0 3;
%       1  3 0;
%       -1 0 0;
%       0  0 0];
%
% start system: (1+i)*(x^2-1);
%               (1+i)*(y^3-1);
%  s= [(1+i)  2  0;
%     (-1-i)  0  0;
%     0       0  0;
%     (1+i)   0  3;
%     (-1-i)  0  0;
%     0       0  0];
% 
% startsolution = structure of that start system solutions
%                 fields: time, multiplicity, err, rco, res, x1,x2,x3....
%
global phcloc
global phctemp

rand('seed',sum(100*clock));
sr = num2str(round(10^8*rand(1,1)));

verbose = 0;
if(verbose) 
     fprintf('tracking..., please wait...\n');
end  
startfile = [phctemp 'start' sr];
targetfile = [phctemp 'target' sr];
startsolsfile = [phctemp 'startsols' sr];
tracktemp = [phctemp 'tracktemp' sr];
trackphcout = [phctemp 'trackphcout' sr];
outfile = [phctemp 'out' sr];
solsfile = [phctemp 'sols' sr];
phczsolsfile = [phctemp 'phczsols' sr];

L = struct('time',{},'multiplicity',{},'err',{},'rco',{},'res',{}); % empty soltuion structure
% check and reset the time value in the solutions 
flag = 0;
num_sols = size(startsolution,2);
for k=1:num_sols
    if(startsolution(k).time == 1)
        flag = 1;
        break;
    else 
        continue;
    end
end
if(flag) % at least one time field is not start with 0
    fprintf('\nBy default, the continuation parameter goes from time = 0.0 to 1.0 \nPlease reset time fields of the start solutions.\n');
    return;
end

% prepare output string and write it into file
write_solution(startsolsfile,startsolution);

% writing start system into a file 
write_system(startfile,s);

% writing target system into a file 
write_system(targetfile,t);

% tempstr = [targetfile '\n' outfile '\ny\n' solsfile '\n' startfile '\n' startsolsfile '\n0\n0\n0\n'];  
tempstr = [phctemp 'target' sr '\n' phctemp 'out' sr '\ny\n' phctemp 'sols' sr '\n' phctemp 'start' sr '\n' phctemp 'startsols' sr '\n0\n0\n0\n'];   
tracktemp_id = fopen(tracktemp,'wt');
fprintf(tracktemp_id, tempstr);
fclose(tracktemp_id);

% call phc executable 
str = [phcloc ' -p <' tracktemp ' >' trackphcout];
system(str);
str = [phcloc ' -z ' outfile ' ' phczsolsfile];
system(str);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% still need to decide which point to extract from the solutions
% right now, only the final solutions at the end of the path are extracted
phczsols_id = fopen(phczsolsfile,'rt');
phczsols = fread(phczsols_id);
fclose(phczsols_id);
char_phczsols = char(phczsols)'; % 1*m string array

sols_id = fopen(solsfile,'rt');
sols = fread(sols_id);
fclose(sols_id);
char_sols = char(sols)';

% extract individual solution
[L] = extract_sols(char_phczsols);
if(verbose)
    fprintf('Total %d paths.\n', size(startsolution,2));
    fprintf('end points:\n');
    for k=1:size(L,2)
        disp(L(k));
    end
end

% remove all intermediate files
delete(startfile);
delete(targetfile);
delete(startsolsfile);
delete(tracktemp);
delete(trackphcout);
delete(outfile);
delete(solsfile);
delete(phczsolsfile);
