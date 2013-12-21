function [F]=phc_filter(S,HD,LD)
%
% [F]=phc_filter(S,HD,LD)
%    removes the points in a witness set LD that belong 
%    to a higher-dimensional component HD.
%
% INPUTs:
%   S  = embedded system
%   HD = a witness set for a component of a higher dimension 
%   LD = a witness set for a component of a lower dimension 
%   LD, HD are both structure arrays.
%
% OUTPUT:
%   F = the filtered witness set of 0-dimensional witness sets
%
global phcloc phctemp
rand('seed',sum(100*clock));
sr = num2str(round(10^8*rand(1,1)));

verbose = 0;
if(verbose) 
    fprintf('calling homotopy membership test ...\n');
end  

embeddingfile = [phctemp 'embedding' sr];
filtertemp = [phctemp 'filtertemp' sr];
filterphcout = [phctemp 'filterphcout' sr];
errorfile = [phctemp 'error' sr]; 

solsfile = cell(2,1);
for j=1:2
    solsfile{j} = [phctemp 'sols' sr '_' num2str(j)];
end
outfile = [phctemp 'out' sr];
memberfile = [phctemp 'member' sr];

% write system into a file
%write_system(embeddingfile,S);
write_cassys(embeddingfile,S);

% write solution into a file
write_solution(solsfile{1}, (LD));
write_solution(solsfile{2}, (HD));

% write input redirection file
temp = ['1\n' phctemp 'embedding' sr '\n' phctemp 'sols' sr '_' num2str(2) '\n' phctemp 'sols' sr '_' num2str(1) '\n' phctemp 'out' sr '\n']; % something to be done here, what about solsfile_2????
filtertemp_id = fopen(filtertemp,'wt');
fprintf(filtertemp_id, temp);
fclose(filtertemp_id);

% call phc executable 
str = [phcloc ' -f <' filtertemp ' >' filterphcout];
system(str);

str = [phcloc ' -z ' outfile ' ' memberfile ' > ' errorfile]; 
system(str);
% handelling the case that the parameters are in the wrong order.
errorid = fopen(errorfile, 'rt');
errmsg = char(fread(errorid))';
fclose(errorid);
if(findstr('Something wrong with solutions', errmsg))
    fprintf('Please check the input parameters of phc_filter, make sure they are in appropriate order.\n');
    % empty solution structure
    F = struct('time',{}, 'multiplicity',{}, 'err',{}, 'rco',{}, 'res',{}); 
else
    % read solutions 
    solf_id = fopen(memberfile, 'rt');
    sols = fread(solf_id);
    fclose(solf_id);
    char_sols = char(sols)'; % 1*m string array
    % extract individual solution
    [F] = extract_sols(char_sols);
    if(verbose) 
       disp(F);
    end
    if(verbose)
       fprintf('Homotopy membership test is done.\n');
    end
    delete(memberfile);
end

% remove intermediate files
delete(embeddingfile);
delete(filtertemp);
delete(filterphcout);
for j=1:2
   delete(solsfile{j});
end
delete(outfile);
delete(errorfile);

