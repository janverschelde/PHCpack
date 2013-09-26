function [W3]=intersection(S1,W1,S2,W2)
%
% [W3]=intersection(S1,W1,S2,W2)
% Applies the diagonal homotopy technique to find a witness set W3 for the
% intersection of the solution sets represented by the witness sets W1 and W2. 
%
% INPUTs:
% S* = polynomial system in cell array format
% W* = witness set of polynomial system in S*
% 
% OUTPUT:
% W3 = a witness set for the intesection 
%      W3{i,1} is (i-1)-dimensional component witness set.
% It means the (i-1)-dimensional component does NOT exist if W3{i,1} is
% empty.
%
global phcloc phctemp
rand('seed',sum(100*clock));
sr = num2str(round(10^8*rand(1,1)));

verbose = 0;
if(verbose)
    fprintf('applying the diagonal homotopy technique ... please wait ...\n');
end
sysfile_1 = [phctemp 'sysfile_1', sr];
sysfile_2 = [phctemp 'sysfile_2', sr];
wsfile_1 = [phctemp 'witnessset_1', sr];
wsfile_2 = [phctemp 'witnessset_2', sr];
intersecttemp = [phctemp 'intersecttemp' sr];
intersectphcout = [phctemp 'intersectphcout' sr];
outfile = [phctemp 'out' sr];
errorfile = [phctemp 'error' sr];

% write system and witness set to files
write_system(sysfile_1, S1);
write_solution(wsfile_1, W1);

write_system(sysfile_2, S2);
write_solution(wsfile_2, W2);

% prepare input redirection file
tempstr = [phctemp 'sysfile_1' sr '\n' phctemp 'witnessset_1' sr '\n' phctemp 'sysfile_2' sr '\n' phctemp 'witnessset_2' sr '\n' phctemp 'out' sr '\n'];   
intersecttemp_id = fopen(intersecttemp,'wt');
fprintf(intersecttemp_id, tempstr);
fclose(intersecttemp_id);

% call phc executable
str = [phcloc ' -w < ' intersecttemp ' > ' intersectphcout];
system(str);

% determine the top dimension of the embeded system
names = fieldnames(W1);
n_field = size(names,1);
n = 0; % top dimension
for j=1:n_field
    if(~isempty(findstr(names{j},'zz')))
        n = n + 1;
    end
end

% for each component, extract solutions into desired format
errordim = cell(n+1,1);
for j=1:n      
    str = [phcloc ' -z ' phctemp 'out' sr '_w' num2str(j-1) ' ' phctemp 'phcz_w' sr num2str(j-1) ' > ' phctemp 'error' sr];
    system(str);
    errorid = fopen(errorfile, 'rt');
    errmsg = char(fread(errorid))';
    fclose(errorid);
    v1 = ~isempty(findstr(errmsg, 'Something wrong with solutions'));
    v2 = ~isempty(findstr(errmsg, 'Could not open the file'));
    if(v1|v2)
      if(verbose) fprintf('No super witness set for dimension %d\n', j-1); end
      errordim{j,1} = j-1;
    else
       if(verbose)
          fprintf('%d dimensional super witness set for intersection exists.\n', j-1);
       end
    end

end
% for each component, read solutions
W3 = cell(0,1); % supter witness sets
n_sw=0;
for j=1:n
    n_sw = n_sw+1;
    flag = 0;
    for k=1:size(errordim,1)
        if((j-1)==errordim{k})
            flag=1;
            break;
         end
    end
    if(~flag) % no solution error for this dimension
       temp = [phctemp 'phcz_w' sr num2str(j-1)];
       solf_id = fopen(temp, 'rt');
       sols = fread(solf_id);
       fclose(solf_id);
       char_sols = char(sols)'; % 1*m string array
       % extract individual solution
       [tempsw] = extract_sols(char_sols);
       W3{n_sw,1} = tempsw;
       if(verbose)
          fprintf('%d dimensional component:\n', j-1);
          W3{n_sw,1}
       end
    else
       continue;
    end 
end
if(isempty(W3))
    if(verbose)
      fprintf('No Supter Witness Sets for intersection exist!\n');
    end
end

% remove intermediate files
delete(sysfile_1);
delete(sysfile_2);
delete(wsfile_1);
delete(wsfile_2);
delete(intersecttemp);
delete(intersectphcout);
delete(outfile);
delete(errorfile);
for j=1:size(W3,1)
    if(isempty(errordim{j,1}))
        str = [outfile '_w' num2str(j-1)];
        delete(str);
    end
    if(~isempty(W3{j,1}))
        str = [phctemp 'phcz_w' sr num2str(j-1)];
        delete(str);
    end
end
if(isempty(W3))
   str = [outfile '_w' num2str(1)];
   delete(str);
end






