function [SW,R] = cascade(S,L)
%
% [SW,R] = cascade(S,L)
% Launching the so-called cascade of homotopies for an embedded system S to
% find a list of witness sets and embedded system for the components of the 
% solution set L in every dimension;
% 
% INPUTs:
% S  = embedded system in cell array format
% L  = structure of solutions for the embedded system S
%
% OUTPUT:
% SW = witness sets (lower-dimensional and higher-dimensional) in n*1 cell
%      array format. 
%      SW{1,1} is 0-dimensional component witness sets
%      SW{2,1} is 1-dimensional component witness sets
%      SW{i,1} is (i-1)-dimensional component witness sets
% It means the (i-1)-dimensional component does NOT exist if SW{i,1} is empty.
%
% R  = embeded systems in n*1 cell array format.
%      R{i,1} is the corresponding embedded system for dimension (i-1)

global phcloc
global phctemp
rand('seed',sum(100*clock));
sr = num2str(round(10^8*rand(1,1)));

verbose = 0;
if(verbose) 
    fprintf('runs the cascade of homotopies..., please wait...\n');
end  
embeddingfile = [phctemp 'embedding' sr];
cascadetemp = [phctemp 'cascadetemp' sr];
cascadephcout = [phctemp 'cascadephcout' sr];
solsfile = [phctemp 'sols' sr];
outfile = [phctemp 'out' sr];

% write system into a file
write_cassys(embeddingfile,S);

% write solution into a file
write_solution(solsfile,L);

% write input redirection file
temp = ['2\n' phctemp 'embedding' sr '\n' phctemp 'sols' sr '\n' phctemp 'out' sr '\n' '0\n'];
cascadetemp_id = fopen(cascadetemp,'wt');
fprintf(cascadetemp_id, temp);
fclose(cascadetemp_id);

% call phc executable 
str = [phcloc ' -c <' cascadetemp ' >' cascadephcout];
system(str);

% determine the top dimension of the embeded system
names = fieldnames(L);
n_field = size(names,1);
n = 0; % top dimension
for j=1:n_field
    if(~isempty(findstr(names{j},'zz')))
        n = n + 1;
    end
end

% for each component, extract solutions into desired format
errordim = cell(n+1,1);
for j=1:(n+1)      
    str = [phcloc ' -z ' phctemp 'embedding' sr '_sw' num2str(j-1) ' ' phctemp 'sw' sr num2str(j-1) ' > ' phctemp 'error' sr];
    system(str);
    str = [phctemp 'error' sr]; 
    errorid = fopen(str, 'rt');
    errmsg = char(fread(errorid))';
    fclose(errorid);
    v1 = ~isempty(findstr(errmsg, 'Something wrong with solutions'));
    v2 = ~isempty(findstr(errmsg, 'Could not open the file'));
    if(v1|v2)
      if(verbose) 
          fprintf('No super witness set for dimension %d.\n', j-1); 
      end
      errordim{j,1} = j-1;
    else
      if(verbose)
         fprintf('super witness set for dimension %d exists.\n', j-1);
      end
    end

end
% for each component, read solutions
SW = cell(0,1); % supter witness sets
R = cell(0,1); % embedded system list
n_sw=0;
for j=1:(n+1)
    flag = 0;
    n_sw=n_sw+1;
    for k=1:size(errordim,1)
        if((j-1)==errordim{k})
            flag=1;
            break;
         end
    end
    if(~flag) % no solution error for this dimension
       temp = [phctemp 'sw' sr num2str(j-1)];
       solf_id = fopen(temp, 'rt');
       sols = fread(solf_id);
       fclose(solf_id);
       char_sols = char(sols)'; % 1*m string array
       % extract individual solution
       [tempsw] = extract_sols(char_sols);
       SW{n_sw,1} = tempsw;
       % read embedded system for each dimension in
       str = [phctemp 'embedding' sr '_sw' num2str(j-1)];
       R{n_sw,1} = read_system(str);
       if(verbose) 
          fprintf('%d dimensional component:\n', j-1);
          disp(SW{n_sw,1});
          fprintf('embedded system for dimension %d:\n', j-1);
          disp(R{n_sw,1});
       end
    else
       continue;
    end 
end

% remove intermediate files
delete(embeddingfile);
delete(cascadetemp);
delete(cascadephcout);
delete(solsfile);
delete(outfile);
for j=1:(n+1)
     if(isempty(errordim{j,1}))
         str = [embeddingfile '_sw' num2str(j-1)];
         delete(str);
         str = [phctemp 'sw' sr num2str(j-1)];
         delete(str);
     else
        continue;
    end
end
str = [phctemp 'error' sr];
delete(str);
if(verbose)
   fprintf('Cascade of homotopies is done.\n');
end