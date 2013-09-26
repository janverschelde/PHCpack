function [WS,R]=eqnbyeqn(S)
%
% [WS,R]=eqnbyeqn(S)
% The equation-by-equation solver of PHCpack to find a list witness sets and 
% embedded system for the components of the solution to the original system 
% in every dimension.
% 
% INPUTs:
% S  = polynomial system in tableau or cell array format
%
% OUTPUT:
% WS = witness sets (lower-dimensional and higher-dimensional) in n*1 cell
%      array format. 
%      WS{i,1} is (i-1)-dimensional component witness sets
% It means the (i-1)-dimensional component does NOT exist if WS{i,1} is empty.
% 
% R  = embeded systems in n*1 cell array format.
%      R{i,1} is the corresponding embedded system for dimension (i-1)
% 
global phcloc phctemp
rand('seed',sum(100*clock));
sr = num2str(round(10^8*rand(1,1)));

verbose = 0;
if(verbose) 
    fprintf('calling equation-by-equation solver..., please wait...\n');
end  
sysfile = [phctemp 'sys' sr];
eqnbyeqntemp = [phctemp 'eqnbyeqntemp' sr];
eqnbyeqnphcout = [phctemp 'eqnbyeqnphcout' sr];
outfile = [phctemp 'out' sr];

% make and write system into a file
if (size(S,2)>1) % system in tableau format
    [f,n_poly,n_var]=make_system(S);
    write_system(sysfile,f);
else 
    write_system(sysfile,S);
end

% write input redirection file
temp = ['y\n' phctemp 'sys' sr '\n' phctemp 'out' sr '\nn\n' '0\n' 'n\n' '0\n' '1\n'];
eqnbyeqntemp_id = fopen(eqnbyeqntemp,'wt');
fprintf(eqnbyeqntemp_id, temp);
fclose(eqnbyeqntemp_id);

% call phc executable 
str = [phcloc ' -a <' eqnbyeqntemp ' >' eqnbyeqnphcout];
system(str);

% determine the top dimension 
id = fopen(outfile, 'rt');
temps = char(fread(id))';
fclose(id);
str = [outfile '_w'];
n = size(findstr(temps, str),2);
if(verbose)
   fprintf('The top dimension is %d.\n',n);
end
% for each component, extract solutions into desired format
errordim = cell(n+1,1);
for j=1:(n+1)    
    str = [phcloc ' -z ' phctemp 'out' sr '_w' num2str(j-1)  ' ' phctemp 'sols_w' sr num2str(j-1) ' > ' phctemp 'error' sr];
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
    end    
end
% for each component, read solutions
WS = cell(0,1); % witness sets
R = cell(0,1); % embedded system list
n_ws = 0;
for j=1:(n+1)
   flag = 0;
   n_ws=n_ws+1;
   for k=1:size(errordim,1)
        if((j-1)==errordim{k})
            flag=1;
            break;
         end
   end
   if(~flag) % no solution error for this dimension 
      temp = [phctemp 'sols_w' sr num2str(j-1)];
      solf_id = fopen(temp, 'rt');
      sols = fread(solf_id);
      fclose(solf_id);
      char_sols = char(sols)'; % 1*m string array
      % extract individual solution
      [tempws] = extract_sols(char_sols);
      WS{n_ws,1} = tempws;
      % read embedded system for each dimension in
      str = [phctemp 'out' sr '_w' num2str(j-1)];
      R{n_ws,1} = read_system(str);
      if(verbose)  
         fprintf('%d dimensional component exists.\n', j-1);
         disp(WS{n_ws,1});
         fprintf('embedded system for dimension %d:\n', j-1);
         disp(R{n_ws,1});
      end
   else
         continue;
   end
end

% remove intermediate files
delete(sysfile);
delete(eqnbyeqntemp);
delete(eqnbyeqnphcout);
delete(outfile);
for j=1:(n+1)
    if(isempty(errordim{j,1}))
       str = [outfile '_w' num2str(j-1)];
       delete(str);
       str = [phctemp 'sols_w' sr num2str(j-1)];
       delete(str);
    else
       continue;
    end
end
str = [phctemp 'error' sr];
delete(str);
if(verbose)
   fprintf('Equation-by-equation solving is done.\n');
end
