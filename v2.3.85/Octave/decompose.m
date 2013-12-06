function  [decom] = decompose(S,fws)
%
% [decom] = decompose(S,fws)
% Performing the numeric irreducible decomposition of a pure-dimensional solution 
% component represented by a witness set fws;
% INPUTs:
% S    = embedded system
% fws  = filtered witness set. It is solution set in structure format.
% 
% OUTPUT:
% decom = irreducible decomposition in n*1 cell array format.
%         decom{j,1} = the jth decomposition
%
global phcloc phctemp
rand('seed',sum(100*clock));
sr = num2str(round(10^8*rand(1,1)));

verbose = 0;
if(verbose) 
     fprintf('calling irreducible decomposition ...\n');
end  

embeddingfile = [phctemp 'embedding' sr];
decomposetemp = [phctemp 'decomposetemp' sr];
decomposephcout = [phctemp 'decomposephcout' sr];
solsfile = [phctemp 'sols' sr];
outfile = [phctemp 'out' sr];
irrfile = [phctemp 'irr' sr];

% write system into a file
write_system(embeddingfile,S);

% write solution into a file
write_solution(solsfile,fws);

% write input redirection file
temp = ['2\n' phctemp 'embedding' sr '\n' phctemp 'sols' sr '\n' phctemp 'out' sr '\n' '1\n0\n' ]; 
decomposetemp_id = fopen(decomposetemp,'wt');
fprintf(decomposetemp_id, temp);
fclose(decomposetemp_id);

% call phc executable 
str = [phcloc ' -f <' decomposetemp ' >' decomposephcout];
system(str);

decom = cell(0,1);
n_sols = 0;
n_irr = 1;
while (n_sols < size(fws,2))
   str = [phcloc ' -z ' phctemp 'embedding' sr '_f' num2str(n_irr) ' ' phctemp 'irr' sr '_' num2str(n_irr)];
   system(str);
  
   % read solutions 
   tempstr = [irrfile '_' num2str(n_irr)];
   solf_id = fopen(tempstr, 'rt');
   sols = fread(solf_id);
   fclose(solf_id);
   char_sols = char(sols)'; % 1*m string array
   % extract individual solution
   [decom{n_irr,1}] = extract_sols(char_sols);
   n_sols = n_sols + size(decom{n_irr,1},2); 
   n_irr = n_irr+1;
end
if(verbose)
    fprintf('total %d irreducible factors.\n', n_irr-1);
    for k=1:(n_irr-1)
        fprintf('irreducible factor %d:\n', k); 
        disp(decom{k,1});
    end
end

% remove intermediate files
delete(embeddingfile);
delete(decomposetemp);
delete(decomposephcout);
delete(solsfile);
delete(outfile);
for k=1:(n_irr-1)
      str = [phctemp 'embedding' sr '_f' num2str(k)];
      delete(str);
      str = [phctemp 'irr' sr '_' num2str(k)];
      delete(str);
  end
