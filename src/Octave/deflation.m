function [L]=deflation(S,Sols)
%
% [L]=deflation(S,Sols)
%    deflates the multiplicity of an isolated singular root
%    while reconditioning the problem
%
% INPUTs:
% S = a polynomial system
% Sols = a (possibly partial) list of isolated solutions of S 
% 
% OUTPUT:
% L{1,1} = system at dimension 0
% L{2,1} = solutions at level 0
% L{3,1} = system ...
% L{4,1} = solutions ...
%
%
global phcloc
global phctemp
rand('seed',sum(100*clock));
% sr = num2str(round(10^8*rand(1,1)));
sr = num2str(int32(10^8*rand(1,1)),10);  % the integer format looks nicer

verbose = 0;
if(verbose) 
    fprintf('calling deflation ...\n');
end  

deflationfile = [phctemp 'deflation' sr];
deflationtemp = [phctemp 'deflationtemp' sr];
deflationphcout = [phctemp 'deflationphcout' sr];
solsfile = [phctemp 'sols' sr];
outfile = [phctemp 'out' sr];

% write system into a file
write_system(deflationfile,S);

% write solution into a file
write_solution(solsfile,Sols);

% write input redirection file, navigating the menu of phc -v
% 6 : Newton's method with deflation for isolated singularities,
% 0 : hardware double precision,
% y : read the input system and solutions from file,
% and the last three characters in answering the menu questions:
% 1 : change algorithmic evaluation into symbolic evaluation
% y : confirm the change
% 0 : no more changes to the settings.
temp = ['6\n0\ny\n' phctemp 'deflation' sr '\n' phctemp 'out' sr '\n' phctemp 'sols' sr '\n' '1\ny\n0\n'];
deflationtemp_id = fopen(deflationtemp,'wt');
fprintf(deflationtemp_id, temp);
fclose(deflationtemp_id);

% call phc executable 
str = [phcloc ' -v <' deflationtemp ' >' deflationphcout];
system(str);

% read output files in and extract files genereted by phc
outfile_id = fopen(outfile, 'rt');
temp = char(fread(outfile_id))';
fclose(outfile_id);
seefile = findstr(temp, 'See the file');
timing = findstr(temp, 'TIMING');
n_files = size(seefile,2);
newfiles = cell(n_files,1);
for j=1:(n_files-1)
    newfiles{j,1} = temp(seefile(j)+13:seefile(j+1)-2); % take care of line feed
end
newfiles{n_files,1} = temp(seefile(n_files)+13:(timing(1)-3));

L = cell(n_files*2,1);
for j=1:2
    id = fopen(newfiles{j,1}, 'rt');
    temps = char(fread(id))';
    fclose(id);
    
    % extract the system 
    semic = findstr(temps, ';');
    newlines = findstr(temps, char(10));
    embed_npoly = size(semic, 2);
    R = cell(embed_npoly,1);
    R{1} = temps((newlines(1)+1):(semic(1)-1)); % starting from the second line
    for k=2:size(R,1)
        R{k} = temps(semic(k-1)+3:semic(k)-1);
    end
    L{2*j-1,1} = R;

    % extract the solutions
    str = [phcloc ' -z ' newfiles{j,1} ' ' phctemp sr num2str(j)];
    system(str);
    str = [phctemp sr num2str(j)];
    id = fopen(str, 'rt');
    temps = char(fread(id))';
    fclose(id);
    % modify variable names to satisfy naming convention of field in
    % structure
    lmpos = findstr(temps, 'lm');
    if(~isempty(lmpos))
      for k=1:size(lmpos,2)
         temps(lmpos(k)+2)='_';
         temps(lmpos(k)+4)='_';
         temps(lmpos(k)+6)=' ';
      end
    end
    L{2*j,1} = extract_sols(temps);
end
if(verbose)
   fprintf('Deflation step is done.\n');
end
% remove intermediate files
delete(deflationfile);
delete(deflationtemp);
delete(deflationphcout);
delete(solsfile);
delete(outfile);
for j=1:n_files
    delete(newfiles{j,1});
    str = [phctemp sr num2str(j)];
    delete(str);
end
