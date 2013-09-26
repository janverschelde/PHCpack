function [R]=embed(S,d)
%
% [R]=embed(S,d)
% Constructs an embedded system for a given system S.
% INPUTs:
% S = polynomial system in tableau or cell array format
% note: if S is in cell array format, the polynomials are expanded.
% d = the expected dimension of the solution set
%
% OUPUTs:
% R = an embedded system for S in cell array format
% 
% Example:
% system sys : x^2+y^2+(z-1)^2-1;
%             (z-0.5)*(z-1)*y;
%             (z-0.5)*(z-1)*x;
% i.e. x^2 + y^2 + z^2 - 2z;
%      yz^2 - 1.5yz + 0.5y;
%      xz^2 - 1.5xz + 0.5x;
% sys =  [1    2 0 0;
%        1    0 2 0;
%        1    0 0 2;
%        -2   0 0 1;
%        0    0 0 0;
%        1    0 1 2;
%        -1.5 0 1 1;
%        0.5  0 1 0;
%        0    0 0 0;
%        1    1 0 2;
%        -1.5 1 0 1;
%        0.5  1 0 0;
%        0    0 0 0];
% [embsys]=embed(sys,1)
%
global phcloc
global phctemp
rand('seed',sum(100*clock));
sr = num2str(round(10^8*rand(1,1)));
verbose = 0;
if(verbose)
     fprintf('embedding the given system ... please wait ...\n');
end
sysfile = [phctemp 'sysfile', sr];
embedsysfile = [phctemp 'embedsys', sr];
embedtemp = [phctemp 'embedtemp' sr];
embedphcout = [phctemp 'embedphcout' sr];

% writing start system into a file 
write_system(sysfile, S);

% prepare input redirection file
tempstr = ['1\ny\n' phctemp 'sysfile' sr '\n' phctemp 'embedsys' sr '\n' num2str(d) '\n' 'n\n'];   
embedtemp_id = fopen(embedtemp,'wt');
fprintf(embedtemp_id, tempstr);
fclose(embedtemp_id);

% call phc executable
str = [phcloc ' -c <' embedtemp ' >' embedphcout];
system(str);

% reading the embed system
embedsys_id = fopen(embedsysfile,'rt');
temps = fread(embedsys_id);
fclose(embedsys_id);
temps = char(temps)'; % 1*m string array

% extract the embed system 
[semic] = findstr(temps, ';');
[plus] = findstr(temps, '+');
embed_npoly = size(semic, 2);
R = cell(embed_npoly,1);
R{1} = temps(plus(1):(semic(1)-1));
for k=2:size(R,1)
    R{k} = temps(semic(k-1)+3:semic(k)-1);
end

if(verbose)
   fprintf('%d polynomials in the embedded system.\n', size(R,1));
   for k=1:size(R,1)
       fprintf('%s;\n\n', R{k});
   end
end

% remove intermediate files
delete(sysfile);
delete(embedsysfile);
delete(embedtemp);
delete(embedphcout);

