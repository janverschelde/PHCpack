function [solutions]=solve_system(T)
%
% [sols,n_var]=solve_system(T)
%    calls the blackbox solver of phc, returns approximations to all
%    complex isolated roots of a square system
%
% INPUT:
% Two format can be used for input
% (1) T   = tableau matrix defining the system, with
%           1) coefficients in first column
%           2) exponents of variables x1,...,xn in columns 2,...,(n+1)
%           3) end of polynomials marked with a row of 0s 
% (2) T = cell array defining the system
%         each cell contains one polynomial
% Example(gaukwa2):
%  w1       + w2       + (-9.98250904334731E-01 + 5.91196413630250E-02*i);
%  w1*x1    + w2*x2    + (-8.92749639148806E-01 + 4.50553084330444E-01*i);
%  w1*x1**2 + w2*x2**2 + ( 1.60088552022675E-01 + 9.87102657027770E-01*i);
%  w1*x1**3 + w2*x2**3 + (-7.25369971319578E-01 + 6.88359211972815E-01*i);
%  tableau format:
%              1                                               1 0 0 0
%              1                                               0 1 0 0
%              -9.98250904334731E-01 + 5.91196413630250E-02*i  0 0 0 0
%              0                                               0 0 0 0
%              1                                               1 0 1 0
%              1                                               0 1 0 1
%              -8.92749639148806E-01 + 4.50553084330444E-01*i  0 0 0 0
%              0                                               0 0 0 0
%              1                                               1 0 2 0
%              1                                               0 1 0 2
%              1.60088552022675E-01 + 9.87102657027770E-01*i   0 0 0 0
%              0                                               0 0 0 0
%              1                                               1 0 3 0
%              1                                               0 1 0 3
%              -7.25369971319578E-01 + 6.88359211972815E-01*i  0 0 0 0
%              0                                               0 0 0 0
% cell array format:
%   T = cell(4,1);
%   T{1} = 'w1       + w2       + (-9.98250904334731E-01 + 5.91196413630250E-02*i)';
%   T{2} = 'w1*x1    + w2*x2    + (-8.92749639148806E-01 + 4.50553084330444E-01*i)';
%   T{3} = 'w1*x1**2 + w2*x2**2 + ( 1.60088552022675E-01 + 9.87102657027770E-01*i)';
%   T{4} = 'w1*x1**3 + w2*x2**3 + (-7.25369971319578E-01 + 6.88359211972815E-01*i)';
%
% OUPUTS:
% solutions = structure of solution
%             fields: time,multiplicity,err,rco,res,x1,x2,x3,...
%
global phcloc % location of the phc executable
global phctemp 
rand('seed',sum(100*clock));
sr = num2str(round(10^8*rand(1,1)));

verbose=0;
if(verbose) 
     fprintf('solving polynomial system ... please wait...\n');
end
infile = [phctemp 'in' sr];
outfile = [phctemp 'out' sr];
solfile = [phctemp 'sol' sr];
% write system
write_system(infile, T);

% call black box solver
str = [phcloc ' -b ' infile ' ' outfile];
system(str);

str = [phcloc ' -z ' infile ' ' solfile];
system(str);

% read solutions in
solf_id = fopen(solfile, 'rt');
sols = fread(solf_id);
fclose(solf_id);
char_sols = char(sols)'; % 1*m string array
% extract individual solution
[solutions] = extract_sols(char_sols);
% print out solutions info.
if(verbose) 
    fprintf('Blackbox solver found %d solutions.\n', size(solutions,2));
    for k=1:size(solutions,2)
       disp(solutions(k));
    end
end

clear sols char_sols n_var
% remove all intermediate files
delete(infile);
delete(outfile);
delete(solfile); 
