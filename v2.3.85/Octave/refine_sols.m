function [refsols]=refine_sols(T,sols,tol_err,tol_sing,tol_res,max)
%
% [refsols]=refine_sols(T,sols,tol_err,tol_sing,tol_res,max)
% Refines the given solutions to any specified precision and makes possible to 
% set error tolerance and other parameters in order to fine-tune the solver;
% INPUTs:
% Two format can be used for input
% (1) T   = tableau matrix defining the system, with
%           1) coefficients in first column
%           2) exponents of variables x1,...,xn in columns 2,...,(n+1)
%           3) end of polynomials marked with a row of 0s 
% (2) T = cell array defining the system
%         each cell contains one polynomial
%
% sols        approximate solutions of given system T
% tol_err     tolerance on magnitude of the correction vector;
% tol_sing    threshold on inverse of condition number;
% tol_res     tolerance on residual
% max         maximum number of iterations in Newton's method
%
% OUTPUT:
% refsols     refined approximate solutions.
%
global phcloc
global phctemp
rand('seed',sum(100*clock));
sr = num2str(round(10^8*rand(1,1)));

verbose = 0;
if(verbose)
    fprintf('refine solutions ... please wait ...\n');
end
 refinfile = [phctemp 'refin' sr];
refoutfile = [phctemp 'refout' sr];
refsolfile = [phctemp 'refsol' sr];
solfile = [phctemp 'sol' sr];
refphcout = [phctemp 'refphcout' sr];
reftemp = [phctemp 'reftemp' sr];
write_system(refinfile,T);        % write system into a file
write_solution(solfile,sols);     % write solutions into a file
reftemp_id = fopen(reftemp,'wt');
tempstr1 = num2str(tol_err);
pos = findstr(tempstr1,'e');
tempstr1(pos)='E';
tempstr2 = num2str(tol_res);
pos = findstr(tempstr2,'e');
tempstr2(pos)='E';
tempstr3 = num2str(tol_sing);
pos = findstr(tempstr3,'e');
tempstr3(pos)='E';
% prepare output string and write it into file
tempstr = ['1\ny\n' phctemp 'refin' sr '\n' phctemp 'refout' sr '\n' phctemp 'sol' sr '\n' 'n\ny\n' phctemp 'refsol' sr '\n' '3\n' tempstr1 '\n4\n' tempstr2 '\n5\n' tempstr3 '\n6\n' num2str(max) '\n0\n'];      
fprintf(reftemp_id, tempstr);
fclose(reftemp_id);
% call phc executable to refine solutions
str = [phcloc ' -v <' reftemp ' >' refphcout];
system(str);
str = [phcloc ' -z ' refoutfile ' ' refsolfile];
system(str);
solf_id = fopen(refsolfile, 'rt');  % read solutions in
sols = fread(solf_id);
fclose(solf_id);
char_sols = char(sols)'; % 1*m string array
% extract individual solution
[refsols] = extract_sols(char_sols);
% print out solutions info.
if(verbose)
    fprintf('%d solutions are refined.\n', size(refsols,2));
    for k=1:size(refsols,2)
       disp(refsols(k));
    end
end
% remove all intermediate files
delete(refinfile);
delete(refoutfile);
delete(solfile);
delete(refsolfile);
delete(refphcout);
delete(reftemp);
