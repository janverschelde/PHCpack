function set_phcpath(pathstr)
%
% set_phcpath(pathstr) sets the path to the location of 
%                      the executable program phc
%
% INPUT :
%   pathstr   full path name for the file name which defines
%             the location for the executable program phc.
%
%  EXAMPLES :
%    set_phcpath('C:/Downloads/PHCpack/phc')   % on Windows
%    set_phcpath('/tmp/phc')                   % on Linux
%
global phcloc
global phctemp
phcloc = pathstr;
phctemp = './temp/';
if(exist(phctemp, 'dir')==0)
    mkdir(phctemp);
end
