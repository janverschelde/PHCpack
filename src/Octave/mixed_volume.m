function [vol,startsys,sols] = mixed_volume(s,flag_1,flag_2)

% [vol, start, sols] = mixed_volume(s, flag_1, flag_2)
% computer mixed volume for polynomial system s
% INPUT:
% s    = polynomial system
% flag_1 = 1 ---- return a random coefficient start system and solutions
%          0 ---- default. Do not return start system and solutions
%          
% flag_2 = 1 ---- default. return stable mixed volume
%        = 0 ---- return common mixed volume
%
% OUTPUT:
% vol   = stable mixed volume by default. If flag_2 is set to be 0, it is
%         common mixed volume.
% startsys = a random coefficient start system
% sols  = solutions fof the start system
% 
% Example: [vol,startsys,sols] = mixed_volume(s,1,1)
% compute stable mixed volume for system s and 
% return a random coefficient start system and solutions.
% 
% NOTE: 
% No startsys and sols are returned if flag is set to be 0.
% sols could be empty in the case that mixed volume of given system s is 0.
% 
global phcloc
global phctemp

set_default_options = false;
if(nargin==1)
    set_default_options = true;
end
if(nargin==2) 
    flag_2=1; %stable mixed volume will be returned.
end
if(set_default_options)
    flag_1=1;
    flag_2=1;
end
if(flag_2==1) 
    mixedvol_str = 'stable mixed volume :';
else mixedvol_str = 'common mixed volume :';
end

rand('seed',sum(100*clock));
sr = num2str(round(10^8*rand(1,1)));

sysfile = [phctemp 'mixedvolsysfile', sr];
outfile = [phctemp 'mixedvolout', sr];
solsfile = [phctemp 'mixedsols', sr];
startfile = [phctemp 'mixedstart', sr];
temp = [phctemp 'mixedvoltemp' sr];
phcout = [phctemp 'mixedvolphcout' sr];

% writing poly. system into a file 
write_system(sysfile, s);

% prepare input redirection file
tempstr = ['y\n' sysfile '\n' outfile '\n' '4\n' 'y\n' 'n\n' 'y\n' startfile '\n'  '0\n0\n'];   
temp_id = fopen(temp,'wt');
fprintf(temp_id, tempstr);
fclose(temp_id);

% call phc executable
str = [phcloc ' -m <' temp ' >' phcout];
system(str);

% extract the mixed volume, start system and solutions
outfile_id = fopen(outfile,'rt');
temps = fread(outfile_id);
fclose(outfile_id);
temps = char(temps)'; % 1*m string array
volstart = findstr(temps, mixedvol_str);
if(flag_2==0)
    volend = findstr(temps, 'stable mixed volume :');
else
    volend = findstr(temps, 'total mixed volume :');
end
vol = str2num(temps(volstart+length(mixedvol_str):volend-1));

startsys = cell(0,0);
sols = struct('time',{},'multiplicity',{},'err',{},'rco',{},'res',{});
if(flag_1==1)
     startsys = read_system(startfile);
     if(vol~=0)
        str = [phcloc ' -z ' outfile ' ' solsfile];
        system(str);
        % extract solutions
        solsfile_id = fopen(solsfile, 'rt');
        temps = fread(solsfile_id);
        fclose(solsfile_id);
        temps = char(temps)';
        sols = extract_sols(temps);
      end
end

% remove intermediate files
delete(sysfile);
delete(outfile);
if(vol~=0 & flag_1)
  delete(solsfile);
end
delete(startfile);
delete(temp);
delete(phcout);

