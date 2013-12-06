function [solutions] = extract_sols(char_sols)
%
% [solutions] = extract_sols(char_sols)
% extract individual solution info. from black-box solver result 
% INPUT: 
% char_sols  = solution string returned by black-box solver
%
% OUTPUT:
% solutions = structure of solution
%             fields: time,multiplicity,err,rco,res,x1,x2,x3,...
%

% define structure to hold solutions
% other fields like x1,x2, etc. will be added dynamically
solutions = struct('time',0, 'multiplicity',0, 'err',0, 'rco',0, 'res',0); 
[start] = findstr(char_sols, '[');
[finish] = findstr(char_sols, '],');
[bracket_comma] = findstr(char_sols, '],');
% replacing I with i
[capi_start] = findstr(char_sols,'I');
for k=1:size(capi_start,2)
    char_sols(capi_start(k))='i';
 end
%n_sols = size(start,2)-1; % number of solutions
n_sols = size(bracket_comma,2)+1;
if(isempty(bracket_comma))
    % only one solution
    begin = start(2)+1;
    ed = findstr(char_sols, ']];')-1;
    solution{1} = char_sols(begin:ed);
else
    % 1st solution
    begin = start(2)+1;
    ed = finish(1)-1; 
    solution{1} = char_sols(begin:ed);
    for k=2:n_sols 
      begin = bracket_comma(k-1)+5;
      if(k~=n_sols)
         ed = bracket_comma(k)-1; 
      else
         ed = findstr(char_sols, ']];')-1;  
      end
      solution{k} = char_sols(begin:ed);
   end   
end

clear start finish begin ed;
% extract field names
[comma] = findstr(solution{1}, ',');
[newline] = findstr(solution{1}, char(10));
[equalsign] = findstr(solution{1},'=');    
n_field = size(newline,2)+3;
tempstr{1,1} = solution{1}(1:(equalsign(1)-2));
for j=2:n_field
    if(j==(n_field-1) || j==n_field) % 3 spaces in front of "res" and "rco"
                      % in the solution string
        start = comma(size(comma,2)-(n_field-j))+3;       
    else 
        start = newline(j-1)+2;
    end
    stop = equalsign(j)-2;
    if(solution{1}(start)==' ')
        start = start+1;
    end
    if(solution{1}(stop)==' ')
        stop = stop-1;
    end
    tempstr{1,j} = solution{1}(start:stop);    
end
    
% send solution information to solution structure
for k=1:n_sols
    % extract individual solution 
    [newline] = findstr(solution{k}, char(10));
    [comma] = findstr(solution{k}, ',');
    [equalsign] = findstr(solution{k},'=');
    field = cell(1,n_field);  
    for j=1:(n_field-1)
        if(j==1)
            temps = solution{k}(1:comma(j));
            field{1,j} = findstr(temps, tempstr{1,j});
        else 
            temps = solution{k}(comma(j-1):comma(j)); % only using part of the solution string
            field{1,j} = findstr(temps, tempstr{1,j});
        end
        start = field{1,j}+ size(tempstr{1,j},2) + 3; 
        stop = size(temps,2)-1;
        solutions(k).(tempstr{1,j}) = str2num((temps(start:stop)));  
    end
    
    % extract "res" value
    field{1,n_field} = findstr(solution{k}, tempstr{1,n_field});
    start = field{1,n_field}+ size(tempstr{1,n_field},2) + 4; 
    stop = size(solution{k},2);
    solutions(k).(tempstr{1,n_field}) = str2num(solution{k}(start:stop));
end
clear solution comma start_x1 start_x2 start_x3 start_x4 start_time start_mult start_err start_rco start_res;


