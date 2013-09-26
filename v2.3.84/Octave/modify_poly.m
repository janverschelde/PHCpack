function [mf]=modify_poly(f)
%
% [mf]=modify_poly(f)
% Taking only the real part of the complex coefficients
% of the given polynomial f
%
start = findstr(f,'(');
finish = findstr(f,')');
start_var = finish+1;
finish_var = start-2;
num_coeff = size(start,2);
coeffs = cell(num_coeff,1);
vars = cell(num_coeff-1,1);
for k=1:num_coeff
   temp = f(start(k)+1:finish(k)-1);
   coeffs{k} = real(str2num(temp)); % take the real part
   if (k<=num_coeff-1)
      vars{k} = f(finish(k)+1:start(k+1)-2);
  end
end
% form the new polynomial
mf = blanks(0);
for k=1:num_coeff
    if(k<=num_coeff-1)
       if(coeffs{k}<0)
          mf = [mf num2str(coeffs{k},'%16.14e') vars{k}];
       else % positive coeff.
          mf = [mf '+' num2str(coeffs{k},'%16.14e') vars{k}];
      end
    else
      if(coeffs{k}<0)
          mf = [mf num2str(coeffs{num_coeff},'%16.14e')];
      else
          mf = [mf '+' num2str(coeffs{num_coeff},'%16.14e')];
      end
    end
end



