function [f]=form_poly(n_var,coeff,exponent)
%
% [f]=form_poly(n_var,coeff,exponent)
% form one polynomial given coefficients and exponents list
% INPUTs:
%   n_var  = number of variables in the system
%   coeff  = mx1 vector of coefficients, m=number of terms
%   exponent = nxm matrix of exponents
%              n should be same as n_var 
%              the mth column is the exponents of the mth term
%
% OUTPUT:
%   f  = sum of terms;
%        Each term is a coefficient times variables raised to
%        their respective exponents
% Example:
% n_var    = 4
% coeff    = [1; 1; 1.6009e-001 +9.8710e-001i]
% exponent = [1 0 0;
%             0 1 0;
%             2 0 0;
%             0 2 0]
% function returns a polynomial x1*x3**2 + x2*x4**2 + (0.16009+0.9871*i)
%
r=size(exponent,1);
c=size(exponent,2);
m=size(coeff,1);
if r~=n_var
  error(['number of unknown does NOT match number of exponents']);
end
% form variable list
x=cell(1,n_var);
for ii=1:n_var
  x{1,ii}=['x' num2str(ii)];
end

% looks like matrix multiplication
monomial=cell(1,c);
for ii=1:1
  for jj=1:c
    % ignore exponent is zero and 1
    if(exponent(1,jj)==1)
      temp=[x{ii,1}];
    elseif (exponent(1,jj)==0)
      temp=blanks(0);
    else 
      temp=[x{ii,1} '**' num2str(exponent(1,jj))];
    end
    for kk=2:n_var
      if (exponent(kk,jj)==1)
        if (isempty(temp))
          temp = [x{ii,kk}]; 
        else 
          temp = [temp '*' x{ii,kk}]; 
        end
      elseif (exponent(kk,jj)~=0)
        if (isempty(temp))
          temp = [x{ii,kk} '**' num2str(exponent(kk,jj))];
        else 
          temp = [temp '*' x{ii,kk} '**' num2str(exponent(kk,jj))];
        end
      else 
        temp=[temp];
      end
    end
    if (coeff(jj,1)~=1)
      if (coeff(jj,1)<0)
        if (isempty(temp))
          monomial{ii,jj}=[' ' num2str(real(coeff(jj,1)),'%16.14e')];
        else
          monomial{ii,jj}=[' ' num2str(real(coeff(jj,1)),'%16.14e') '*' temp];
        end
      else % >0 case, coeff. is not possible to be 0
        if (isempty(temp))
          monomial{ii,jj}=[' + ' num2str(real(coeff(jj,1)),'%16.14e')];
        else
          monomial{ii,jj}=[' + ' num2str(real(coeff(jj,1)),'%16.14e') '*' temp];
        end
      end
      if (imag(coeff(jj,1))~=0) % complex coeff.
        if (isempty(temp))
          if(jj~=1)
            if (imag(coeff(jj,1))<0)
              if (abs(imag(coeff(jj,1)))~=1)
                monomial{ii,jj}=[' + (' num2str(real(coeff(jj,1)),'%16.14e') num2str(imag(coeff(jj,1)),'%16.14e') '*i)'];
              else
                monomial{ii,jj}=[' + (' num2str(real(coeff(jj,1)),'%16.14e') '-i)'];
              end
            else
              if (imag(coeff(jj,1))~=1)
                monomial{ii,jj}=[' + (' num2str(real(coeff(jj,1)),'%16.14e') '+' num2str(imag(coeff(jj,1)),'%16.14e') '*i)'];
              else
                monomial{ii,jj}=[' + (' num2str(real(coeff(jj,1)),'%16.14e') '+'  'i)'];
              end
            end    
          else 
            monomial{ii,jj}=['(' num2str(real(coeff(jj,1)),'%16.14e') '+' num2str(imag(coeff(jj,1)),'%16.14e') '*i)'];
          end
        elseif(jj~=1)
          if (imag(coeff(jj,1))<0)
            monomial{ii,jj}=[' + (' num2str(real(coeff(jj,1)),'%16.14e') ' ' num2str(imag(coeff(jj,1)),'%16.14e') '*i)' '*' temp];
          else
            monomial{ii,jj}=[' + (' num2str(real(coeff(jj,1)),'%16.14e') '+' num2str(imag(coeff(jj,1)),'%16.14e') '*i)' '*' temp];
          end
        else 
          monomial{ii,jj}=['(' num2str(real(coeff(jj,1)),'%16.14e') '+' num2str(imag(coeff(jj,1)),'%16.14e') '*i)' '*' temp];
        end
      end
    else % coeff ==1
      if (isempty(temp))
        monomial{ii,jj}=[' + 1'];
      elseif(jj~=1)
        monomial{ii,jj}=[' + ' temp];
      else
        monomial{ii,jj}=[temp];
      end 
    end
  end
end

% sum all terms
sum=monomial{1,1};
for jj=2:c
    sum=[sum monomial{1,jj}];
end
f=[sum];
