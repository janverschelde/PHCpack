function s = cmplx2str(x)

% Writes complex (or real) number x as a string
% 
% cmplx2str(1+2i) = (1.00000000000000e+00+2.00000000000000e+00*i)
% cmplx2str(1-2i) = (1.00000000000000e+00-2.00000000000000e+00*i)
% cmpl2str(3)     = 3.00000000000000e+00
% cmpl2str(-3)    = (-3.00000000000000e+00)
% cmpl2str(2i)    = (2.00000000000000e+00*i)
% cmpl2str(-2i)   = (-2.00000000000000e+00*i)

% Bor Plestenjak 19.5.2014
% Jan Verschelde 12.9.2014, patched bug with negative real coefficients,
%   the -3 became + -3 in the output.  The patch now writes (-3).

if x==0
    s = blanks(0);
    return
end

if real(x)==0
    s = blanks(0);
else
    s = num2str(real(x),'%16.14e');
end

if imag(x)~=0
    if imag(x)>0
        if real(x)~=0
            s = ['(',s,'+',num2str(imag(x),'%16.14e'),'*i)'];
        else
            s = ['(',num2str(imag(x),'%16.14e'),'*i)'];
        end
    else
       s = ['(',s,num2str(imag(x),'%16.14e'),'*i)'];
    end
else % added round brackets for negative real coefficients
   if real(x) < 0
       s = ['(',s,')'];
   end
end
