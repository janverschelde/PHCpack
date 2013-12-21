function [n_poly,n_var,n_monomials,coeffs,exponents]=extract_info(tableau)
%
% [n_poly,n_var,n_monomials,coeffs,exponents]=extract_info(tableau)
% INPUT:
% tableau   = tableau matrix defining the system, with
%     1) coefficients in first column
%     2) exponents of variables x1,...,xn in columns 2,...,(n+1)
%     3) end of polynomials marked with a row of 0s 
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
%
% OUTPUTs:
% n_poly = number of polynomials in the system
% n_var  = number of variables in the system
% n_monomials = n_monomials(i) is the number of monomials in polynomial i
% coeffs = coefficient list; coeffs(i) is the coefficient of the ith monomial
% exponents = exponents list; ith row is the list of exponents of ith monomial 

% parse tableau into its individual polynomials
n_col = size(tableau,2);
n_var = n_col-1;
% obtain number of polynomials
monos = tableau(:,1)~=0; % zero coeff. in the 1st column indicates the end of poly.
% n_poly = nnz(~monos); % number of 0s in the 1st column would be # of polynomials
% nnz is not supported by Octave
n_poly = 0;
for k=1:size(~monos,1)
    if(~monos(k)~=0)
        n_poly = n_poly+1;
    end
end
% obtain number of monomials
% n_monomials(i) is the number of monomials in polynomial i
end_poly = [0 find(~monos)']; % indices of end of polynomials
n_monomials = end_poly(2:(n_poly+1))-end_poly(1:n_poly)-1; 

% obtain coeff. list and exponents
save = find(monos);
coeffs = tableau(save,1);
exponents = tableau(save,2:n_col)'; % ith column is the list of exponents of ith monomial 

