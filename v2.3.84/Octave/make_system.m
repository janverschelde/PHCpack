function [f,n_poly,n_var]=make_system(tableau)
%
% [f,n_poly,n_var]=make_system(tableau)
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
% OUTPUTs: 
% f      = polynomial system 
% n_poly = number of polynomials in the system
% n_var  = number of variables in the system
%

[n_poly,n_var,n_monomials,coeffs,exponents]=extract_info(tableau);

% if n_poly~=n_var
%     error(['number of unknowns must be equal to number of polynomials.']);
%     error(['number of unknowns: ', num2str(n_var), '\n']);
%     error(['number of polynomials: ', num2str(n_poly), '\n']);
% end

% form polynomial system
f=cell(n_poly,1);
j=0;
for k=1:n_poly
  i=j+1;
  j=j+n_monomials(k);
  [f{k,1}]=form_poly(n_var,coeffs(i:j),exponents(:,i:j));
end

