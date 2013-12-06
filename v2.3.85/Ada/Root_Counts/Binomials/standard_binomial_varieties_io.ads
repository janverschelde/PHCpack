with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Integer_Matrices;
with Standard_Complex_Laurentials;      use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;

package Standard_Binomial_Varieties_io is

-- DESCRIPTION :
--   This package offers routines to write and read algebraic sets
--   defined by binomial systems.

-- OUTPUT ROUTINES :

  procedure Fill_Symbol_Table ( n : in natural32 );

  -- DESCRIPTION :
  --   Fills the symbol table with n symbols in case it is empty,
  --   n is the dimension of the ambient space. 
  --   The n symbols are x1, x2, etc.

  procedure Write_Header ( file : in file_type; n,d : in natural32 );

  -- DESCRIPTION :
  --   Writes the beginning of a polynomial system of n equations
  --   in n+d variables and the start of the first polynomial to 
  --   initialize the symbol table with the d parameters t1, t2, etc.

  procedure Write_Solution
               ( file : in file_type; d : in natural32;
                 M : in Standard_Integer_Matrices.Matrix;
                 c : in Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the monomial parametrization of the toric binomial variety
  --   to file.

  -- ON ENTRY :
  --   file      file opened for output or standard_output;
  --   d         the dimension of the algebraic set;
  --   M         unimodular coordinate transformation with in its
  --             first d rows the generators of the cone of tropisms;
  --   c         solution to a binomial system in y-coordinates,
  --             where x = y^M.

  procedure Write_System
               ( d : in natural32;
                 M : in Standard_Integer_Matrices.Matrix;
                 c : in Standard_Complex_Vectors.Vector;
                 p : out Standard_Complex_Laur_Systems.Laur_Sys );

  -- DESCRIPTION :
  --   Writes the monomial parametrization of the toric binomial variety  
  --   to a polynomial system of n equations in n+d variables.

  -- REQUIRED : p'range is M'range(1), or equivalently,
  --   p'range = 1..n where n is the number of variables.

  -- ON ENTRY :
  --   file      file opened for output or standard_output;
  --   d         the dimension of the algebraic set;
  --   M         unimodular coordinate transformation with in its
  --             first d rows the generators of the cone of tropisms;
  --   c         solution to a binomial system in y-coordinates,
  --             where x = y^M.
 
  -- ON RETURN :
  --   p         a polynomial system with equations of the form
  --             x(j) - tc(j)*tm, where tm is a monomial in d parameters.
  --             The d parameters come first in the symbol table.

  procedure Write_Free_Affine_Solution
               ( file : in file_type;
                 s,f : in Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes a "free" affine solution which arises after setting all
  --   variables as indicated in s to zero, where the remaining variables
  --   indicated by f are unconstrained because all equations vanished.

  procedure Write_Free_Affine_System
               ( d : in natural32;
                 f : in Standard_Integer_Vectors.Vector;
                 p : out Standard_Complex_Laur_Systems.Laur_Sys );

  -- DESCRIPTION :
  --   Writes a "free" affine solution of dimension d to a polynomial system.

  -- REQUIRED : d equals the number of ones in f and p'range = f'range.

  -- ON ENTRY :
  --   d         the number of free parameters;
  --   f         f(i) = 1 if the i-th variable is free
  --             otherwise the variable will be considered as zero.

  -- ON RETURN :
  --   p         a polynomial system with equations of the form x(j) - t(k),
  --             for free variables, or simply x(j) if zero.
  --             The first d parameters appear first in the symbol table.

  procedure Write_Affine_Solution
               ( file : in file_type; d : in natural32;
                 s,f : in Standard_Integer_Vectors.Vector; 
                 M : in Standard_Integer_Matrices.Matrix;
                 c : in Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the monomial parametrization of the affine binomial variety
  --   to file.  The toric part is given by (d,M,c).
  --   The dimension of the set equals d plus the number of free variables.

  -- REQUIRED : d > 0.

  -- ON ENTRY :
  --   file      file opened for output or standard_output;
  --   d         the dimension of the toric algebraic set;
  --   s         selection of variables to be zero:
  --             s(k) = 1 means that the k-th variable is zero;
  --   f         free variables which dropped out when variables
  --             in s were set to zero;
  --   M         unimodular coordinate transformation with in its
  --             first d rows the generators of the cone of tropisms;
  --   c         solution to a binomial system in y-coordinates,
  --             where x = y^M.

  procedure Write_Affine_System
               ( d,e : in natural32;
                 s,f : in Standard_Integer_Vectors.Vector; 
                 M : in Standard_Integer_Matrices.Matrix;
                 c : in Standard_Complex_Vectors.Vector;
                 p : out Standard_Complex_Laur_Systems.Laur_Sys );

  -- DESCRIPTION :
  --   Writes the monomial parametrization of the affine binomial variety
  --   to a polynomial system.  The toric part is given by (d,M,c).
  --   The dimension of the set equals d plus the number of free variables.

  -- REQUIRED : d > 0.

  -- ON ENTRY :
  --   file      file opened for output or standard_output;
  --   d         the dimension of the toric algebraic set;
  --   e         the number of free variables;
  --   s         selection of variables to be zero:
  --             s(k) = 1 means that the k-th variable is zero;
  --   f         free variables which dropped out when variables
  --             in s were set to zero;
  --   M         unimodular coordinate transformation with in its
  --             first d rows the generators of the cone of tropisms;
  --   c         solution to a binomial system in y-coordinates,
  --             where x = y^M.

  -- ON RETURN :
  --   p         a polynomial system with three types of equations:
  --             p(j) = x(j), if x(j) equals zero;
  --             p(j) = x(j) - t(k), if x(j) is free;
  --             p(j) = x(j) - tc(j)*tm, where tc(j) is the transformed
  --             coefficient and tm a monomial in the parameters.

  procedure Write_Solution
               ( file : in file_type; d : in natural32;
                 M : in Standard_Integer_Matrices.Matrix;
                 w : in Standard_Integer_Vectors.Vector;
                 c : in Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the monomial parametrization f the binomial variety,
  --   where denominators are stored in w.

  -- ON ENTRY :
  --   file      file opened for output or standard_output;
  --   d         the dimension of the algebraic set;
  --   M         denominators of the coordinate transformation with in its
  --             first d rows the generators of the cone of tropisms;
  --   w         denominators for the first d rows of M;
  --   c         solution to a binomial system in y-coordinates,
  --             where x = y^M.

-- INPUT ROUTINES :

  function Variable_Term ( t : Term; n,d : natural32 ) return boolean;

  -- DESCRIPTION :
  --   A variable term has a nonzero entry in t.dg(d+1..d+n).
  --   Returns true if t is a variable term, false if otherwise.

  procedure Extract_Binomial_Variety
               ( p : in Poly; n,d,i : in natural32;
                 T : out Standard_Integer_Matrices.Matrix;
                 c : out Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Extracts the i-th component of the tropisms of the i-th
  --   Laurent polynomial p of a binomial system, and the
  --   corresponding coefficient.
 
  -- REQUIRED : the first d variables in p are parameters and 
  --   p is of the form x - c*t^m, for some coefficient c and
  --   exponent vector m in the d parameters t.
 
  procedure Extract_Binomial_Variety
               ( s : in Laur_Sys; n,d : in natural32;
                 T : out Standard_Integer_Matrices.Matrix;
                 c : out Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Extracts the tropisms of the n binomial equations in s.

  procedure Parse_Binomial_Variety
               ( p : in Laur_Sys;
                 T : out Standard_Integer_Matrices.Link_to_Matrix;
                 c : out Standard_Complex_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Parses a Laurent polynomial system p into a binomial variety.
  
  -- REQUIRED : p consists of n binomial equations in n+d variables
  --   of the form x - c*t^m, where m is defined by the tropisms.
  --   The first d unknowns in every monomial are the t-parameters.

  -- ON ENTRY :
  --   p         a Laurent polynomial system in the correct format,
  --             of n equations in n+d variables, where n is the ambient
  --             dimension and d the dimension of the binomial variety.

  -- ON RETURN :
  --   T         an n-by-d matrix containing the tropisms;
  --   c         a vector of length n with the coefficients.

end Standard_Binomial_Varieties_io;
