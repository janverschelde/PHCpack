with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Integer_Matrices;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;

package Standard_Affine_Binomials is

-- DESCRIPTION :
--   This package offers tools to compute all affine solutions of
--   a system of binomial equations.  This computation corresponds to
--   computing a generalized permanent of the incidence matrix of the
--   bipartite graph linking monomials to variables.

  procedure Extract_Two_Terms
              ( p : in Standard_Complex_Laurentials.Poly;
                t1,t2 : out Standard_Complex_Laurentials.Term;
                fail : out boolean );

  -- DESCRIPTION :
  --   Extracts two terms of the binomial p.

  -- ON ENTRY :
  --   p        a Laurent polynomial, assumed to have exactly 2 terms.

  -- ON RETURN :
  --   t1       first term of the binomial, if not fail;
  --   t2       second term of the binomial, if not fail;
  --   fail     is true if p is not a binomials, false otherwise.

  procedure Incidence_Matrix 
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                A : out Standard_Integer_Matrices.Matrix;
                fail : out boolean );

  -- DESCRIPTION :
  --   We construct a bipartite graph relating every variable to
  --   each monomial in which it occurs with a positive exponent.
  --   This procedure creates the incidence matrix of this
  --   bipartite graph: M[x^a,x_k] = 1 if a_k > 0, 0 otherwise.
 
  -- REQUIRED : A'last(2) = Number_of_Unknowns(p(p'first)).
 
  -- ON ENTRY :
  --   p       a system of binomial equations;
 
  -- ON RETURN :
  --   A       the incidence matrix of the bipartite graph,
  --           with rows indexed by the monomials occurring in p
  --           and columns indexed by the variables in p;
  --   fail    true if the system p is not binomial.

  procedure Nonzero_Binomials
             ( A : in Standard_Integer_Matrices.Matrix;
               s : in Standard_Integer_Vectors.Vector;
               e : out Standard_Integer_Vectors.Vector;
               cnt : out integer32; valid : out boolean );

  -- DESCRIPTION :
  --   Given in A the incidence matrix of a binomial system
  --   and in s0 the selection of variables to be set to zero,
  --   the vector on return indicates which binomials remain nonzero.
 
  -- REQUIRED : e'range = 1..A'last(1)/2.
 
  -- ON ENTRY :
  --   A       the incidence matrix of a binomial system;
  --   s       s(i) = 1 if the i-th variable equals zero, 0 otherwise.

  -- ON RETURN :
  --   e       e(i) = 1 if the i-th binomial remains nonzero,
  --           e(i) = 0 if the i-th binomial turns zero by s,
  --           e(i) = -1 if the i-th binomial becomes a monomial,
  --           which means that s is an invalid selection of variables;
  --   cnt     number of remaining binomials (if selection s is valid);
  --   valid   true if the system obtained after setting the variables
  --           to zero as defined by the selection in s is binomial,
  --           false if only one monomial turned zero by s.

  function Free_Variables
             ( A : Standard_Integer_Matrices.Matrix;
               s : Standard_Integer_Vectors.Vector )
             return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   A variable is free if there is no monomial left that contains it.

  -- ON ENTRY :
  --   A       the incidence matrix of a binomial system;
  --   s       selection of variables to be set to zero.

  -- ON RETURN :
  --   a vector of s'range where the i-th component equals 1
  --   if the variable does not occur in a nonzero binomial,
  --   and 0 otherwise.

  function Subsystem
             ( p : Standard_Complex_Laur_Systems.Laur_Sys;
               nq : integer32;
               eq : Standard_Integer_Vectors.Vector ) 
             return Standard_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Returns a selection of nq polynomials in p, defined by eq.

  -- REQUIRED : nq equals the number of ones in eq.

  -- ON ENTRY :
  --   p       a polynomial system;
  --   nq      number of equations to select from p;
  --   eq      if eq(k) = 1, then the k-th equation of p is selected.

  function Eliminate_Variables
             ( x,s : Standard_Integer_Vectors.Vector; s_cnt : integer32 )
             return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Eliminates from the vector x those variables flagged in s with 1.
  --   The number of ones in s equals s_cnt, so the vector on return
  --   has range x'first..x'last - s_cnt.

  function Eliminate_Variables
             ( p : Standard_Complex_Laurentials.Poly;
               s : Standard_Integer_Vectors.Vector; s_cnt : integer32 )
             return Standard_Complex_Laurentials.Poly;

  -- DESCRIPTION :
  --   Returns the polynomial p with all variables indexed by s removed.
  --   The number of variables set to zero equals s_cnt.

  function Eliminate_Variables
             ( p : Standard_Complex_Laur_Systems.Laur_Sys;
               s : Standard_Integer_Vectors.Vector; s_cnt : integer32 )
             return Standard_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Returns the system p with all variables indexed by s removed.
  --   The number of variables set to zero equals s_cnt.

  function Insert_Zero_Values
             ( c : Standard_Complex_Vectors.Vector;
               s : Standard_Integer_Vectors.Vector )
             return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Inserts to the leading coefficient vector as many zeroes as s_cnt,
  --   the number of ones in s, at those k for which s(k) = 1.

  -- REQUIRED : c'range = s'first..s'last - s_cnt.

  function Insert_Zero_Values
             ( c : Standard_Complex_Solutions.Solution;
               s : Standard_Integer_Vectors.Vector; s_cnt : integer32 )
             return Standard_Complex_Solutions.Solution;

  -- DESCRIPTION :
  --   Inserts to the leading coefficient vector as many zeroes as s_cnt
  --   in the locations as specified by the selection s.

  -- REQUIRED : c.v'range = s'first..s'last - s_cnt.
 
  function Insert_Zero_Values
             ( c : Standard_Complex_Solutions.Solution_List;
               s : Standard_Integer_Vectors.Vector; s_cnt : integer32 )
             return Standard_Complex_Solutions.Solution_List;
               
  -- DESCRIPTION :
  --   Inserts to the leading coefficient vectors in c as many zeroes
  --   as s_cnt in locations as specified by the selection s.

  -- REQUIRED : for every solution x, x.v'range = s'first..s'last - s_cnt.

  function Insert_Unit_Vectors
             ( M : Standard_Integer_Matrices.Matrix;
               s : Standard_Integer_Vectors.Vector )
             return Standard_Integer_Matrices.Matrix;

  -- DESCRIPTION :
  --   Inserts standard unit vectors in every row k and column k of M
  --   where s(k) = 1, for as many as s_cnt = number of ones in s.

  -- REQUIRED : M'range(1) = M'range(2) = s'first..s'last - s_cnt.

end Standard_Affine_Binomials;
