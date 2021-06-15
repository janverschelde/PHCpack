with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_Poly_Matrices;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;

package Standard_Multiple_Deflation is

-- DESCRIPTION :
--   The operations in this package enable higher-order deflation,
--   also known as multiple deflation because its close link to
--   the determination of the multiplity structure.

  function Symbolic_Deflate
              ( nq,nv : natural32;
                a : Standard_Complex_Poly_Matrices.Matrix;
                h : Standard_Complex_Matrices.Matrix ) return Poly_Sys;

  -- DESCRIPTION :
  --   Given the nullity matrix a and random matrix h, the symbolic
  --   deflation of the polynomial system in the first n rows of the
  --   first column of a will be returned.
  --   The number of multipliers equals a'last(2)-1.

  -- REQUIRED : a'last(2) = h'last(2).

  -- ON ENTRY :
  --   nq       number of equations in the original system;
  --   nv       number of variables in the original system;
  --   a        nullity matrix, with in its first n rows of its first column
  --            the original polynomial system;
  --   h        random matrix with in its rows the coefficients of linear
  --            equations in the multipliers.

  function Symbolic_Deflate
              ( nq,nv,r : natural32; 
                a : Standard_Complex_Poly_Matrices.Matrix ) return Poly_Sys;

  -- DESCRIPTION :
  --   Given the nullity matrix a of a system of nq equations in nv unknowns,
  --   and given the corank r, a deflated system is returned.

  function Symbolic_Deflate
              ( p : Poly_Sys; d,r : natural32 ) return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the symbolic deflation of the system p of order d,
  --   where r is the corank of the nullity matrix at the zero.

  -- ON ENTRY :
  --   p        a polynomial system;
  --   d        order of the deflation;
  --   r        corank of the nullity matrix at the zero.

  procedure Predict_Order
              ( p : in Eval_Poly_sys;
                A : in out Standard_Complex_Matrices.Matrix;
                z : in Vector; tol : in double_float; max_d : in integer32;
                corank : out natural32; d : out natural32 );
  procedure Predict_Order
              ( file : in file_type; p : in Eval_Poly_Sys;
                A : in out Standard_Complex_Matrices.Matrix;
                z : in Vector; tol : in double_float; max_d : in integer32;
                corank : out natural32; d : out natural32 );
  procedure Predict_Order
              ( p : in Poly_Sys; z : in Vector; tol : in double_float;
                corank : out natural32; d : out natural32 );
  procedure Predict_Order
              ( file : in file_type;
                p : in Poly_Sys; z : in Vector; tol : in double_float;
                corank : out natural32; d : out natural32 );

  -- DESCRIPTION :
  --   Predicts the order of the deflation using vectors from the
  --   kernel of the Jacobian matrix at the vector z.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   p        a polynomial system;
  --   A        Jacobian matrix of p evaluated at z;
  --   z        an approximate zero;
  --   tol      tolerance to determine the numerical rank;
  --   max_d    maximal bound on the order.

  -- ON RETURN :
  --   corank   corank of the Jacobian matrix at z;
  --   d        prediction for the order of the deflation.

  procedure Numeric_Deflate
              ( p : in Poly_Sys; z : in Vector;
                d : in natural32; tol : in double_float;
                r : out natural32; dp : out Link_to_Poly_Sys );
  procedure Numeric_Deflate
              ( file : in file_type; p : in Poly_Sys; z : in Vector;
                d : in natural32; tol : in double_float;
                r : out natural32; dp : out Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Uses the corank of the nullity matrix of order d of p at z
  --   to set up the augmented polynomial system to deflate p at z.

  -- ON ENTRY : 
  --   file     for intermediate output and diagnostics;
  --   p        a polynomial system;
  --   z        an approximate root of the system;
  --   d        order of the deflation;
  --   tol      tolerance to determine the numerical rank.

  -- ON RETURN :
  --   r        corank of the nullity matrix of order d of p at z;
  --   dp       deflated polynomial system.

  procedure Interactive_Symbolic_Deflation
              ( file : in file_type; p : in Poly_Sys;
                sol : in out Solution; tol : in double_float );
  procedure Interactive_Symbolic_Deflation
              ( file : in file_type; p : in Poly_Sys;
                sols : in out Solution_List; tol : in double_float );

  -- DESCRIPTION :
  --   Interactive implementation of Newton's method with multiple deflation.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   p        a polynomial system;
  --   sol(s)   approximate solutions to p;
  --   tol      tolerance to decide rank.

  -- ON RETURN :
  --   sol(s)   refined solutions.

end Standard_Multiple_Deflation;
