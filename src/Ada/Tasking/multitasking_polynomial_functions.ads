with Standard_Integer_Numbers;       use Standard_Integer_Numbers;
with Standard_Natural_VecVecs;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;

package Multitasking_Polynomial_Functions is

-- DESCRIPTION :
--   To evaluate a polynomial system with multiple tasks,
--   we use a flattened representation, via a lexicographically
--   increasingly ordered vector of exponents and a coefficient matrix.

-- EVALUATION OF THE MONOMIALS :

  function Silent_Eval
             ( n : integer32;
               v : Standard_Integer_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector;
  function Silent_Eval
             ( n : integer32;
               v : Standard_Integer_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector;
  function Silent_Eval
             ( n : integer32;
               v : Standard_Integer_VecVecs.VecVec;
               x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Applies n tasks to evaluate the monomial vector defined by v
  --   at the complex vector x.  There is no intermediate output.

  function Reporting_Eval
             ( n : integer32;
               v : Standard_Integer_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector;
  function Reporting_Eval
             ( n : integer32;
               v : Standard_Integer_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector;
  function Reporting_Eval
             ( n : integer32;
               v : Standard_Integer_VecVecs.VecVec;
               x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Applies n tasks to evaluate the monomial vector defined by v at
  --   the complex vector x.  Diagnostics show progress of evaluation.

-- EVALUATION OF THE DENSE FLATTENED SYSTEM :

  function Silent_Eval
             ( n : integer32;
               A : Standard_Complex_Matrices.Matrix;
               v : Standard_Integer_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector;
  function Silent_Eval
             ( n : integer32;
               A : DoblDobl_Complex_Matrices.Matrix;
               v : Standard_Integer_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector;
  function Silent_Eval
             ( n : integer32;
               A : QuadDobl_Complex_Matrices.Matrix;
               v : Standard_Integer_VecVecs.VecVec;
               x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Applies n tasks to evaluate the flattened polynomial system
  --   with coefficient matrix A and monomials defined by v, at x.
  --   There is no intermediate output.

  function Reporting_Eval
             ( n : integer32;
               A : Standard_Complex_Matrices.Matrix;
               v : Standard_Integer_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector;
  function Reporting_Eval
             ( n : integer32;
               A : DoblDobl_Complex_Matrices.Matrix;
               v : Standard_Integer_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector;
  function Reporting_Eval
             ( n : integer32;
               A : QuadDobl_Complex_Matrices.Matrix;
               v : Standard_Integer_VecVecs.VecVec;
               x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Applies n tasks to evaluate the flattened polynomial system
  --   with coefficient matrix A and monomials defined by v, at x.
  --   Intermediate output shows the progress of the evaluation.

-- EVALUATION OF THE SPARSE FLATTENED SYSTEM :

  function Silent_Eval 
             ( n : integer32;
               c : Standard_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector;
  function Silent_Eval 
             ( n : integer32;
               c : DoblDobl_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector;
  function Silent_Eval 
             ( n : integer32;
               c : QuadDobl_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Applies n tasks to evaluate the sparse flattened system at x
  --   with coefficients in c, monomials defined in v and indices in k.
  --   There is no intermediate output.

  function Reporting_Eval 
             ( n : integer32;
               c : Standard_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector;
  function Reporting_Eval 
             ( n : integer32;
               c : DoblDobl_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector;
  function Reporting_Eval 
             ( n : integer32;
               c : QuadDobl_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Applies n tasks to evaluate the sparse flattened system at x
  --   with coefficients in c, monomials defined in v and indices in k.
  --   Intermediate output shows the progress of the evaluation.

-- EVALUATION OF THE SPARSE FLATTENED SYSTEM with LOOPING Workers :

  function Silent_Looping_Eval 
             ( n : integer32;
               c : Standard_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector;
  function Silent_Looping_Eval 
             ( n : integer32;
               c : DoblDobl_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector;
  function Silent_Looping_Eval 
             ( n : integer32;
               c : QuadDobl_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Applies n tasks to evaluate the sparse flattened system at x
  --   with coefficients in c, monomials defined in v and indices in k.
  --   There is no intermediate output.

  function Reporting_Looping_Eval 
             ( n : integer32;
               c : Standard_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector;
  function Reporting_Looping_Eval 
             ( n : integer32;
               c : DoblDobl_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector;
  function Reporting_Looping_Eval 
             ( n : integer32;
               c : QuadDobl_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Applies n tasks to evaluate the sparse flattened system at x
  --   with coefficients in c, monomials defined in v and indices in k.
  --   Intermediate output shows the progress of the evaluation.

end Multitasking_Polynomial_Functions;
