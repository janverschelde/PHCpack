with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with DoblDobl_Complex_Poly_Systems;     use DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_SysFun;      use DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;    use DoblDobl_Complex_Jaco_Matrices;

package DoblDobl_Systems_Pool is

-- DESCRIPTION :
--   Provides storage facility for multiple polynomial systems,
--   with double double precision complex coefficients.
--   For each system we store four different data structures:
--   (1) a system is an array of polynomials, as lists of terms;
--   (2) an evaluator is an array of nester Horner schemes;
--   (3) the Jacobian matrix as matrix of polynomials;
--   (4) the Jacobian evaluator is a matrix of evaluators.
--   This package is modeled after the package Systems_Container.

-- CREATORS :

  procedure Initialize ( n : in integer32 );

  -- DESCRIPTION :
  --   Allocates memory to store n systems in the pool.

  procedure Initialize ( k : in integer32; p : in Poly_Sys );

  -- DESCRIPTION :
  --   Initializes the k-th system in the pool with p,
  --   making a deep copy of the complete data structure for p.

-- CONSTRUCTORS :

  procedure Create ( k : in integer32; p : in Poly_Sys );

  -- DESCRITPION :
  --   Initializes the k-th systems in the pool with p
  --   and creates an evaluator, a Jacobian matrix, and
  --   a Jacobian matrix for the system.

  procedure Create_Evaluator ( k : in integer32 );

  -- DESCRIPTION :
  --   Creates an evaluator for the k-th system in the pool.

  procedure Create_Jacobian_Matrix ( k : in integer32 );

  -- DESCRIPTION :
  --   Creates a Jacobian matrix for the k-th system in the pool.

  procedure Create_Jacobian_Evaluator ( k : in integer32 );

  -- DESCRIPTION :
  --   Creates a Jacobian evaluator for the k-th system in the pool.

-- SELECTORS :

  function Size return natural32;

  -- DESCRIPTION :
  --   Returns the size of the systems pool.

  function Retrieve ( k : integer32 ) return Link_to_Poly_Sys;

  -- DESCRIPTION :
  --   Retrieves the k-th system stored in the pool.

  function Evaluator ( k : integer32 ) return Link_to_Eval_Poly_Sys;

  -- DESCRIPTION :
  --   Retrieves the k-th evaluator stored in the pool.

  function Jacobian_Matrix ( k : integer32 ) return Link_to_Jaco_Mat;

  -- DESCRIPTION :
  --   Retrieves the k-th Jacobian matrix stored in the pool.

  function Jacobian_Evaluator ( k : integer32 ) return Link_to_Eval_Jaco_Mat;

  -- DESCRIPTION :
  --   Retrieves the k-th Jacobian evaluator stored in the pool.

-- DESTRUCTORS :

  procedure Clear_System ( k : in integer32 );

  -- DESCRIPTION :
  --   Deallocation of the k-th system in the pool.

  procedure Clear_Evaluator ( k : in integer32 );

  -- DESCRIPTION :
  --   Deallocation of the k-th evaluator in the pool.

  procedure Clear_Jacobian_Matrix ( k : in integer32 );

  -- DESCRIPTION :
  --   Deallocation of the k-th Jacobian matrix in the pool.

  procedure Clear_Jacobian_Evaluator ( k : in integer32 );

  -- DESCRIPTION :
  --   Deallocation of the k-th Jacobian evaluator in the pool.

  procedure Clear ( k : in integer32 );

  -- DESCRIPTION :
  --   Deallocation of all four data structures for
  --   the k-th system in the pool.

  procedure Clear; 

  -- DESCRIPTION :
  --   Deallocation of all memory for the n systems.

end DoblDobl_Systems_Pool;
