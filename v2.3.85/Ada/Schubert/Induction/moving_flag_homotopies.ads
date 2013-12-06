with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Vecvecs;
with Standard_Natural_Matrices;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_VecMats;
with Standard_Complex_Poly_Matrices;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Brackets;                          use Brackets;

package Moving_Flag_Homotopies is

-- DESCRIPTION :
--   A moving flag homotopy introduces new coefficients in a general flag
--   as defined by a generalizing sequence in a checker game.

  function Random_Flag
             ( n : integer32 ) return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a random upper triangular matrix as a
  --   representation of a general flag.

  function One_Flag
             ( n : integer32 ) return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a flag with the same support as a random flag,
  --   filled with ones instead of random complex numbers.

  function Identity ( n : integer32 ) return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the n-by-n identity matrix.

  function Numeric_Transformation
             ( t : Standard_Natural_Matrices.Matrix; g : Complex_Number )
             return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the transformation t in its numeric form, replacing
  --   the free coefficient by the gamma constant g.  
  procedure Add_t_Symbol;

  -- DESCRIPTION :
  --   Enlarges the symbol table with t and adds the symbol 't'
  --   as last symbol to represent the continuation parameter.

  procedure Initialize_Homotopy_Symbols
             ( dim : in natural32;
               locmap : in Standard_Natural_Matrices.Matrix );

  -- DESCRIPTION :
  --   Initializes the symbol table to print a homotopy in the moving flag
  --   with the given localization pattern, using t as added symbol.

  function Symbolic_Transformation
             ( n,v : integer32; t : Standard_Natural_Matrices.Matrix )
             return Standard_Complex_Poly_Matrices.Matrix;
  function Symbolic_Transformation
             ( n,v : integer32; gamma : Complex_Number;
               t : Standard_Natural_Matrices.Matrix )
             return Standard_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a symbolic representation of a transformation, as a matrix of
  --   polynomials in n variables, with the pattern given in t.
  --   By default, gamma is one as multiplier to the variable t.

  -- REQUIRED : 1 <= v <= n.

  -- ON ENTRY :
  --   n       total number of variables in the polynomials on return;
  --   v       index of the variable for the continuation parameter t;
  --   gamma   multiplier for the continuation parameter t;
  --   t       localization pattern for the transformation as a matrix
  --           of zeroes, ones, and a '2'.  The '2' locates t.

  function Inverse_Symbolic_Transformation
             ( n,v : integer32; t : Standard_Natural_Matrices.Matrix )
             return Standard_Complex_Poly_Matrices.Matrix;
  function Inverse_Symbolic_Transformation
             ( n,v : integer32; gamma : Complex_Number;
               t : Standard_Natural_Matrices.Matrix )
             return Standard_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a symbolic representation of the inverse of the transformation
  --   defined as a matrix of polynomials in n variables, with the pattern 
  --   given in t, used to define the original symbolic transformation.
  --   By default, gamma is one as multiplier to the variable t.

  -- REQUIRED : 1 <= v <= n.

  -- ON ENTRY :
  --   n       total number of variables in the polynomials on return;
  --   v       index of the variable for the continuation parameter t;
  --   gamma   multiplier for the continuation parameter t;
  --   t       localization pattern for the transformation as a matrix
  --           of zeroes, ones, and a '2'.  The '2' locates t.

  function Evaluate_Transformation
             ( t : Standard_Complex_Poly_Matrices.Matrix; v : Complex_Number )
             return Standard_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the value of the transformation, substituting the value
  --   for the continuation parameter by v.
  --   The polynomials in the matrix on return have one variable less.

  function Moving_Flag
             ( f : Standard_Complex_Matrices.Matrix;
               t : Standard_Complex_Poly_Matrices.Matrix )
             return Standard_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns f*t as a matrix of polynomials, used to give coordinates
  --   for a k-plane in a basis moving to the opposite flag.

  procedure Move_to_Opposite_Flag
             ( f : in out Standard_Complex_Matrices.Matrix;
               q,p : in Standard_Natural_Vectors.Vector;
               mf : in Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Multiplies f with the numerical transformation matrix
  --   to make one move towards the opposite flag mf, using
  --   q as the parent permutation of p in the specializing poset.

  function Symbolic_Plane
             ( n,k : integer32; p,rows,cols : Standard_Natural_Vectors.Vector )
             return Standard_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the symbolic representation of a k-plane in n-space
  --   as determined by the localization pattern defined by the permutation p
  --   and the rows and columns of the white checkers.
  --   The matrix on return has n rows and k columns.

  function Filter_Zero_Equations ( p : Poly_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   The system on return has the same polynomials as in p,
  --   except for the zero polynomials.

  function Square ( n : integer32; p : Poly_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   Returns an n-by-n system, adding random linear combinations of 
  --   the polynomials p(i), i > n, to the first n polynomials.

  procedure One_Flag_Homotopy
             ( n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Bracket;
               vf,mf,nf : in Standard_Complex_Matrices.Matrix;
               h : out Link_to_Poly_Sys );
  procedure One_Flag_Homotopy
             ( file : in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Bracket;
               vf,mf,nf : in Standard_Complex_Matrices.Matrix;
               h : out Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Builds the polynomial equations in the homotopy with one flag.

  -- ON ENTRY :
  --   file    for intermediate output and diagnostics;
  --   n       dimension of the ambient space, number of black checkers;
  --   k       dimension of the plane, number of white checkers;
  --   q       parent permutation of p in the specializing poset;
  --   p       permutation of the first n numbers;
  --   rows    row positions of the white checkers;
  --   cols    column position of the white checkers;
  --   ic      intersection condition for the fixed flag;
  --   vf      coordinates of a general flag to keep fixed;
  --   mf      coordinates of the moving flag;
  --   nf      accumulation of numeric form of moving flag,
  --           starts at the identity matrix.

  -- ON RETURN :
  --   h       homotopy to move one flag towards the final mf.

  function Concatenate ( s : Array_of_Poly_Sys ) return Link_to_Poly_Sys;

  -- DESCRIPTION :
  --   Concatenates all systems in s into one system on return.

  procedure Flag_Conditions
             ( n,k : in integer32;
               x : in Standard_Complex_Poly_Matrices.Matrix;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               f : out Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Returns in f the system of polynomial equations to express the
  --   conditions inposed by the input planes on a k-plane in n-space x.

  -- ON ENTRY :
  --   n       dimension of the ambient space, number of black checkers;
  --   k       dimension of the plane, number of white checkers;
  --   x       represents a moving k-plane in n-space;
  --   ic      intersection conditions for the fixed flags;
  --   vf      coordinates of general flags to keep fixed.

  -- ON RETURN :
  --   f       conditions of the flag on a k-plane in n-space
  --           meeting the input planes in vf.

  procedure Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               f : out Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Returns in f the system of polynomial equations to express the
  --   conditions inposed by the input planes on a k-plane in n-space,
  --   for a degenerate moving flag, equal to the identity matrix.

  -- ON ENTRY :
  --   n       dimension of the ambient space, number of black checkers;
  --   k       dimension of the plane, number of white checkers;
  --   p       permutation of n numbers in the specializing poset;
  --   rows    row positions of the white checkers;
  --   cols    column position of the white checkers;
  --   ic      intersection conditions for the fixed flags;
  --   vf      coordinates of general flags to keep fixed.

  -- ON RETURN :
  --   f       conditions of the flag on a k-plane in n-space
  --           meeting the input planes in vf.

  procedure Flag_Conditions
             ( file : in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               mf,nf : in Standard_Complex_Matrices.Matrix;
               f : out Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Returns in f the system of polynomial equations to express the
  --   conditions inposed by the input planes on a k-plane in n-space,
  --   for a flag intermediate between the identity and opposite flag.

  -- ON ENTRY :
  --   file    for diagnostic output;
  --   n       dimension of the ambient space, number of black checkers;
  --   k       dimension of the plane, number of white checkers;
  --   q       parent permutation of p in the specializing poset;
  --   p       permutation of the first n numbers;
  --   rows    row positions of the white checkers;
  --   cols    column position of the white checkers;
  --   ic      intersection conditions for the fixed flags;
  --   vf      coordinates of general flags to keep fixed;
  --   mf      coordinates of the opposite flag;
  --   nf      current accumulated flag, moving to the opposite flag.

  -- ON RETURN :
  --   f       conditions of the flag on a k-plane in n-space
  --           meeting the input planes in vf.

  procedure Moving_Flag_Homotopy
             ( n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               mf,nf : in Standard_Complex_Matrices.Matrix;
               h : out Link_to_Poly_Sys );
  procedure Moving_Flag_Homotopy
             ( file : in file_type; n,k : in integer32;
               q,p,rows,cols : in Standard_Natural_Vectors.Vector;
               ic : in Standard_Natural_VecVecs.VecVec;
               vf : in Standard_Complex_VecMats.VecMat;
               mf,nf : in Standard_Complex_Matrices.Matrix;
               h : out Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Builds the polynomial equations in the homotopy for a moving flag,
  --   capable to keep track of multiple intersection conditions.

  -- ON ENTRY :
  --   file    for intermediate output and diagnostics;
  --   n       dimension of the ambient space, number of black checkers;
  --   k       dimension of the plane, number of white checkers;
  --   q       parent permutation of p in the specializing poset;
  --   p       permutation of the first n numbers;
  --   rows    row positions of the white checkers;
  --   cols    column position of the white checkers;
  --   ic      intersection conditions for the fixed flags;
  --   vf      coordinates of general flags to keep fixed;
  --   mf      coordinates of the moving flag;
  --   nf      accumulation of numeric form of moving flag,
  --           starts at the identity matrix.

  -- ON RETURN :
  --   h       homotopy to move one flag towards the final mf.

  function Inconsistent ( p : Poly_Sys ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the system contains at least one nonzero
  --   polynomial of degree zero.

  function Linear_Equations ( p : Poly_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   Returns only those polynomials in p that are linear.
  --   The polynomials on return are shared pointers, not deep copies.

  procedure Coefficients
             ( p : in Poly_Sys; 
               A : out Standard_Complex_Matrices.Matrix;
               b : out Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Given a linear system encoded in the polynomials in p,
  --   returns the coefficient matrix in A and the righthandside vector
  --   of the linear system in b.

  -- REQUIRED : A'range(1) = p'range = b'range
  --   and A'last(2) equals the number of variables in p.

  procedure First_Solution
             ( f : in Poly_Sys; fail : out boolean;
               x : out Standard_Complex_Vectors.Vector;
               res : out double_float );

  -- DESCRIPTION :
  --   Computes the first solution for flag conditions with an
  --   initial localization pattern that turn into a linear system.

  -- ON ENTRY :
  --   f       flag conditions represented as a polynomial system.

  -- ON RETURN :
  --   fail    true if no solution could be found, false otherwise;
  --   x       start solution for h if not fail;
  --   res     residual of the start solution.

  procedure Start_Solution
             ( h : in Poly_Sys; fail : out boolean;
               x : out Standard_Complex_Vectors.Vector;
               res : out double_float );

  -- DESCRIPTION :
  --   Returns a start solution to the homotopy or reports failure.

  -- REQUIRED : x'range = 1..n, n+1 = number of variables in h.

  -- ON ENTRY :
  --   h       a moving flag homotopy, t is the last variable.

  -- ON RETURN :
  --   fail    true if no solution could be found, false otherwise;
  --   x       start solution for h if not fail;
  --   res     residual of the start solution.

  function Cheater_Homotopy_Flag 
             ( nv : in integer32;
               start,target : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns start*(1-t) + target*t for continuation parameter t
  --   as a matrix of nv+1 unknowns, where t is at the last position.

  procedure Flag_Conditions
             ( n,k : in integer32;
               p,rows,cols : in Standard_Natural_Vectors.Vector;
               cond : in Bracket;
               start,target : in Standard_Complex_Matrices.Matrix;
               f : out Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Returns in fthe flag conditions for a k-plane in n-space with
  --   localization pattern determined by (p,rows,cols) that has
  --   to meet the Cheater homotopy flag under the conditions cond.

end Moving_Flag_Homotopies;
