with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Matrices;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Matrices;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Matrices;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Matrices;

package Setup_Flag_Homotopies is

-- DESCRIPTION :
--   This package collects useful functions and procedures to help
--   define homotopies that move flags to solve Schubert problems,
--   in double, double double, and quad double precision.

  function Random_Flag
             ( n : integer32 ) return Standard_Complex_Matrices.Matrix;
  function Random_Flag
             ( n : integer32 ) return DoblDobl_Complex_Matrices.Matrix;
  function Random_Flag
             ( n : integer32 ) return QuadDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a random upper triangular matrix as a
  --   representation of a general flag.

  function One_Flag
             ( n : integer32 ) return Standard_Complex_Matrices.Matrix;
  function One_Flag
             ( n : integer32 ) return DoblDobl_Complex_Matrices.Matrix;
  function One_Flag
             ( n : integer32 ) return QuadDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a flag with the same support as a random flag,
  --   filled with ones instead of random complex numbers.

  function Identity ( n : integer32 ) return Standard_Complex_Matrices.Matrix;
  function Identity ( n : integer32 ) return DoblDobl_Complex_Matrices.Matrix;
  function Identity ( n : integer32 ) return QuadDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the n-by-n identity matrix.

  function Moved_Flag
             ( n : integer32 ) return Standard_Complex_Matrices.Matrix;
  function Moved_Flag
             ( n : integer32 ) return DoblDobl_Complex_Matrices.Matrix;
  function Moved_Flag
             ( n : integer32 ) return QuadDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the coordinates for the moved flag in n-space.

  function Numeric_Transformation
             ( t : Standard_Natural_Matrices.Matrix;
               g : Standard_Complex_Numbers.Complex_Number )
             return Standard_Complex_Matrices.Matrix;
  function Numeric_Transformation
             ( t : Standard_Natural_Matrices.Matrix;
               g : DoblDobl_Complex_Numbers.Complex_Number )
             return DoblDobl_Complex_Matrices.Matrix;
  function Numeric_Transformation
             ( t : Standard_Natural_Matrices.Matrix;
               g : QuadDobl_Complex_Numbers.Complex_Number )
             return QuadDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the transformation t in its numeric form, replacing
  --   the free coefficient by the gamma constant g.  

  function Numeric_Transformation
             ( t : Standard_Natural_Matrices.Matrix )
             return Standard_Complex_Matrices.Matrix;
  function Numeric_Transformation
             ( t : Standard_Natural_Matrices.Matrix )
             return DoblDobl_Complex_Matrices.Matrix;
  function Numeric_Transformation
             ( t : Standard_Natural_Matrices.Matrix )
             return QuadDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the transformation t in its numeric form.
  --   The free coefficients is replaced by one and the coefficient
  --   just below it gets the value minus one.

  procedure Write_Standard_Moving_Flag
              ( file : in file_type;
                flag : in Standard_Complex_Matrices.Matrix );
  procedure Write_DoblDobl_Moving_Flag
              ( file : in file_type;
                flag : in DoblDobl_Complex_Matrices.Matrix );
  procedure Write_QuadDobl_Moving_Flag
              ( file : in file_type;
                flag : in QuadDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Writes the coordinates of the moving flag as an integer matrix.

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

  procedure Insert_Scaling_Symbol ( i,j : in natural32 );

  -- DESCRIPTION :
  --   Inserts the symbol sij to the symbol table,
  --   where sij is used to recondition the swap homotopy.
  --   The last symbol in the symbol table is the continuation parameter.

  function Symbolic_Transformation
             ( n,v : integer32; t : Standard_Natural_Matrices.Matrix )
             return Standard_Complex_Poly_Matrices.Matrix;
  function Symbolic_Transformation
             ( n,v : integer32; t : Standard_Natural_Matrices.Matrix )
             return DoblDobl_Complex_Poly_Matrices.Matrix;
  function Symbolic_Transformation
             ( n,v : integer32; t : Standard_Natural_Matrices.Matrix )
             return QuadDobl_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a symbolic representation of a transformation, as a matrix of
  --   polynomials in n variables, with the pattern given in t.
  --   By default, gamma is one as multiplier to the variable t.

  -- REQUIRED : 1 <= v <= n.

  -- ON ENTRY :
  --   n       total number of variables in the polynomials on return;
  --   v       index of the variable for the continuation parameter t;
  --   t       localization pattern for the transformation as a matrix
  --           of zeroes, ones, and a '2'.  The '2' locates t.

  function Symbolic_Transformation
             ( n,v : integer32;
               gamma : Standard_Complex_Numbers.Complex_Number;
               t : Standard_Natural_Matrices.Matrix )
             return Standard_Complex_Poly_Matrices.Matrix;
  function Symbolic_Transformation
             ( n,v : integer32;
               gamma : DoblDobl_Complex_Numbers.Complex_Number;
               t : Standard_Natural_Matrices.Matrix )
             return DoblDobl_Complex_Poly_Matrices.Matrix;
  function Symbolic_Transformation
             ( n,v : integer32;
               gamma : QuadDobl_Complex_Numbers.Complex_Number;
               t : Standard_Natural_Matrices.Matrix )
             return QuadDobl_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a symbolic representation of a transformation, as a matrix of
  --   polynomials in n variables, with the pattern given in t.

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
             ( n,v : integer32; t : Standard_Natural_Matrices.Matrix )
             return DoblDobl_Complex_Poly_Matrices.Matrix;
  function Inverse_Symbolic_Transformation
             ( n,v : integer32; t : Standard_Natural_Matrices.Matrix )
             return QuadDobl_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a symbolic representation of the inverse of the transformation
  --   defined as a matrix of polynomials in n variables, with the pattern 
  --   given in t, used to define the original symbolic transformation.

  -- REQUIRED : 1 <= v <= n.

  -- ON ENTRY :
  --   n       total number of variables in the polynomials on return;
  --   v       index of the variable for the continuation parameter t;
  --   t       localization pattern for the transformation as a matrix
  --           of zeroes, ones, and a '2'.  The '2' locates t.

  function Inverse_Symbolic_Transformation
             ( n,v : integer32;
               gamma : Standard_Complex_Numbers.Complex_Number;
               t : Standard_Natural_Matrices.Matrix )
             return Standard_Complex_Poly_Matrices.Matrix;
  function Inverse_Symbolic_Transformation
             ( n,v : integer32;
               gamma : DoblDobl_Complex_Numbers.Complex_Number;
               t : Standard_Natural_Matrices.Matrix )
             return DoblDobl_Complex_Poly_Matrices.Matrix;
  function Inverse_Symbolic_Transformation
             ( n,v : integer32;
               gamma : QuadDobl_Complex_Numbers.Complex_Number;
               t : Standard_Natural_Matrices.Matrix )
             return QuadDobl_Complex_Poly_Matrices.Matrix;

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
             ( t : Standard_Complex_Poly_Matrices.Matrix;
               v : Standard_Complex_Numbers.Complex_Number )
             return Standard_Complex_Poly_Matrices.Matrix;
  function Evaluate_Transformation
             ( t : DoblDobl_Complex_Poly_Matrices.Matrix;
               v : DoblDobl_Complex_Numbers.Complex_Number )
             return DoblDobl_Complex_Poly_Matrices.Matrix;
  function Evaluate_Transformation
             ( t : QuadDobl_Complex_Poly_Matrices.Matrix;
               v : QuadDobl_Complex_Numbers.Complex_Number )
             return QuadDobl_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the value of the transformation, substituting the value
  --   for the continuation parameter by v.
  --   The polynomials in the matrix on return have one variable less.

  function Moving_Flag
             ( f : Standard_Complex_Matrices.Matrix;
               t : Standard_Complex_Poly_Matrices.Matrix )
             return Standard_Complex_Poly_Matrices.Matrix;
  function Moving_Flag
             ( f : DoblDobl_Complex_Matrices.Matrix;
               t : DoblDobl_Complex_Poly_Matrices.Matrix )
             return DoblDobl_Complex_Poly_Matrices.Matrix;
  function Moving_Flag
             ( f : QuadDobl_Complex_Matrices.Matrix;
               t : QuadDobl_Complex_Poly_Matrices.Matrix )
             return QuadDobl_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns f*t as a matrix of polynomials, used to give coordinates
  --   for a k-plane in a basis moving to the opposite flag.

  procedure Move_to_Opposite_Flag
             ( f : in out Standard_Complex_Matrices.Matrix;
               q,p : in Standard_Natural_Vectors.Vector;
               mf : in Standard_Complex_Matrices.Matrix );
  procedure Move_to_Opposite_Flag
             ( f : in out DoblDobl_Complex_Matrices.Matrix;
               q,p : in Standard_Natural_Vectors.Vector;
               mf : in DoblDobl_Complex_Matrices.Matrix );
  procedure Move_to_Opposite_Flag
             ( f : in out QuadDobl_Complex_Matrices.Matrix;
               q,p : in Standard_Natural_Vectors.Vector;
               mf : in QuadDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Multiplies f with the numerical transformation matrix
  --   to make one move towards the opposite flag mf, using
  --   q as the parent permutation of p in the specializing poset.

  function Symbolic_Plane
             ( n,k : integer32; p,rows,cols : Standard_Natural_Vectors.Vector )
             return Standard_Complex_Poly_Matrices.Matrix;
  function Symbolic_Plane
             ( n,k : integer32; p,rows,cols : Standard_Natural_Vectors.Vector )
             return DoblDobl_Complex_Poly_Matrices.Matrix;
  function Symbolic_Plane
             ( n,k : integer32; p,rows,cols : Standard_Natural_Vectors.Vector )
             return QuadDobl_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the symbolic representation of a k-plane in n-space
  --   as determined by the localization pattern defined by the permutation p
  --   and the rows and columns of the white checkers.
  --   The matrix on return has n rows and k columns.

  function Filter_Zero_Equations
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function Filter_Zero_Equations
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Filter_Zero_Equations
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   The system on return has the same polynomials as in p,
  --   except for the zero polynomials.

  function Square ( n : integer32;
                    p : Standard_Complex_Poly_Systems.Poly_Sys )
                  return Standard_Complex_Poly_Systems.Poly_Sys;
  function Square ( n : integer32;
                    p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
                  return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Square ( n : integer32;
                    p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
                  return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns an n-by-n system, adding random linear combinations of 
  --   the polynomials p(i), i > n, to the first n polynomials,
  --   in double, double double, or quad double precision.

  function Concatenate
             ( s : Standard_Complex_Poly_Systems.Array_of_Poly_Sys )
             return Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
  function Concatenate
             ( s : DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys )
             return DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
  function Concatenate
             ( s : QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys )
             return QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  -- DESCRIPTION :
  --   Concatenates all systems in s into one system on return.

  procedure Append
              ( s : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                p : in Standard_Complex_Polynomials.Poly );
  procedure Append
              ( s : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                p : in DoblDobl_Complex_Polynomials.Poly );
  procedure Append
              ( s : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                p : in QuadDobl_Complex_Polynomials.Poly );

  -- DESCRIPTION :
  --   Appends the polynomial p to the system s.

end Setup_Flag_Homotopies;
