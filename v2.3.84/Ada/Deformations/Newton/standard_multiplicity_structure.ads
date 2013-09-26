with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
 
package Standard_Multiplicity_Structure is

-- DESCRIPTION :
--   This package calculates the dual space at an isolated zero.
--   The dimension of the dual space is the multiplicity.

  procedure Orthogonal_Basis
              ( nr,nc : in natural32; a : in Matrix; b : out Matrix );

  -- DESCRIPTION :
  --   Computes an orthogonal basis for the column space of a,
  --   which has nr rows and nc columns.  The result is in b.

  procedure Multiplicity_Structure
             ( f : in Poly_Sys; z : in Standard_Complex_Vectors.Vector;
               tol : in double_float; max_h : in natural32;
               h : out Standard_Natural_Vectors.Vector; m : out natural32 );
  procedure Multiplicity_Structure
             ( file : in file_type;
               f : in Poly_Sys; z : in Standard_Complex_Vectors.Vector;
               tol : in double_float; max_h : in natural32;
               h : out Standard_Natural_Vectors.Vector; m : out natural32 );

  -- DESCRIPTION :
  --   Computes the multiplicity structure of the zero z of f(x) = 0.

  -- REQUIRED : h'range = 0..max_h.

  -- ON ENTRY :
  --   file    for writing diagnostics;
  --   f       a polynomial system;
  --   z       an isolated singular solution of f(x) = 0;
  --   tol     tolerance used to determine the numerical rank;
  --   max_h   upper bound on the order of the differentials,
  --           used to form the nullity matrices.

  -- ON RETURN :  
  --   h       vector of range 0..max_h with values of the Hilbert function,
  --           h(i) gives increase in nullity from i-th order differentials;
  --   m       multiplicity of zero z, equals the sum of the elements in h.

  procedure Multiplicity_Structure
             ( file : in file_type; output : in boolean;
               f : in Poly_Sys; sols : in out Solution_List;
               tol : in double_float; max_h : in natural32 );

  -- DESCRIPTION :
  --   Computes the multiplicity structure of the solutions of f(x) = 0.

  -- ON ENTRY :
  --   file    for writing diagnostics and results;
  --   output  flag to indicate if extra diagnostics are desired;
  --   f       a polynomial system;
  --   sols    isolated singular solutions of f(x) = 0;
  --   tol     tolerance used to determine the numerical rank;
  --   max_h   upper bound on the order of the differentials,
  --           used to form the nullity matrices.

  -- ON RETURN :  
  --   sols    solution with multiplicity field adjusted with the
  --           sum of the entrie in the Hilbert function. 

  procedure Driver_to_Multiplicity_Structure
             ( file : in file_type;
               f : in Poly_Sys; sols : in out Solution_List );

  -- DESCRIPTION :
  --   Prompts the user for the tolerance on the numerical rank and the
  --   maximal differential order, and asks whether extra output is needed.
  --   Then the multiplicity structure is computed for the solutions.

  -- ON ENTRY :
  --   file    for writing diagnostics and results;
  --   f       a polynomial system;
  --   sols    isolated singular solutions of f(x) = 0.

  -- ON RETURN :  
  --   sols    solution with multiplicity field adjusted with the
  --           sum of the entries in the Hilbert function. 

end Standard_Multiplicity_Structure;
