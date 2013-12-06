with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecMats;          use Standard_Complex_VecMats;

package Pieri_Homotopy is

-- DESCRIPTION :
--   To call the Pieri homotopies from a function in C,
--   the Pieri homotopies are organized as a state machine.

  procedure Initialize_Dimensions ( m,p,q : in natural32 );

  -- DESCRIPTION :
  --   The problem is to compute all curves of degree q,
  --   which produce p-planes which meet mp + q(m+p) given
  --   m-planes at specified interpolation points.

  procedure Initialize_Input_Planes ( planes : in VecMat );

  -- DESCRIPTION :
  --   Stores the input m-planes in the state machine.
  --   The range of planes should be 1..mp+q(m+p).

  procedure Initialize_Interpolation_Points
              ( points : in Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Stores the interpolation points in the state machine.
  --   The range of points should be 1..mp+q(m+p).

  procedure Initialize ( m,p,q : in natural32; planes : in VecMat;
                         points : in Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Executes the above three initialization routines in sequence.

  procedure Store_Start_Pivots
              ( top,bottom : in Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   This procedure defines the top and bottom pivots of the 
  --   localization pattern of the solution p-plane producing curve
  --   of degree q at the start of the homotopy.

  procedure Store_Target_Pivots
              ( top,bottom : in Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   This procedure defines the top and bottom pivots of the 
  --   localization pattern of the solution p-plane producing curve
  --   of degree q at the end of the homotopy.

  function Degree_of_Freedom
             ( top,bottom : Standard_Natural_Vectors.Vector )
             return natural32;

  -- DESCRIPTION :
  --   Returns the number of free variables in a localization pattern
  --   defined by top and bottom pivots.  This is a useful routine
  --   from prompting the user to supply a start solution vector.

  procedure Store_Start ( x : in Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   This vector provides actual numbers for the solution curve
  --   at the start of the homotopy.
  --   The length of the vector must equal the number of free
  --   positions defined by the localization pattern at the start.

  procedure Retrieve_Target ( x : out Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   This vector provides actual numbers for the solution curve
  --   at the end of the homotopy.
  --   The length of the vector must equal the number of free
  --   positions defined by the localization pattern at the end.

  procedure Track_Path;
  procedure Track_Path ( file : file_type );

  -- DESCRIPTION :
  --   Runs one path starting at the stored start solution
  --   and provides one target solution, ready for retrieval.
  --   When a file is provided, diagnostics are written.

  function Verify_Determinants return double_float;
  function Verify_Determinants ( file : in file_type ) return double_float;

  -- DESCRIPTION :
  --   Returns the sum of the absolute values of the determinants at 
  --   the intersection conditions at the current target solution.
  --   When a file is provided, diagnostics are written.

  procedure Clear;

  -- DESCRIPTION :
  --   Destroys the state machine.

end Pieri_Homotopy;
