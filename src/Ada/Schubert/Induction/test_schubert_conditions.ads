with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Natural_Matrices;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Poly_Matrices;
with DoblDobl_Complex_Poly_Matrices;
with QuadDobl_Complex_Poly_Matrices;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Brackets;                          use Brackets;
with Bracket_Monomials;                 use Bracket_Monomials;
with Standard_Bracket_Polynomials;      use Standard_Bracket_Polynomials;
with Checker_Posets;                    use Checker_Posets;

package Test_Schubert_Conditions is

-- DESCRIPTION :
--   Test on the symbolic and numeric Schubert intersection conditions.

  procedure Symbolic_Plane ( n,k : in integer32 );

  -- DESCRIPTION :
  --   Shows a symbolic representation for a general k-plane in n-space.

  function Identity ( n,nv : integer32 )
             return Standard_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the n-by-n identity matrix where the ones are polynomials
  --   in nv variables.

  procedure Elaborate_Flag_Minors
               ( n,k,f,i : in natural32; fm : in Bracket_Polynomial;
                 A : in Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Elaborates the flag minors as defined by the bracket monomial,
  --   for k-planes in n-space meeting a dimension f-plane
  --   at intersection dimension equal to i.

  procedure Impose_Schubert_Condition ( n,k : in natural32 );

  -- DESCRIPTION :
  --   Generates the polynomial equations imposed by one Schubert
  --   condition on a k-plane in n-space.

  procedure Standard_Test_Minor_Expansions ( n,k : in integer32 );

  -- DESCRIPTION :
  --   Shows all k-by-k minors of a general k-plane in n-space,
  --   represented by a general symbolic localization pattern,
  --   as a matrix of polynomials with standard double coefficients.

  procedure DoblDobl_Test_Minor_Expansions ( n,k : in integer32 );

  -- DESCRIPTION :
  --   Shows all k-by-k minors of a general k-plane in n-space,
  --   represented by a general symbolic localization pattern,
  --   as a matrix of polynomials with double double coefficients.

  procedure QuadDobl_Test_Minor_Expansions ( n,k : in integer32 );

  -- DESCRIPTION :
  --   Shows all k-by-k minors of a general k-plane in n-space,
  --   represented by a general symbolic localization pattern,
  --   as a matrix of polynomials with quad double coefficients.

  procedure Test_Minor_Expansions ( n,k : integer32 );

  -- DESCRIPTION :
  --   Prompts the user for the precision and then calls the
  --   appropriate testing procedure.

  function Read_Localization_Map ( n,k : integer32 )
             return Standard_Natural_Matrices.Matrix;

  -- DESCRIPTION :
  --   Prompts the user to choose for the default general localization map,
  --   or to enter an n-by-k matrix of natural numbers to represent
  --   a localization map for a k-plane in n-space.

  procedure Standard_Flag_Conditions ( n,k : in integer32; b : in Bracket );

  -- DESCRIPTION :
  --   Enumerates all conditions imposed by the flag,
  --   according to the k-bracket, in standard double precision.

  procedure DoblDobl_Flag_Conditions ( n,k : in integer32; b : in Bracket );

  -- DESCRIPTION :
  --   Enumerates all conditions imposed by the flag,
  --   according to the k-bracket, in double double precision.

  procedure Flag_Conditions ( n,k : in integer32; bm : in Bracket_Monomial );

  -- DESCRIPTION :
  --   Enumerates all conditions imposed by as many flags as factors in bm,
  --   according to the k-brackets in the bracket monomial.

  procedure Test_Schubert_Conditions ( n,k : in integer32 );

  -- DESCRIPTION :
  --   Prompts for a bracket for a k-plane in n-space,
  --   and then generates the conditions imposed by the bracket.

  function Generate_Point
             ( n,k : integer32; b : Bracket;
               flag : Standard_Complex_Matrices.Matrix ) 
             return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Generates a random point in the variety of k-planes in n-space
  --   that satisfies the Schubert condition imposed by the bracket b
  --   and the flag, in standard double precision.

  function Generate_Point
             ( n,k : integer32; b : Bracket;
               flag : DoblDobl_Complex_Matrices.Matrix ) 
             return DoblDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Generates a random point in the variety of k-planes in n-space
  --   that satisfies the Schubert condition imposed by the bracket b
  --   and the flag, in double double precision.

  function Generate_Point
             ( n,k : integer32; b : Bracket;
               flag : QuadDobl_Complex_Matrices.Matrix ) 
             return QuadDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Generates a random point in the variety of k-planes in n-space
  --   that satisfies the Schubert condition imposed by the bracket b
  --   and the flag, in double double precision.

  procedure Eliminate ( x : in out Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Eliminates elements on top and bottom of x to fit it into
  --   a general localization pattern.

  -- REQUIRED : x is a generic k-plane and no pivoting is needed.

  procedure Divide_Pivots ( x : in out Standard_Complex_Matrices.Matrix;
                            b : in Bracket );

  -- DESCRIPTION :
  --   Divides the columns of x by x(b(i),i), for i in b'range.

  procedure Eliminate_Pivots
              ( x : in out Standard_Complex_Matrices.Matrix;
                b : in Bracket );

  -- DESCRIPTION :
  --   Makes the numbers at the right of the pivots equal to zero
  --   by making column combinations so x fits the pattern imposed by b.
 
  function Generate_Standard_Point
             ( n,k : integer32; b : Bracket;
               flag : Standard_Complex_Matrices.Matrix ) 
             return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   After generating a random k-plane in n-space that satisfies the
  --   intersection conditions imposed by the bracket b and the flag,
  --   eliminates to fit into a standard localization pattern.

  function Full_Localization_Map
             ( n,k : in integer32 ) return Standard_Natural_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a localization map that allows all entries to be variable.

  function Stiefel_Localization_Map
             ( n,k : in integer32; b : in Bracket )
             return Standard_Natural_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a localization map for a k-plane in n-space that
  --   represents the Stiefel coordinates for the Schubert variety
  --   defined by the bracket b.

  function Full_Flatten
             ( x : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Stacks the columns of x one after the other into a big vector.
  --   Assumed is that x fits a full localization map.

  function Full_Flatten
             ( x : DoblDobl_Complex_Matrices.Matrix )
             return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Stacks the columns of x one after the other into a big vector.
  --   Assumed is that x fits a full localization map.

  function Full_Flatten
             ( x : QuadDobl_Complex_Matrices.Matrix )
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Stacks the columns of x one after the other into a big vector.
  --   Assumed is that x fits a full localization map.

  function General_Flatten
             ( x : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Assuming that x fits a general localization pattern,
  --   the columns of x are stacked one after the other into
  --   one long vector.

  procedure Standard_Point_Test_at_Conditions ( n,k : in integer32 );

  -- DESCRIPTION :
  --   Prompts the user for a k-bracket and generates the corresponding
  --   polynomial system that expresses the Schubert conditions,
  --   in standard double precision.

  procedure DoblDobl_Point_Test_at_Conditions ( n,k : in integer32 );

  -- DESCRIPTION :
  --   Prompts the user for a k-bracket and generates the corresponding
  --   polynomial system that expresses the Schubert conditions,
  --   in double double precision.

  procedure QuadDobl_Point_Test_at_Conditions ( n,k : in integer32 );

  -- DESCRIPTION :
  --   Prompts the user for a k-bracket and generates the corresponding
  --   polynomial system that expresses the Schubert conditions,
  --   in quad double precision.

  procedure Point_Test_at_Conditions ( n,k : in integer32 );

  -- DESCRIPTION :
  --   Prompts the user for the working precision,
  --   and then does the point test at the conditions.

  procedure Truncate_Triangular_Part
              ( A : in out Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Sets all elements below the diagonal equal to zero.
  --   The diagonal elements are set to one.

  procedure Point_Test_at_Minimal_Conditions ( n,k : in integer32 );

  -- DESCRIPTION :
  --   Prompts the user for a k-bracket and generates the corresponding
  --   polynomial system that expresses the Schubert conditions,
  --   using the more efficient formulation for the Schubert problem.

  procedure Cheater_Homotopy ( n,k,nq : in integer32; b : in Bracket );

  -- DESCRIPTION :
  --   Generates a random flag and generates a k-plane that satisfies
  --   all the Schubert conditions imposed by the k-bracket.

  procedure Test_Cheater_Homotopy ( n,k : in integer32 );

  -- DESCRIPTION :
  --   Sets up a cheater's homotopy for one point.

  procedure Generalizing_Moving_Flags ( n,dim : in integer32 );

  -- DESCRIPTION :
  --   Shows all boards, flags and transformations in a generalizing
  --   moving flag homotopy, using as many indeterminates as dim.

  -- REQUIRED : the symbol table is initialized properly.

  procedure Symbolic_Moving_Flags ( n : in integer32 );

  -- DESCRIPTION :
  --   Allows the user to view all moving flags in n-space,
  --   as defined by descending and rising black checkers.
  --   In particular, the descending black checker determines the
  --   column of the new variable and its row is determined by the
  --   rising black checker.

  procedure Symbolic_Localization_Patterns ( n,k : in integer32 );

  -- DESCRIPTION :
  --   Allows the user to see all localization patterns
  --   along a path in an intersection poset.

  procedure Define_Moving_Flag_Homotopy ( n,k : in integer32 );

  -- DESCRIPTION :
  --   Defines a moving flag homotopy for a k-plane in n-space.


  procedure Test_One_Flag_Homotopy ( n,k : in integer32 );

  -- DESCRIPTION:
  --   Prompts for a permutation and a configuration of white checkers
  --   to test a homotopy to move one flag for a k-plane in n-space.

  procedure Write_Coefficients ( ps : in Poset );

  -- DESCRIPTION :
  --   Writes the Littlewood-Richardson coefficients of root
  --   and at the leaves of the poset ps.

  procedure Create_Schubert_Poset
               ( n,k,nb : in integer32; b : Bracket_Monomial );

  -- DESCRIPTION :
  --   A Schubert poset defines the generalizations of the moving flag,
  --   in resolving an intersection condition, starting at the leaves.

  procedure Define_Schubert_Systems ( n,k : in integer32 );

  -- DESCRIPTION :
  --   Defines all polynomial systems occurring in the resolution
  --   of a general Schubert problem.

  procedure Verify_a_Solution
              ( n,k : in integer32; cond : in Bracket;
                flag : in Standard_Complex_Matrices.Matrix;
                locmap : in Standard_Natural_Matrices.Matrix;
                sol : in Solution; tol : in double_float );

  -- DESCRIPTION :
  --   Verifies the solution with respect to the flag and the condition.

  procedure Verify_Solutions
              ( n,k : in integer32; cond : in Bracket;
                flag : in Standard_Complex_Matrices.Matrix;
                sols : in Solution_List; tol : in double_float );

  -- DESCRIPTION :
  --   Verifies the solutions with respect to the flag and the condition.

  procedure Evaluate_Solutions
              ( n,k : in integer32; cond : in Bracket;
                flag : in Standard_Complex_Matrices.Matrix;
                sols : in Solution_List );

  -- DESCRIPTION :
  --   Evaluates the solutions at the expanded polynomial equations.

  procedure Verify_Solutions_of_Schubert_Problem ( n,k : in integer32 );

  -- DESCRIPTION :
  --   Prompts the user for the file name for the coordinates of a flag,
  --   represented as an complex n-by-n matrix, then asks for a list of
  --   solutions, and then finally for a bracket condition.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu a prompts for a test.

end Test_Schubert_Conditions;
