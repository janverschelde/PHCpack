with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Double_Double_Vectors;
with Quad_Double_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Standard_Integer_Matrices;
with Standard_Integer64_Matrices;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_Systems;
with Exponent_Vectors;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Polyhedral_Start_Systems is

-- DESCRIPTION :
--   For the implementation of polyhedral continuation on shared memory
--   parallel computers with threads, we need to be more careful about
--   memory allocation and deallocation.  Therefore, this package offers
--   (1) tableau data structure for selection of start systems;
--   (2) inplace binomial system solvers for start solutions;
--   (3) allocating work space for exponents and coefficients.

-- (1) TABLEAU DATA STRUCTURES FOR START SYSTEM SELECTION :

  function Is_Equal
             ( x : Standard_Integer_Vectors.Link_to_Vector;
               y : Standard_Floating_Vectors.Link_to_Vector )
             return boolean;

  -- DESCRIPTION :
  --   Returns true for all i in x'range : x(i) = integer(y(i)).

  function Coefficient
              ( cff : Standard_Complex_Vectors.Link_to_Vector;
                exp : Standard_Integer_VecVecs.Link_to_VecVec;
                pt : Standard_Floating_Vectors.Link_to_Vector )
              return Standard_Complex_Numbers.Complex_Number;
  function Coefficient
              ( cff : DoblDobl_Complex_Vectors.Link_to_Vector;
                exp : Standard_Integer_VecVecs.Link_to_VecVec;
                pt : Standard_Floating_Vectors.Link_to_Vector )
              return DoblDobl_Complex_Numbers.Complex_Number;
  function Coefficient
              ( cff : QuadDobl_Complex_Vectors.Link_to_Vector;
                exp : Standard_Integer_VecVecs.Link_to_VecVec;
                pt : Standard_Floating_Vectors.Link_to_Vector )
              return QuadDobl_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Returns the coefficient in cff matching with the exponent
  --   with value in the point pt.

  procedure Select_Coefficients
              ( q_c : in Standard_Complex_VecVecs.VecVec;
                q_e : in Standard_Integer_VecVecs.Array_of_VecVecs;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                cff : out Standard_Complex_Vectors.Vector );
  procedure Select_Coefficients
              ( q_c : in Standard_Complex_VecVecs.VecVec;
                q_e : in Exponent_Vectors.Exponent_Vectors_Array;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                cff : out Standard_Complex_Vectors.Vector );
  procedure Select_Coefficients
              ( q_c : in DoblDobl_Complex_VecVecs.VecVec;
                q_e : in Standard_Integer_VecVecs.Array_of_VecVecs;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                cff : out DoblDobl_Complex_Vectors.Vector );
  procedure Select_Coefficients
              ( q_c : in DoblDobl_Complex_VecVecs.VecVec;
                q_e : in Exponent_Vectors.Exponent_Vectors_Array;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                cff : out DoblDobl_Complex_Vectors.Vector );
  procedure Select_Coefficients
              ( q_c : in QuadDobl_Complex_VecVecs.VecVec;
                q_e : in Standard_Integer_VecVecs.Array_of_VecVecs;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                cff : out QuadDobl_Complex_Vectors.Vector );
  procedure Select_Coefficients
              ( q_c : in QuadDobl_Complex_VecVecs.VecVec;
                q_e : in Exponent_Vectors.Exponent_Vectors_Array;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                cff : out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Selects the coefficients cff corresponding to the points in pts,
  --   given the tableau data structure of a system q.

  -- REQUIRED : cff'length = sum of Length_Of(pts(i)), for i in pts'range.

  -- ON ENTRY :
  --   q_c      coefficients of a polynomial system;
  --   q_e      corresponding exponents of the coefficients in q_c;
  --   pts      coordinates of a cell.

  -- ON RETURN :
  --   cff      coefficients corresponding to the points in pts.

  procedure Write_Tableau
              ( c : in Standard_Complex_VecVecs.VecVec;
                e : in Standard_Integer_VecVecs.Array_of_VecVecs );

  -- DESCRIPTION :
  --   Writes the tableau representation of a system,
  --   defined by coefficients in c and corresponding exponents in e.

  procedure Write_Tableau
              ( c : in Standard_Complex_Vectors.Vector;
                e : in Arrays_of_Floating_Vector_Lists.Array_of_Lists );

  -- DESCRIPTION :
  --   Writes the tableau representation for a system with coefficients in c
  --   and corresponding exponents in e.

  -- REQUIRED : c'length = sum of Length_Of(e(i)) for i in e'range.

  procedure Fully_Mixed_to_Binomial_Format
              ( c : in Standard_Complex_Vectors.Vector;
                e : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                A : out Standard_Integer_Matrices.Matrix;
                b : out Standard_Complex_Vectors.Vector );
  procedure Fully_Mixed_to_Binomial_Format
              ( c : in DoblDobl_Complex_Vectors.Vector;
                e : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                A : out Standard_Integer64_Matrices.Matrix;
                b : out DoblDobl_Complex_Vectors.Vector );
  procedure Fully_Mixed_to_Binomial_Format
              ( c : in QuadDobl_Complex_Vectors.Vector;
                e : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                A : out Standard_Integer64_Matrices.Matrix;
                b : out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Converts the coefficients and supports from the tableau format
  --   into the standard binomial system, as x^A = b.

  -- REQUIRED :
  --   e contains the supports of a fully mixed cell
  --   with corresponding coefficients in c, and
  --   if n is the ambient dimension (before lifting),
  --   then: A'range(1) = A'range(2) = 1..n = b'range = e'range.

  -- ON ENTRY :
  --   c        coefficients of a fully mixed cell subsystem,
  --            of range 1..2*n;
  --   e        supports of a fully mixed cell, e'range = 1..n.

  -- ON RETURN :
  --   A        matrix with exponent vectors;
  --   b        right hand side vector.

  procedure Select_Subsystem_to_Matrix_Format 
              ( q_c : in Standard_Complex_VecVecs.VecVec;
                q_e : in Standard_Integer_VecVecs.Array_of_VecVecs;
                mix : in Standard_Integer_Vectors.Vector;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                A : out Standard_Integer_Matrices.Matrix;
                C : out Standard_Complex_Matrices.Matrix;
                b : out Standard_Complex_Vectors.Vector );
  procedure Select_Subsystem_to_Matrix_Format 
              ( q_c : in Standard_Complex_VecVecs.VecVec;
                q_e : in Exponent_Vectors.Exponent_Vectors_Array;
                mix : in Standard_Integer_Vectors.Vector;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                A : out Standard_Integer_Matrices.Matrix;
                C : out Standard_Complex_Matrices.Matrix;
                b : out Standard_Complex_Vectors.Vector );
  procedure Select_Subsystem_to_Matrix_Format 
              ( q_c : in DoblDobl_Complex_VecVecs.VecVec;
                q_e : in Standard_Integer_VecVecs.Array_of_VecVecs;
                mix : in Standard_Integer_Vectors.Vector;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                A : out Standard_Integer64_Matrices.Matrix;
                C : out DoblDobl_Complex_Matrices.Matrix;
                b : out DoblDobl_Complex_Vectors.Vector );
  procedure Select_Subsystem_to_Matrix_Format 
              ( q_c : in DoblDobl_Complex_VecVecs.VecVec;
                q_e : in Exponent_Vectors.Exponent_Vectors_Array;
                mix : in Standard_Integer_Vectors.Vector;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                A : out Standard_Integer64_Matrices.Matrix;
                C : out DoblDobl_Complex_Matrices.Matrix;
                b : out DoblDobl_Complex_Vectors.Vector );
  procedure Select_Subsystem_to_Matrix_Format 
              ( q_c : in QuadDobl_Complex_VecVecs.VecVec;
                q_e : in Standard_Integer_VecVecs.Array_of_VecVecs;
                mix : in Standard_Integer_Vectors.Vector;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                A : out Standard_Integer64_Matrices.Matrix;
                C : out QuadDobl_Complex_Matrices.Matrix;
                b : out QuadDobl_Complex_Vectors.Vector );
  procedure Select_Subsystem_to_Matrix_Format 
              ( q_c : in QuadDobl_Complex_VecVecs.VecVec;
                q_e : in Exponent_Vectors.Exponent_Vectors_Array;
                mix : in Standard_Integer_Vectors.Vector;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                A : out Standard_Integer64_Matrices.Matrix;
                C : out QuadDobl_Complex_Matrices.Matrix;
                b : out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Selects the coefficients of q_c corresponding to the points in pts
  --   for a semi-mixed system.

  -- REQUIRED :
  --   pts contains the supports of a fine mixed cell, ordered
  --   along its type of mixture in mix and q is sorted as well;
  --   A'range(1) = A'range(2) = q_c'range = q_e'range
  --   = C'range(1) = C'range(2) = b'range.

  -- ON ENTRY :
  --   q_c      coefficients of a random coefficient system;
  --   q_e      exponents corresponding to the coefficients;
  --   mix      type of mixture;
  --   pts      points in a fine mixed cell.

  -- ON RETURN :
  --   A        the columns of A contain the points in pts,
  --            with the first point subtracted from it;
  --   C        coefficient corresponding to the points in pts;
  --   b        right hand side vector of the linear system C y = b,
  --            where b is minus the coefficient of the first point
  --            in each component of the support.

-- (2) INPLACE BINOMIAL SYSTEM SOLVERS FOR START SOLUTIONS :

  function Create ( n : integer32 ) return Standard_Complex_Solutions.Solution;
  function Create ( n : integer32 ) return DoblDobl_Complex_Solutions.Solution;
  function Create ( n : integer32 ) return QuadDobl_Complex_Solutions.Solution;

  -- DESCRIPTION :
  --   Creates an n-dimensional solution vector.

  function Create ( n,m : integer32 )
                  return Standard_Complex_Solutions.Solution_List;
  function Create ( n,m : integer32 )
                  return DoblDobl_Complex_Solutions.Solution_List;
  function Create ( n,m : integer32 )
                  return QuadDobl_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Returns a list of m solutions of dimension n.

  procedure Allocate
              ( first,last : in out Standard_Complex_Solutions.Solution_List;
                n,m : in integer32 );
  procedure Allocate
              ( first,last : in out DoblDobl_Complex_Solutions.Solution_List;
                n,m : in integer32 );
  procedure Allocate
              ( first,last : in out QuadDobl_Complex_Solutions.Solution_List;
                n,m : in integer32 );

  -- DESCRIPTION :
  --   Appends to the list with head first and tail last m solutions
  --   of dimension n.

  function Product_of_Diagonal
             ( A : Standard_Integer_Matrices.Matrix ) return integer32;
  function Product_of_Diagonal
             ( A : Standard_Integer64_Matrices.Matrix ) return integer64;

  -- DESCRIPTION :
  --   Returns the product of the elements on the diagonal of A.

  function Volume_of_Diagonal
             ( A : Standard_Integer_Matrices.Matrix ) return natural32;
  function Volume_of_Diagonal
             ( A : Standard_Integer64_Matrices.Matrix ) return natural64;

  -- DESCRIPTION :
  --   Returns the product of the elements on the diagonal of A,
  --   multiplied by -1 if negative, so what is returned respresents
  --   the volume of the vectors spanned by A.

  function Volume_of_Cell
             ( A : Standard_Integer_Matrices.Matrix ) return natural32;
  function Volume_of_Cell
             ( A : Standard_Integer64_Matrices.Matrix ) return natural64;

  -- DESCRIPTION :
  --   Returns the determinant of the matrix A,
  --   multiplied by -1 if negative, so what is returned respresents
  --   the volume of the vectors spanned by A.

  procedure Fully_Mixed_Start_Systems
              ( q : in Standard_Complex_Laur_Systems.Laur_Sys;
                mcc : in Mixed_Subdivision );
  procedure Fully_Mixed_Start_Systems
              ( q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mcc : in Mixed_Subdivision );
  procedure Fully_Mixed_Start_Systems
              ( q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mcc : in Mixed_Subdivision );

  -- DESCRIPTION :
  --   Solves the start systems of q defined by the mixed cells in mcc,
  --   using one single thread with inplace solvers.

  -- REQUIRED : the system is fully mixed and all cells are fine.

  procedure Semi_Mixed_Start_Systems
              ( q : in Standard_Complex_Laur_Systems.Laur_Sys;
                m : in natural32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Mixed_Subdivision );
  procedure Semi_Mixed_Start_Systems
              ( q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                m : in natural32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Mixed_Subdivision );
  procedure Semi_Mixed_Start_Systems
              ( q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                m : in natural32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Mixed_Subdivision );

  -- DESCRIPTION :
  --   Solves the start systems of q defined by the mixed cells in mcc,
  --   using one single thread with inplace solvers.

  -- REQUIRED : the system is semi-mixed and sorted just as the supports
  --   of the cells in mcc.  Moreover, all cells are fine.

  procedure Check_Solutions
              ( q : in Standard_Complex_Laur_Systems.Laur_Sys;
                mcc : in Mixed_Subdivision;
                sols : in Standard_Complex_Solutions.Array_of_Solution_Lists;
                res : out Standard_Floating_Vectors.Vector );
  procedure Check_Solutions
              ( cff : in Standard_Complex_VecVecs.VecVec;
                exp : in Standard_Integer_VecVecs.Array_of_VecVecs;
                mcc : in Mixed_Subdivision;
                sols : in Standard_Complex_Solutions.Array_of_Solution_Lists;
                res : out Standard_Floating_Vectors.Vector );
  procedure Check_Solutions
              ( q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mcc : in Mixed_Subdivision;
                sols : in DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                res : out Double_Double_Vectors.Vector );
  procedure Check_Solutions
              ( cff : in DoblDobl_Complex_VecVecs.VecVec;
                exp : in Standard_Integer_VecVecs.Array_of_VecVecs;
                mcc : in Mixed_Subdivision;
                sols : in DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                res : out Double_Double_Vectors.Vector );
  procedure Check_Solutions
              ( q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mcc : in Mixed_Subdivision;
                sols : in QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                res : out Quad_Double_Vectors.Vector );
  procedure Check_Solutions
              ( cff : in QuadDobl_Complex_VecVecs.VecVec;
                exp : in Standard_Integer_VecVecs.Array_of_VecVecs;
                mcc : in Mixed_Subdivision;
                sols : in QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                res : out Quad_Double_Vectors.Vector );

  -- DESCRIPTION :
  --   Computes the residuals for the solution computed in sols.

  -- ON ENTRY :
  --   q        a random coefficient system;
  --   cff      coefficients of a random coefficient start system q;
  --   exp      corresponding exponents of the tableau structure for q;
  --   mcc      a regular mixed-cell configuration;
  --   sols     sols(i) are solution computed by task i:
  --            cell k mod #tasks = i-1.

  -- ON RETURN :
  --   res      res(i) contains the sum of all residuals for sols(i).

  procedure Check_Solutions
              ( q : in Standard_Complex_Laur_Systems.Laur_Sys;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Mixed_Subdivision;
                sols : in Standard_Complex_Solutions.Array_of_Solution_Lists;
                res : out Standard_Floating_Vectors.Vector );
  procedure Check_Solutions
              ( q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Mixed_Subdivision;
                sols : in DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                res : out Double_Double_Vectors.Vector );
  procedure Check_Solutions
              ( q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Mixed_Subdivision;
                sols : in QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                res : out Quad_Double_Vectors.Vector );

  -- DESCRIPTION :
  --   For semi-mixed systems q, the residuals of sols are computed.

  -- ON ENTRY :
  --   q        a random coefficient system;
  --   mix      type of mixture;
  --   cff      coefficients of a random coefficient start system q;
  --   exp      corresponding exponents of the tableau structure for q;
  --   mcc      a regular mixed-cell configuration;
  --   sols     sols(i) are solution computed by task i:
  --            cell k mod #tasks = i-1.

  -- ON RETURN :
  --   res      res(i) contains the sum of all residuals for sols(i).

-- (3) ALLOCATING WORK SPACE FOR EXPONENTS AND COEFFICIENTS :

  procedure Allocate_Workspace_for_Exponents
              ( epv : in Exponent_Vectors.Exponent_Vectors_Array;
                dpw : in out Standard_Floating_VecVecs.Array_of_VecVecs );

  -- DESCRIPTION :
  --   Allocates space for the powers dpw in the polyhedral homotopy,
  --   using the dimensions of the exponent vectors array epv.
  --   The array of vecvecs gives every task its own work space.

  procedure Allocate_Workspace_for_Coefficients
              ( cff : in Standard_Complex_VecVecs.VecVec;
                cft : in out Standard_Complex_VecVecs.Array_of_VecVecs );
  procedure Allocate_Workspace_for_Coefficients
              ( cff : in DoblDobl_Complex_VecVecs.VecVec;
                cft : in out DoblDobl_Complex_VecVecs.Array_of_VecVecs );
  procedure Allocate_Workspace_for_Coefficients
              ( cff : in QuadDobl_Complex_VecVecs.VecVec;
                cft : in out QuadDobl_Complex_VecVecs.Array_of_VecVecs );

  -- DESCRIPTION :
  --   Allocates space for the coefficients cft in the polyhedral homotopy,
  --   using the dimensions of the coefficients in cff.
  --   The array of vecvecs gives every task its own work space.

end Polyhedral_Start_Systems;
