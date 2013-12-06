with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Integer_Matrices;
with Standard_Complex_Matrices;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_SysFun;       use Standard_Complex_Laur_SysFun;
with Standard_Complex_Laur_Jacomats;     use Standard_Complex_Laur_Jacomats;
with Exponent_Vectors;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Multitasking_Polyhedral_Trackers is

-- DESCRIPTION :
--   This package offers routines to implement polyhedral continuation
--   methods for multithreading (shared memory parallel computing).
--   To control memory allocation and deallocation:
--   (1) tableau data structure for selection of start systems;
--   (2) inplace binomial system solvers for start solutions;
--   (3) tracking paths defined by polyhedral homotopies.

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
              return Complex_Number;

  -- DESCRIPTION :
  --   Returns the coefficient in cff matching with the exponent
  --   with value in the point pt.

  procedure Select_Coefficients
              ( q_c : in Standard_Complex_VecVecs.VecVec;
                q_e : in Standard_Integer_VecVecs.Array_of_VecVecs;
                pts : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                cff : out Standard_Complex_Vectors.Vector );

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

  function Create ( n : integer32 ) return Solution;

  -- DESCRIPTION :
  --   Creates an n-dimensional solution vector.

  function Create ( n,m : integer32 ) return Solution_List;

  -- DESCRIPTION :
  --   Returns a list of m solutions of dimension n.

  function Product_of_Diagonal
             ( A : Standard_Integer_Matrices.Matrix ) return integer32;

  -- DESCRIPTION :
  --   Returns the product of the elements on the diagonal of A.

  function Volume_of_Diagonal
             ( A : Standard_Integer_Matrices.Matrix ) return natural32;

  -- DESCRIPTION :
  --   Returns the product of the elements on the diagonal of A,
  --   multiplied by -1 if negative, so what is returned respresents
  --   the volume of the vectors spanned by A.

  function Volume_of_Cell
             ( A : Standard_Integer_Matrices.Matrix ) return natural32;

  -- DESCRIPTION :
  --   Returns the determinant of the matrix A,
  --   multiplied by -1 if negative, so what is returned respresents
  --   the volume of the vectors spanned by A.

  procedure Fully_Mixed_Start_Systems
              ( q : in Laur_Sys; mcc : in Mixed_Subdivision );

  -- DESCRIPTION :
  --   Solves the start systems of q defined by the mixed cells in mcc,
  --   using one single thread with inplace solvers.

  -- REQUIRED : the system is fully mixed and all cells are fine.

  procedure Semi_Mixed_Start_Systems
              ( q : in Laur_Sys; m : in natural32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Mixed_Subdivision );

  -- DESCRIPTION :
  --   Solves the start systems of q defined by the mixed cells in mcc,
  --   using one single thread with inplace solvers.

  -- REQUIRED : the system is semi-mixed and sorted just as the supports
  --   of the cells in mcc.  Moreover, all cells are fine.

  procedure Check_Solutions
              ( q : in Laur_Sys;
                mcc : in Mixed_Subdivision;
                sols : in Array_of_Solution_Lists;
                res : out Standard_Floating_Vectors.Vector );
  procedure Check_Solutions
              ( cff : in Standard_Complex_VecVecs.VecVec;
                exp : in Standard_Integer_VecVecs.Array_of_VecVecs;
                mcc : in Mixed_Subdivision;
                sols : in Array_of_Solution_Lists;
                res : out Standard_Floating_Vectors.Vector );

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
              ( q : in Laur_Sys;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Mixed_Subdivision;
                sols : in Array_of_Solution_Lists;
                res : out Standard_Floating_Vectors.Vector );

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

  procedure Reporting_Multithreaded_Solve_Start_Systems
              ( nt : in integer32; q : in Laur_Sys;
                r : in integer32; mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision;
                mixvol : out natural32;
                sols : out Array_of_Solution_Lists;
                res : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Solves the start systems of q defined by the mixed cells in mcc,
  --   using nt threads of execution.
  --   Intermediate output is written to screen.

  -- REQUIRED : the system is fully mixed and all cells are fine, and
  --   sols is an array of range 1..nt = res'range.

  -- ON ENTRY :
  --   nt       the number of tasks;
  --   q        a random coefficient start system;
  --   r        number of distinct supports;
  --   mix      type of mixture;
  --   mcc      a regular mixed-cell configuration.

  -- ON RETURN :
  --   mixvol   sum of the mixed volumes of the cells in mcc;
  --   sols     the start solutions computed by the tasks,
  --            in particular sols(i) is computed by task i,
  --            following a static workload assignment.
  --   res      res(i) contains the sums of all residuals of
  --            the start solutions computed by task i.

  procedure Silent_Multithreaded_Solve_Start_Systems
              ( nt : in integer32; q : in Laur_Sys;
                r : in integer32; mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision;
                mixvol : out natural32;
                sols : out Array_of_Solution_Lists;
                res : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Solves the start systems of q defined by the mixed cells in mcc,
  --   using nt threads of execution.  There is no intermediate output.

  -- REQUIRED : the system is fully mixed and all cells are fine, and
  --   sols is an array of range 1..nt = res'range.

  -- ON ENTRY :
  --   nt       the number of tasks;
  --   q        a random coefficient start system;
  --   r        number of distinct supports;
  --   mix      type of mixture;
  --   mcc      a regular mixed-cell configuration.

  -- ON RETURN :
  --   mixvol   sum of the mixed volumes of the cells in mcc;
  --   sols     the start solutions computed by the tasks,
  --            in particular sols(i) is computed by task i,
  --            following a static workload assignment.
  --   res      res(i) contains the sums of all residuals of
  --            the start solutions computed by task i.

-- (3) TRACKING PATHS DEFINED BY A POLYHEDRAL HOMOTOPY :

  procedure Silent_Multitasking_Path_Tracker
              ( q : in Laur_Sys; nt,n,r : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mcc : in Mixed_Subdivision; h : in Eval_Coeff_Laur_Sys;
                c : in Standard_Complex_VecVecs.VecVec;
                e : in Exponent_Vectors.Exponent_Vectors_Array;
                j : in Eval_Coeff_Jaco_Mat; mf : in Mult_Factors;
                sols : out Solution_List );
  procedure Reporting_Multitasking_Path_Tracker
              ( q : in Laur_Sys; nt,n,r : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mcc : in Mixed_Subdivision; h : in Eval_Coeff_Laur_Sys;
                c : in Standard_Complex_VecVecs.VecVec;
                e : in Exponent_Vectors.Exponent_Vectors_Array;
                j : in Eval_Coeff_Jaco_Mat; mf : in Mult_Factors;
                sols : out Solution_List );

  -- DESCRIPTION :
  --   Uses nt tasks to solve a random coefficient system q.
  --   The reporting version allows the user to monitor the progress 
  --   of the computations on screen.

  -- ON ENTRY :
  --   q        a random coefficient system;
  --   nt       number of tasks;
  --   n        ambient dimension;
  --   r        number of different supports;
  --   mix      type of mixture;
  --   lif      supports of q, with a lifting as occurred in mcc;
  --   mcc      a regular mixed-cell configuration;
  --   h        coefficient homotopy;
  --   c        coefficients of the homotopy h;
  --   j        coefficient Jacobian matrix;
  --   mf       multiplication factors in coefficient Jacobian matrix.
 
  -- ON RETURN :
  --   sols     all solutions of the system q.

  procedure Silent_Multitasking_Path_Tracker
              ( q : in Laur_Sys; nt,n,m : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Mixed_Subdivision; sols : out Solution_List );
  procedure Reporting_Multitasking_Path_Tracker
              ( file : in file_type; q : in Laur_Sys; nt,n,m : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Mixed_Subdivision; sols : out Solution_List );

  -- DESCRIPTION :
  --   Uses nt tasks to solve a random coefficient system q.
  --   The reporting version allows monitoring the computations.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   q        a random coefficient system;
  --   nt       number of tasks;
  --   n        ambient dimension, before the lifting;
  --   m        number of different supports;
  --   mix      type of mixture;
  --   mcc      a regular mixed-cell configuration.
 
  -- ON RETURN :
  --   sols     all solutions of the system q.

end Multitasking_Polyhedral_Trackers;
