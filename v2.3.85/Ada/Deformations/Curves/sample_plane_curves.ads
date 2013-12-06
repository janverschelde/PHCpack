with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Polynomials;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Polynomials;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Polynomials;
with Multprec_Complex_Numbers;
with Multprec_Complex_Vectors;
with Multprec_Complex_Polynomials;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Multprec_Complex_Solutions;

package Sample_Plane_Curves is

-- DESCRIPTION :
--   Sampling points on a a plane algebraic curve reduces to the
--   solutions of a polynomial in one variable when the line that
--   intersects the polynomial defining the algebraic curve is vertical,
--   i.e.: if the first coordinates takes a fixed value.
--   The routines in this package are organized in three parts:
--   (1) tools to store the points in lists of solutions;
--   (2) wrappers to the univariate polynomial root finders;
--   (3) making a grid of points for testing purposes.

-- PART I : points on curves as solutions with same first coordinate

  function Insert_First_Coordinate
              ( s : Standard_Complex_Solutions.Solution;
                x : Standard_Complex_Numbers.Complex_Number )
              return Standard_Complex_Solutions.Solution;
  function Insert_First_Coordinate
              ( s : DoblDobl_Complex_Solutions.Solution;
                x : DoblDobl_Complex_Numbers.Complex_Number )
              return DoblDobl_Complex_Solutions.Solution;
  function Insert_First_Coordinate
              ( s : QuadDobl_Complex_Solutions.Solution;
                x : QuadDobl_Complex_Numbers.Complex_Number )
              return QuadDobl_Complex_Solutions.Solution;
  function Insert_First_Coordinate
              ( s : Multprec_Complex_Solutions.Solution;
                x : Multprec_Complex_Numbers.Complex_Number )
              return Multprec_Complex_Solutions.Solution;

  -- DESCRIPTION :
  --   Returns a solution with x inserted to s at its first position.

  -- REQUIRED : s.n = 1.

  function Insert_First_Coordinate
              ( s : Standard_Complex_Solutions.Solution_List;
                x : Standard_Complex_Numbers.Complex_Number )
              return Standard_Complex_Solutions.Solution_List;
  function Insert_First_Coordinate
              ( s : DoblDobl_Complex_Solutions.Solution_List;
                x : DoblDobl_Complex_Numbers.Complex_Number )
              return DoblDobl_Complex_Solutions.Solution_List;
  function Insert_First_Coordinate
              ( s : QuadDobl_Complex_Solutions.Solution_List;
                x : QuadDobl_Complex_Numbers.Complex_Number )
              return QuadDobl_Complex_Solutions.Solution_List;
  function Insert_First_Coordinate
              ( s : Multprec_Complex_Solutions.Solution_List;
                x : Multprec_Complex_Numbers.Complex_Number )
              return Multprec_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Inserts x as first coordinate to every solution in the list s.

  function Second_Coordinate
              ( s : Standard_Complex_Solutions.Solution_List )
              return Standard_Complex_Vectors.Vector;
  function Second_Coordinate
              ( s : DoblDobl_Complex_Solutions.Solution_List )
              return DoblDobl_Complex_Vectors.Vector;
  function Second_Coordinate
              ( s : QuadDobl_Complex_Solutions.Solution_List )
              return QuadDobl_Complex_Vectors.Vector;
  function Second_Coordinate
              ( s : Multprec_Complex_Solutions.Solution_List )
              return Multprec_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Extracts the second coordinate of every solution in s.
  --   The vector on return has range 1..Length_Of(s).

  procedure Update_Samples
              ( s : in out Standard_Complex_Solutions.Solution_List;
                x : in Standard_Complex_Numbers.Complex_Number;
                y : in Standard_Complex_Vectors.Vector );
  procedure Update_Samples
              ( s : in out DoblDobl_Complex_Solutions.Solution_List;
                x : in DoblDobl_Complex_Numbers.Complex_Number;
                y : in DoblDobl_Complex_Vectors.Vector );
  procedure Update_Samples
              ( s : in out QuadDobl_Complex_Solutions.Solution_List;
                x : in QuadDobl_Complex_Numbers.Complex_Number;
                y : in QuadDobl_Complex_Vectors.Vector );
  procedure Update_Samples
              ( s : in out Multprec_Complex_Solutions.Solution_List;
                x : in Multprec_Complex_Numbers.Complex_Number;
                y : in Multprec_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Updates the list of solutions in s with x as first coordinate
  --   for all solutions and the corresponding coordinates in y as
  --   the second coordinate.

  -- REQUIRED : Length_Of(s) = y'length.

-- PART II : calling univariate polynomial solvers

  function Sample ( p : Standard_Complex_Polynomials.Poly;
                    x : Standard_Complex_Numbers.Complex_Number )
                  return Standard_Complex_Solutions.Solution_List;
  function Sample ( p : DoblDobl_Complex_Polynomials.Poly;
                    x : DoblDobl_Complex_Numbers.Complex_Number )
                  return DoblDobl_Complex_Solutions.Solution_List;
  function Sample ( p : QuadDobl_Complex_Polynomials.Poly;
                    x : QuadDobl_Complex_Numbers.Complex_Number )
                  return QuadDobl_Complex_Solutions.Solution_List;
  function Sample ( p : Multprec_Complex_Polynomials.Poly;
                    x : Multprec_Complex_Numbers.Complex_Number;
                    size : natural32 )
                  return Multprec_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Replaces the first variable of p by the value in x
  --   and computes the roots of the resulting univariable polynomial.
  --   Every solution in the list on return has two coordinates,
  --   the first coordinate equal to the value of x.
  --   In the multiprecision version, the size of the numbers is
  --   given by the value of size.

  -- REQUIRED : Number_of_Unknowns(p) = 1.  

  procedure Sample ( p : in Standard_Complex_Polynomials.Poly;
                     x : in Standard_Complex_Numbers.Complex_Number; 
                     sols : in out Standard_Complex_Solutions.Solution_List;
                     fail : out boolean );
  procedure Sample ( p : in DoblDobl_Complex_Polynomials.Poly;
                     x : in DoblDobl_Complex_Numbers.Complex_Number; 
                     sols : in out DoblDobl_Complex_Solutions.Solution_List;
                     fail : out boolean );
  procedure Sample ( p : in QuadDobl_Complex_Polynomials.Poly;
                     x : in QuadDobl_Complex_Numbers.Complex_Number; 
                     sols : in out QuadDobl_Complex_Solutions.Solution_List;
                     fail : out boolean );
  procedure Sample ( p : in Multprec_Complex_Polynomials.Poly;
                     x : in Multprec_Complex_Numbers.Complex_Number; 
                     size : in natural32;
                     sols : in out Multprec_Complex_Solutions.Solution_List;
                     fail : out boolean );

  -- DESCRIPTION :
  --   The solutions in sols are used as starting point in the
  --   iterative method to compute the univariate roots.
  --   The working precision in the multiprecision version is determined
  --   by the value of size.

  procedure Standard_Evaluate
              ( p : in Standard_Complex_Polynomials.Poly;
                s : in Standard_Complex_Solutions.Solution_List );
  procedure DoblDobl_Evaluate
              ( p : in DoblDobl_Complex_Polynomials.Poly;
                s : in DoblDobl_Complex_Solutions.Solution_List );
  procedure QuadDobl_Evaluate
              ( p : in QuadDobl_Complex_Polynomials.Poly;
                s : in QuadDobl_Complex_Solutions.Solution_List );
  procedure Multprec_Evaluate
              ( p : in Multprec_Complex_Polynomials.Poly;
                s : in Multprec_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Writes the values of the solutions to screen.

-- PART III : building a grid of points

  function Perturb ( x : Standard_Complex_Numbers.Complex_Number;
                     m : double_float )
                   return Standard_Complex_Numbers.Complex_Number;
  function Perturb ( x : DoblDobl_Complex_Numbers.Complex_Number;
                     m : double_float )
                   return DoblDobl_Complex_Numbers.Complex_Number;
  function Perturb ( x : QuadDobl_Complex_Numbers.Complex_Number;
                     m : double_float )
                   return QuadDobl_Complex_Numbers.Complex_Number;
  function Perturb ( x : Multprec_Complex_Numbers.Complex_Number;
                     m : double_float )
                   return Multprec_Complex_Numbers.Complex_Number;

  -- DESCRIPTION :
  --   Returns x + m*z, where z is a random number with |z| = 1.

  function Branch ( s : Standard_Complex_Solutions.Array_of_Solution_Lists;
                    k : natural32 )
                  return Standard_Complex_Solutions.Solution_List;
  function Branch ( s : DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                    k : natural32 )
                  return DoblDobl_Complex_Solutions.Solution_List;
  function Branch ( s : QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                    k : natural32 )
                  return QuadDobl_Complex_Solutions.Solution_List;
  function Branch ( s : Multprec_Complex_Solutions.Array_of_Solution_Lists;
                    k : natural32 )
                  return Multprec_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   The list on return contains the k-th solution of every list in s.

  -- REQUIRED : k is in the range 1..Length_Of(s(i)), for i in s'range.

  procedure Standard_Slice
              ( p : in Standard_Complex_Polynomials.Poly;
                d : in natural32;
                s : out Standard_Complex_Solutions.Array_of_Solution_Lists );
  procedure DoblDobl_Slice
              ( p : in DoblDobl_Complex_Polynomials.Poly;
                d : in natural32;
                s : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists );
  procedure QuadDobl_Slice
              ( p : in QuadDobl_Complex_Polynomials.Poly;
                d : in natural32;
                s : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists );
  procedure Multprec_Slice
              ( p : in Multprec_Complex_Polynomials.Poly;
                d,size : in natural32;
                s : out Multprec_Complex_Solutions.Array_of_Solution_Lists );

  -- DESCRIPTION :
  --   Prepares d + 1 slices of the plane algebraic curve defined
  --   by the polynomial in p, where d = Degree(p).
  --   The working precision in the multiprecision version is in size.

  -- REQUIRED : s'range = 0..d.

end Sample_Plane_Curves;
