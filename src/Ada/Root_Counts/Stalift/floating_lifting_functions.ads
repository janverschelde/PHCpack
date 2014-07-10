with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Vectors;          use Standard_Floating_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;
with Standard_Complex_Poly_Functions;
with Standard_Complex_Laur_Functions;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_SysFun;       use Standard_Complex_Laur_SysFun;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Lists_of_Floating_Vectors;          use Lists_of_Floating_Vectors;
with Arrays_of_Floating_Vector_Lists;    use Arrays_of_Floating_Vector_Lists;
with Floating_Faces_of_Polytope;         use Floating_Faces_of_Polytope;

package Floating_Lifting_Functions is

-- DESCRIPTION :
--   This package provides a suite of floating-point lifting functions.

-- RANDOM FLOATING-POINT LIFTING :

  function Random_Lift ( lflow,lfupp : double_float ) return double_float;
  function Random_Lift ( v : Vector; lflow,lfupp : double_float ) return Vector;
  function Random_Lift ( L : List; lflow,lfupp : double_float ) return List;
  function Random_Lift ( L : Array_of_Lists; lflow,lfupp : Vector )
                       return Array_of_Lists;
  -- DESCRIPTION :
  --   Random lifting values between lflow and lfupp.

-- LINEAR LIFTING FUNCTIONS :

  function Linear_Lift ( x,v : Vector ) return Vector;
  function Linear_Lift ( f : Face; v : Vector ) return Face;
  function Linear_Lift ( L : List; v : Vector ) return List;
  function Linear_Lift ( f : Faces; v : Vector ) return Faces;

  -- DESCRIPTION :
  --   Returns a linearly lifted vector, list or faces.

-- RANDOM FLOATING-POINT LINEAR LIFTING FUNCTIONS :

  function Random ( n : integer32; lflow,lfupp : double_float ) return Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 1..n with randomly generated numbers,
  --   in [lflow,lfupp].  Random linear lifting functions are provided
  --   by using this randomly generated vector.

-- POLYNOMIAL LIFTING FUNCTIONS :

  function Polynomial_Lift
             ( lf : Standard_Complex_Polynomials.Poly;
               x : Standard_Floating_Vectors.Vector )
             return Standard_Floating_Vectors.Vector;
  function Polynomial_Lift
             ( lf : Standard_Complex_Poly_Functions.Eval_Poly;
               x : Standard_Floating_Vectors.Vector )
             return Standard_Floating_Vectors.Vector;
  function Polynomial_Lift ( lf : Standard_Complex_Polynomials.Poly;
                             L : List ) return List;
  function Polynomial_Lift ( lf : Standard_Complex_Poly_Functions.Eval_Poly;
                             L : List ) return List;
  function Polynomial_Lift ( lf : Standard_Complex_Poly_Systems.Poly_Sys;
                             L : Array_of_Lists )
                           return Array_of_Lists;
  function Polynomial_Lift ( lf : Standard_Complex_Laur_Systems.Laur_Sys;
                             L : Array_of_Lists )
                           return Array_of_Lists;
  function Polynomial_Lift ( lf : Eval_Poly_Sys; L : Array_of_Lists )
                           return Array_of_Lists;
  function Polynomial_Lift ( lf : Eval_Laur_Sys; L : Array_of_Lists )
                           return Array_of_Lists;

-- FOR STABLE MIXED VOLUMES :

  function Max_Degree ( p : in Standard_Complex_Poly_Systems.Poly_Sys )
                      return integer32;
  function Max_Degree ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys )
                      return integer32;
  function Max_Degree ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys )
                      return integer32;
  function Max_Degree ( p : in Standard_Complex_Laur_Systems.Laur_Sys )
                      return integer32;
  function Max_Degree ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys )
                      return integer32;
  function Max_Degree ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys )
                      return integer32;

  -- DESCRIPTION :
  --   Returns the maximal degree of all polynomials in the system p.

  function Lifting_Bound ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                           max : double_float := 1.0E+8 )
                         return double_float;
  function Lifting_Bound ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                           max : double_float := 1.0E+8 )
                         return double_float;
  function Lifting_Bound ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                           max : double_float := 1.0E+8 )
                         return double_float;
  function Lifting_Bound ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                           max : double_float := 1.0E+8 )
                         return double_float;
  function Lifting_Bound ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                           max : double_float := 1.0E+8 )
                         return double_float;
  function Lifting_Bound ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                           max : double_float := 1.0E+8 )
                         return double_float;

  -- DESCRIPTION :
  --   Returns a conservative lower bound for the lifting: n*(n+1)*d^n, 
  --   where d is the maximal degree of the polynomials in the system p.
  --   For larger polynomial systems, this bound could be larger than
  --   the largest integer.

  function Stable_Lift ( L : List; b : double_float ) return List;
  function Stable_Lift ( L : Array_of_Lists; b : double_float )
                       return Array_of_Lists;

  -- DESCRIPTION :
  --   Applies a random lifting between 0.0 and 1.0 for the points in l.
  --   In case the origin does not belong to l, then it will be added
  --   with lifting value equal to b.

  function Lifting_Bound
              ( lifted : Arrays_of_Floating_Vector_Lists.Array_of_Lists ) 
              return double_float;

  -- DESCRIPTION :
  --   Returns the maximal lifting value of the origins to detect whether
  --   a user-given subdivision was used for stable mixed volumes.

end Floating_Lifting_Functions;
