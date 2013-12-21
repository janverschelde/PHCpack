with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_VecVecs;           use Standard_Integer_VecVecs;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;
with Standard_Complex_Poly_Functions;
with Standard_Complex_Laur_Functions;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Standard_Complex_Laur_SysFun;       use Standard_Complex_Laur_SysFun;

package Integer_Lifting_Functions is

-- DESCRIPTION :
--   This package offers three different types of integer-valued
--   lifting functions: linear, polynomial, and point-wise.

-- LINEAR LIFTING :

  function Linear_Lift ( lf,v : Vector ) return Vector;
  function Linear_Lift ( lf : Vector; L : List ) return List;
  function Linear_Lift ( lf : VecVec;
                         L : Array_of_Lists ) return Array_of_Lists;

  -- DESCRIPTION :
  --   The last entry of the enlarged vector on return equals <lf,v>,
  --   for every vector v in the lists.

  -- REQUIRED : lf'range = v'range.

-- POLYNOMIAL LIFTING :

  function Polynomial_Lift
             ( lf : Standard_Complex_Polynomials.Poly; 
               v : Standard_Integer_Vectors.Vector )
             return Standard_Integer_Vectors.Vector;
  function Polynomial_Lift
             ( lf : Standard_Complex_Poly_Functions.Eval_Poly;
               v : Standard_Integer_Vectors.Vector )
             return Standard_Integer_Vectors.Vector;
  function Polynomial_Lift
             ( lf : Standard_Complex_Polynomials.Poly; L : List ) return List;
  function Polynomial_Lift
             ( lf : Standard_Complex_Laurentials.Poly; L : List ) return List;
  function Polynomial_Lift
             ( lf : Standard_Complex_Poly_Functions.Eval_Poly;
               L : List ) return List;
  function Polynomial_Lift
             ( lf : Standard_Complex_Laur_Functions.Eval_Poly;
               L : List ) return List;
  function Polynomial_Lift
             ( lf : Poly_Sys; L : Array_of_Lists ) return Array_of_Lists;
  function Polynomial_Lift
             ( lf : Laur_Sys; L : Array_of_Lists ) return Array_of_Lists;
  function Polynomial_Lift
             ( lf : Eval_Poly_Sys; L : Array_of_Lists ) return Array_of_Lists;
  function Polynomial_Lift
             ( lf : Eval_Laur_Sys; L : Array_of_Lists ) return Array_of_Lists;

  -- DESCRIPTION :
  --   As lf is complex valued, the rounded real part is used as lifting.

-- RANDOM LIFTING :

  function Random_Lift ( lflow,lfupp : integer32; v : Vector ) return Vector;
  function Random_Lift ( lflow,lfupp : integer32; L : List ) return List;
  function Random_Lift ( lflow,lfupp : Vector; L : Array_of_Lists )
                       return Array_of_Lists;

  -- DESCRIPTION :
  --   The lifting values are random integers between lflow and lfupp.

-- RANDOM LINEAR LIFTING :

  function Random_Linear_Lift
             ( lflow,lfupp : integer32; v : Vector ) return Vector;
  function Random_Linear_Lift
             ( lflow,lfupp : integer32; L : List ) return List;
  function Random_Linear_Lift
             ( lflow,lfupp : Vector; L : Array_of_Lists )
             return Array_of_Lists;

  -- DESCRIPTION :
  --   Linear lifting with random vectors of integers between lflow and lfupp.

-- POINT-WISE LIFTING :

  function Point_Lift ( lf : integer32; v : Vector ) return Vector;
  function Point_Lift ( lf : Vector; L : List ) return List;
  function Point_Lift ( lf : VecVec;
                        L : Array_of_Lists ) return Array_of_Lists;

  -- DESCRIPTION :
  --   The enlarged vector on return has as last component the lifting lf(k),
  --   where k is the position of the vector v in the lists.

  -- REQUIRED : Length_Of(l) = lf'length.

end Integer_Lifting_Functions;
