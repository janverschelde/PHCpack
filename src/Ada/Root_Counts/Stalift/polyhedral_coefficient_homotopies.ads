with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_VecVecs;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Lists_of_Integer_Vectors;
with Lists_of_Floating_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Exponent_Vectors;

package Polyhedral_Coefficient_Homotopies is

-- DESCRIPTION :
--   In a polyhedral coefficient homotopy, the coefficients are
--   random complex numbers times a power of the continuation parameter.
--   This package offers tools to evaluate and work with polyhedral
--   coefficient homotopies.

  function Minimum ( v : Standard_Integer_VecVecs.VecVec ) return integer32;
  function Minimum ( v : Standard_Floating_VecVecs.VecVec )
                   return double_float;

  -- DESCRIPTION :
  --   Returns the smallest nonzero element in v, taken with absolute value.

  function Scale ( v : Standard_Integer_Vectors.Vector; s : integer32 )
                 return Standard_Floating_Vectors.Vector;
  function Scale ( v : Standard_Floating_Vectors.Vector; s : double_float )
                 return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector v divided by s.

  procedure Scale ( v : in out Standard_Floating_Vectors.Vector;
                    s : in double_float );

  -- DESCRIPTION :
  --   Divides every entry in v by s, for abs(s) > 1.0e-8.

  function Scale ( v : Standard_Integer_VecVecs.VecVec )
                 return Standard_Floating_VecVecs.VecVec;
  function Scale ( v : Standard_Floating_VecVecs.VecVec )
                 return Standard_Floating_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns an array of scaled vectors, after application of Scale
  --   to each vector in the vector of vectors v, with s = Minimum(v).
  --   After this operation, the smallest positive element in the
  --   vector on return should be one.

  procedure Scale ( v : in out Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Applies the Scale to each vector in v using the same smallest
  --   positive element in v as divisor.

  procedure Shift ( v : in out Standard_Integer_Vectors.Vector );
  procedure Shift ( v : in out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Assuming all elements in v are positive, the minimum in v is
  --   subtracted from all elements in v, so that on return, the smallest
  --   element in v equals zero.

  procedure Search_Lifting
              ( L : in Lists_of_Floating_Vectors.List;
                x : in Standard_Integer_Vectors.Vector;
                found : out boolean; y : out double_float );

  -- DESCRIPTION :
  --   Searches the lifting of the point x in the lifted list L.
  --   If found, then y equals the lifting value,
  --   otherwise, the value for y on return has no meaning.

  function Power_Transform
              ( e : Standard_Integer_VecVecs.VecVec;
                s : Lists_of_Integer_Vectors.List;
                normal : Standard_Integer_Vectors.Vector )
              return Standard_Integer_Vectors.Vector;
  function Power_Transform
              ( e : Standard_Integer_VecVecs.VecVec;
                s : Lists_of_Floating_Vectors.List;
                normal : Standard_Floating_Vectors.Vector )
              return Standard_Floating_Vectors.Vector;
  
  -- DESCRIPTION :
  --   Performs the coordinate transformation of the polyhedral homotopy,
  --   using an inner normal of a mixed cell.

  -- ON ENTRY :
  --   e        exponent vectors of a polynomial (before lifting);
  --   s        lifted supports of the polynomial;
  --   normal   inner normal to a mixed cell.

  -- ON RETURN :
  --   Vector with all inner products of the normal with the exponents
  --   in e so that points in the support of the mixed cell will all
  --   have a zero exponent (and thus can be used in a start system).

  procedure Power_Transform
              ( e : in Standard_Integer_VecVecs.VecVec;
                s : in Lists_of_Floating_Vectors.List;
                normal : in Standard_Floating_Vectors.Vector;
                L : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Procedural version of the coordinate transformation in the 
  --   polyhedral homotopy with an inner normal of a mixed cell.
  --   L contains the inner products of the normal with the exponents in e.

  -- REQUIRED : L'range = e'range.

  -- ON ENTRY :
  --   e        exponent vectors of a polynomial (before lifting);
  --   s        lifted supports of the polynomial;
  --   normal   inner normal to a mixed cell.

  -- ON RETURN :
  --   L        equals Power_Transform(e,s,normal).

  function Power_Transform
              ( e : Exponent_Vectors.Exponent_Vectors_Array;
                s : Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                mix,normal : Standard_Integer_Vectors.Vector )
              return Standard_Integer_VecVecs.VecVec;
  function Power_Transform
              ( e : Exponent_Vectors.Exponent_Vectors_Array;
                s : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mix : Standard_Integer_Vectors.Vector;
                normal : Standard_Floating_Vectors.Vector )
              return Standard_Floating_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Performs the coordinate transformation on all equations in
  --   the polyhedral homotopy.

  -- ON ENTRY :
  --   e        all exponent vectors of a polynomial system;
  --   s        lifted supports of a polynomial system;
  --   mix      type of mixture;
  --   normal   inner normal to a mixed cell.

  -- ON RETURN :
  --   For each exponent vector in e, there is a corresponding power
  --   in the vector of vectors on return.  Points in the supports of
  --   the mixed cell with the given normal have a power equal to zero.

  procedure Power_Transform
              ( e : in Exponent_Vectors.Exponent_Vectors_Array;
                s : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mix : in Standard_Integer_Vectors.Vector;
                normal : in Standard_Floating_Vectors.Vector;
                L : in out Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   This is the procedural version of the functino Power_Transform,
  --   on return L = Power_Transform(e,s,mix,normal).

  -- REQUIRED : L'range = e'range and sufficient memory is allocated
  --   for each L(i), i.e.: L(i)'range = e(i)'range.

-- EVALUATION ROUTINES :
--   for integer and floating point exponents m,
--   implemented as procedures or as functions below.
-- for standard double precision :
  procedure Eval ( c : in Standard_Complex_Vectors.Vector;
                   t : in double_float;
                   m : in Standard_Integer_Vectors.Vector;
                   ctm : in out Standard_Complex_Vectors.Vector );
  procedure Eval ( c : in Standard_Complex_Vectors.Vector;
                   t : in double_float;
                   m : in Standard_Floating_Vectors.Vector;
                   ctm : in out Standard_Complex_Vectors.Vector );
  procedure Eval ( c : in Standard_Complex_VecVecs.VecVec;
                   t : in double_float; 
                   m : in Standard_Integer_VecVecs.VecVec;
                   ctm : in out Standard_Complex_VecVecs.VecVec );
  procedure Eval ( c : in Standard_Complex_VecVecs.VecVec;
                   t : in double_float;
                   m : in Standard_Floating_VecVecs.VecVec;
                   ctm : in out Standard_Complex_VecVecs.VecVec );
-- for double double precision :
  procedure Eval ( c : in DoblDobl_Complex_Vectors.Vector;
                   t : in double_double;
                   m : in Standard_Integer_Vectors.Vector;
                   ctm : in out DoblDobl_Complex_Vectors.Vector );
  procedure Eval ( c : in DoblDobl_Complex_Vectors.Vector;
                   t : in double_double;
                   m : in Standard_Floating_Vectors.Vector;
                   ctm : in out DoblDobl_Complex_Vectors.Vector );
  procedure Eval ( c : in DoblDobl_Complex_VecVecs.VecVec;
                   t : in double_double; 
                   m : in Standard_Integer_VecVecs.VecVec;
                   ctm : in out DoblDobl_Complex_VecVecs.VecVec );
  procedure Eval ( c : in DoblDobl_Complex_VecVecs.VecVec;
                   t : in double_double;
                   m : in Standard_Floating_VecVecs.VecVec;
                   ctm : in out DoblDobl_Complex_VecVecs.VecVec );
-- for quad double precision :
  procedure Eval ( c : in QuadDobl_Complex_Vectors.Vector;
                   t : in quad_double;
                   m : in Standard_Integer_Vectors.Vector;
                   ctm : in out QuadDobl_Complex_Vectors.Vector );
  procedure Eval ( c : in QuadDobl_Complex_Vectors.Vector;
                   t : in quad_double;
                   m : in Standard_Floating_Vectors.Vector;
                   ctm : in out QuadDobl_Complex_Vectors.Vector );
  procedure Eval ( c : in QuadDobl_Complex_VecVecs.VecVec;
                   t : in quad_double; 
                   m : in Standard_Integer_VecVecs.VecVec;
                   ctm : in out QuadDobl_Complex_VecVecs.VecVec );
  procedure Eval ( c : in QuadDobl_Complex_VecVecs.VecVec;
                   t : in quad_double;
                   m : in Standard_Floating_VecVecs.VecVec;
                   ctm : in out QuadDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Returns ctm = c*t**m in ctm.  The implementation as a procedure
  --   seems unnatural, but functions returns vectors of vectors would
  --   lead to too many memory allocations and deallocations.

  -- WARNING :
  --   For t close to zero, the evaluation procedures crash.

end Polyhedral_Coefficient_Homotopies;
