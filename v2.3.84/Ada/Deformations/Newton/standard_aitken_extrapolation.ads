with text_io;                            use text_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_VecVecs;           use Standard_Complex_VecVecs;
with Standard_Complex_VecLists;          use Standard_Complex_VecLists;

package Standard_Aitken_Extrapolation is

  function Conjugate_Product ( x,y : Vector ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the inner product of x with y,
  --   taking complex conjugates of the components of y.

  -- REQUIRED : x'range = y'range.

  function Norm2 ( x : Vector ) return Complex_Number;

  -- DESCRIPTION :
  --   The square of the 2-norm of x is the result of
  --   the conjugate product of x with itself.

  function Accelerate ( v0,v1,v2 : Vector ) return Vector;
  function Accelerate ( file : file_type; v0,v1,v2 : Vector ) return Vector;

  -- DESCRIPTION :
  --   Uses Aitken extrapolation to accelerate the sequence v0,v1,v2.
  --   Extra diagnostics are written to file when file provided as parameter.

  function Extrapolate ( v : VecVec ) return VecVec;
  function Extrapolate ( file : file_type; v : VecVec ) return VecVec;

  -- DESCRIPTION :
  --   Applies a version of Aitken extrapolation in several variables
  --   to the sequence of vectors in v.  
  --   The sequence on return has range v'first..v'last-2.

  -- REQUIRED : v'length > 2.

  function Extrapolate ( v : List ) return List;
  function Extrapolate ( file : file_type; v : List ) return List;

  -- DESCRIPTION :
  --   Applies the Aitken extrapolation to the list v,
  --   returning a list of length(v)-2.

  function Extrapolate ( v : List; k : positive ) return List;

  -- DESCRIPTION :
  --   Applies Aitken extrapolation k times.
  --   We have Extrapolate(v,1) = Extrapolate(v).

end Standard_Aitken_Extrapolation;
