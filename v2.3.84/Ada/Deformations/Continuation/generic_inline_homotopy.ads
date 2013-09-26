with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;

generic
 
  with procedure Homotopy_Constants ( a : in Complex_Number; k : in positive );
  with function Eval_Homotopy ( x : Vector; t : Complex_Number ) return Vector;
  with function Diff_Homotopy ( x : Vector; t : Complex_Number ) return Matrix;
  with function Diff_Homotopy ( x : Vector; t : Complex_Number ) return Vector;

package Generic_Inline_Homotopy is

-- DESCRIPTION :
--   Generic implementation of artificial-parameter homotopy.

  procedure Init ( a : in Complex_Number; k : in positive );
  
  function Eval ( x : Vector; t : Complex_Number ) return Vector;
  function Diff ( x : Vector; t : Complex_Number ) return Matrix;
  function Diff ( x : Vector; t : Complex_Number ) return Vector;

end Generic_Inline_Homotopy;
