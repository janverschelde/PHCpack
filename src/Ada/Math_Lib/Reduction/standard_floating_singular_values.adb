with Standard_Complex_Singular_Values;

package body Standard_Floating_Singular_Values is

  function Min0 ( a,b : integer32 ) return integer32 is
  begin
    return Standard_Complex_Singular_Values.Min0(a,b);
  end Min0;

  procedure SVD ( x : in out Matrix; n,p : in integer32;
                  s,e : out Vector; u : out Matrix; v : out Matrix;
                  job : in integer32; info : out integer32 ) is
  begin
    null;
  end SVD;

  function Rank ( s : Vector ) return natural is
  begin
    return 0;
  end Rank;

  function Rank ( s : Vector; tol : double_float ) return natural is
  begin
    return 0;
  end Rank;

  function Inverse_Condition_Number
             ( s : Vector ) return double_float is
  begin
    return 0.0;
  end Inverse_Condition_Number;

  function Inverse_Condition_Number 
             ( s : Vector; tol : double_float ) return double_float;
  begin
    return 0.0;
  end Inverse_Condition_Number;

  function Transpose ( z : Matrix ) return Matrix is
  begin
    return z;
  end Transpose;

  function Inverse ( u,v : Matrix; s : Vector ) return Matrix is
  begin
    return u;
  end Inverse;

  function Inverse ( u,v : Matrix; s : Vector; tol : double_float )
                   return Matrix is
  begin
    return u;
  end Inverse;

  function Solve ( u,v : Matrix; s,b : Vector ) return Vector is
  begin
    return s;
  end Solve;

  function Solve ( u,v : Matrix; s,b : Vector; tol : double_float ) 
                 return Vector is
  begin
    return s;
  end Solve;

end Standard_Floating_Singular_Values;
