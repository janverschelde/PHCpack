package body Generic_Inline_Homotopy is

  procedure Init ( a : in Complex_Number; k : in positive ) is
  begin
    Homotopy_Constants(a,k);
  end Init;
  
  function Eval ( x : Vector; t : Complex_Number ) return Vector is
  begin
    return Eval_Homotopy(x,t);
  end Eval;

  function Diff ( x : Vector; t : Complex_Number ) return Matrix is
  begin
    return Diff_Homotopy(x,t);
  end Diff;

  function Diff ( x : Vector; t : Complex_Number ) return Vector is
  begin
    return Diff_Homotopy(x,t);
  end Diff;

end Generic_Inline_Homotopy;
