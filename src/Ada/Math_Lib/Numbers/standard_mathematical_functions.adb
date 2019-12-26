with Ada.Numerics.Generic_Elementary_Functions;
with Ada.Numerics;

package body Standard_Mathematical_Functions is

  package Double_Elementary_Functions is
    new Ada.Numerics.Generic_Elementary_Functions (double_float);

-- ADDED FOR PATCH around bug in gcc version 3.2.2

 -- function C_COS ( x : double_float ) return double_float is

 --   function cos ( x : double_float ) return double_float;
 --   pragma interface(C,cos);

 -- begin
 --   return cos(x);
 -- end C_COS;

 -- function C_SIN ( x : double_float ) return double_float is

 --   function sin ( x : double_float ) return double_float;
 --   pragma interface(C,sin);

 -- begin
 --   return sin(x);
 -- end C_SIN;

-- END OF ADDITION FOR PATCH

  function EXP ( x : double_float ) return double_float is
  begin
    return Double_Elementary_Functions.exp(x);
  end EXP;

  function LN ( x : double_float ) return double_float is
  begin
    return Double_Elementary_Functions.log(x);
  end LN;

  function "**" ( x,y : double_float ) return double_float is
  begin
    return Double_Elementary_Functions."**"(x,y);
  end "**";

  function LOG2 ( x : double_float ) return double_float is
  begin
    return Double_Elementary_Functions.LOG(x,2.0);
  end LOG2;

  function LOG10 ( x : double_float ) return double_float is
  begin
    return Double_Elementary_Functions.LOG(x,10.0);
  end LOG10;

  function SQRT ( x : double_float ) return double_float is
  begin
    return Double_Elementary_Functions.SQRT(x);
  end SQRT;

  function SIN ( x : double_float ) return double_float is
  begin
    return Double_Elementary_Functions.SIN(x); -- PATCH for gcc 3.2.2
   -- return C_SIN(x);
  end SIN;

  function COS ( x : double_float ) return double_float is
  begin
    return Double_Elementary_Functions.COS(x); -- PATCH for gcc 3.2.2
   -- return C_COS(x);
  end COS;

  function TAN ( x : double_float ) return double_float is
  begin
    return Double_Elementary_Functions.TAN(x);
  end TAN;

  function ARCSIN ( x : double_float ) return double_float is
  begin
    return Double_Elementary_Functions.ARCSIN(x);
  end ARCSIN;

  function ARCCOS ( x : double_float ) return double_float is
  begin
    return Double_Elementary_Functions.ARCCOS(x);
  end ARCCOS;

  function ARCTAN ( x : double_float ) return double_float is
  begin
    return Double_Elementary_Functions.ARCTAN(x);
  end ARCTAN;

  function Radius ( x,y : double_float ) return double_float is
  begin
    return Sqrt(x**2+y**2);
  end Radius;

  function Angle ( x,y : double_float ) return double_float is
  begin
    return Double_Elementary_Functions.ARCTAN(y=>x,x=>y);
  end Angle;

end Standard_Mathematical_Functions;
