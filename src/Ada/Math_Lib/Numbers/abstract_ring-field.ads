generic

  with function "<" ( a,b : number ) return boolean;
  with function ">" ( a,b : number ) return boolean;

  with function AbsVal ( a : number ) return number;           -- returns |a|

  with function "/" ( a,b : number ) return number;

  with procedure Div ( a : in out number; b : in number );     -- a := a/b

package Abstract_Ring.Field is end;

-- DESCRIPTION :
--   Defines a field as an extension of a ring with an order and division.
