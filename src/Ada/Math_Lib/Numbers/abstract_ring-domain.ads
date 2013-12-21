generic

  with function "<" ( a,b : number ) return boolean;
  with function ">" ( a,b : number ) return boolean;

  with function "/" ( a,b : number ) return number;

  with function  Rmd ( a,b : number ) return number;            -- a mod b
  with procedure Rmd ( a : in out number; b : in number );      -- a := a mod b

  with procedure Div ( a : in out number; b : in number );      -- a := a/b
  with procedure Div ( a,b : in number;                         -- a = b*q+r
                       q : out number; r : out number );        -- q := a/b
  with procedure Div ( a : in out number; b : in number;        -- a := a/b
                       r : out number );                        -- r := a mod b

package Abstract_Ring.Domain is end;

-- DESCRIPTION :
--   An Euclidean domain is a ring where for any pair (a,b), b /= 0,
--   there exists a quotient q and a remainder r, such that a = b*q + r.
--   Note that here we only have the remainder and division operations.
--   Any package that makes usage of the full euclidean domain must as
--   well have the ring, instantiated with the same number type.
