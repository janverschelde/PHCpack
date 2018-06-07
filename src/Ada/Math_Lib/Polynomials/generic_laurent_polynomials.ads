with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Generic_Lists;
with Abstract_Ring;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;

generic

  with package Ring is new Abstract_Ring(<>);

package Generic_Laurent_Polynomials is

-- DESCRIPTION :
--   This package represents Laurent polynomials in several variables with
--   coefficients over any ring, to be specified by instantiation.
--   The exponents can be negative.

  use Ring;

-- DATA STRUCTURES :

  type Degrees is new Standard_Integer_Vectors.Link_to_Vector;

  type Term is record
    cf : number;     -- coefficient of the term
    dg : Degrees;    -- the degrees of the term
  end record;

  type Poly is private;

  Null_Poly : constant Poly;    -- represents zero in the polynomial ring
  One_Poly : constant Poly;     -- represents one in the polynomial ring

-- CONSTRUCTORS :

  function Create ( i : integer ) return Poly;
  function Create ( n : number ) return Poly;

  function Create ( t : Term ) return Poly;

  procedure Copy ( t1 : in Term; t2 : in out Term );    -- makes a deep copy
  procedure Copy ( p: in Poly; q : in out Poly );

-- SELECTORS :

  function Equal ( t1,t2 : Term )  return boolean;
  function Equal ( p,q : Poly )  return boolean;

  function Number_of_Unknowns ( p : Poly ) return natural32;
  function Number_of_Terms    ( p : Poly ) return natural32;

  function Size_of_Support ( t : Term ) return natural32;
  -- returns number of variables with nonzero exponent
  function Size_of_Support ( p : Poly ) return natural32;
  -- returns number of variables with nonzero exponent

  function Variables_in_Support
             ( t : Term ) return Standard_Natural_Vectors.Vector;
  -- returns 0/1 vector with 1 at place of a variable with nonzero exponent
  function Variables_in_Support
             ( p : Poly ) return Standard_Natural_Vectors.Vector;
  -- returns 0/1 vector with 1 at place of a variable with nonzero exponent

  function Degree ( p : Poly ) return integer32;            -- return deg(p);

  function Maximal_Degree ( p : Poly; i : integer32 ) return integer32;
             -- returns maximal degree of xi in p;
  function Maximal_Degrees ( p : Poly ) return Degrees;
             -- Maximal_Degrees(p)(i) = Maximal_Degree(p,i)
  function Minimal_Degree ( p : Poly; i : integer32 ) return integer32;
             -- returns minimal degree of xi in p;
  function Minimal_Degrees ( p : Poly ) return Degrees;
             -- Minimal_Degrees(p)(i) = Minimal_Degree(p,i)

  function "<" ( d1,d2 : Degrees ) return boolean;          -- return d1 < d2
  function ">" ( d1,d2 : Degrees ) return boolean;          -- return d1 > d2

  function Coeff ( p : Poly; d : Degrees ) return number;
   -- Ex.: Coeff(c1*x^2+c2*x*y^3,(1 2))=c2;  Coeff(c1*x^2+c2,(1 0))=zero;

  function Head ( p : Poly ) return Term;               -- return head term

-- ARITHMETICAL OPERATIONS :

  function "+" ( p : Poly; t : Term ) return Poly;      -- return p+t;
  function "+" ( t : Term; p : Poly ) return Poly;      -- return t+p;
  function "+" ( p : Poly ) return Poly;                -- returns copy of p;
  function "+" ( p,q : Poly ) return Poly;              -- return p+q;
  function "-" ( p : Poly; t : Term ) return Poly;      -- return p-t;
  function "-" ( t : Term; p : Poly ) return Poly;      -- return t-p;
  function "-" ( p : Poly ) return Poly;                -- return -p;
  function "-" ( p,q : Poly ) return Poly;              -- return p-q;
  function "*" ( p : Poly; a : number ) return Poly;    -- return a*p;
  function "*" ( a : number; p : Poly ) return Poly;    -- return p*a;
  function "*" ( p : Poly; t : Term ) return Poly;      -- return p*t;
  function "*" ( t : Term; p : Poly ) return Poly;      -- return t*p;
  function "*" ( p,q : Poly ) return Poly;              -- return p*q;
  function "**" ( t : Term; k : natural32 ) return Term;  -- return t**k;
  function "**" ( p : Poly; k : natural32 ) return Poly;  -- return p**k;

  procedure Add ( p : in out Poly; t : in Term );       -- p := p + t;
  procedure Add ( p : in out Poly; q : in Poly );       -- p := p + q;
  procedure Sub ( p : in out Poly; t : in Term );       -- p := p - t;
  procedure Min ( p : in out Poly );                    -- p := -p;
  procedure Sub ( p : in out Poly; q : in Poly );       -- p := p - q;
  procedure Mul ( p : in out Poly; a : in number );     -- p := p * a;
  procedure Mul ( p : in out Poly; t : in Term );       -- p := p * t;
  procedure Mul ( p : in out Poly; q : in Poly );       -- p := p * q;
  procedure Pow ( t : in out Term; k : in natural32 );  -- t := t**k;
  procedure Pow ( p : in out Poly; k : in natural32 );  -- p := p**k;

  function  Diff ( p : Poly; i : integer32 ) return Poly; 
  procedure Diff ( p : in out Poly; i : in integer32 );
    -- symbolic differentiation w.r.t. the i-th unknown of p

-- ITERATORS : run through all terms of p and apply the generic procedure.

  generic
    with procedure process ( t : in out Term; continue : out boolean );
  procedure Changing_Iterator ( p : in out Poly );  -- t can be changed
  generic
    with procedure process ( t : in Term; continue : out boolean );
  procedure Visiting_Iterator ( p : in Poly );      -- t can only be read

-- DESTRUCTORS : deallocate memory.

  procedure Clear ( t : in out Term );
  procedure Clear ( p : in out Poly );

private

  package Term_List is new Generic_Lists(Term);
  type Poly_Rep is new Term_List.List;

  type Poly is access Poly_Rep;

  Null_Poly : constant Poly := null;

  One_Term : constant Term := (one,null);
  One_Poly_Rep : constant Poly_Rep := Create(One_Term);
  One_Poly : constant Poly := new Poly_Rep'(One_Poly_Rep);

end Generic_Laurent_Polynomials;
