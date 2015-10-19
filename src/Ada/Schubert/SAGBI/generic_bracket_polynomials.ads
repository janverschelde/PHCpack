with Abstract_Ring;
with Generic_Lists;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Bracket_Monomials;                  use Bracket_Monomials;

generic

  with package Ring is new Abstract_Ring(<>);

package Generic_Bracket_Polynomials is

-- DESCRIPTION :
--   This package offers a data representation of polynomials with
--   coefficients over any coefficient ring and brackets as unknowns.

  use Ring;

  type Bracket_Term is record
    coeff : number;
    monom : Bracket_Monomial;
  end record;

  type Bracket_Polynomial is private;

  Null_Bracket_Poly : constant Bracket_Polynomial;

-- CONSTRUCTORS :

  function Create ( m : Bracket_Monomial ) return Bracket_Polynomial;
  function Create ( t : Bracket_Term ) return Bracket_Polynomial;

  -- DESCRIPTION :
  --   Creates a bracket polynomial consisting out of one term.
  --   When only the monomial is given, the coefficient equals one.

  procedure Copy_Multiply ( t1 : in Bracket_Term; t2 : in out Bracket_Term );
  procedure Copy_Append ( t1 : in Bracket_Term; t2 : in out Bracket_Term );

  -- DESCRIPTION :
  --   Copies the term t1 to the term t2, either by multiplying the monomials
  --   (which leads to a sort), or by just appending them (no sort).

  procedure Copy ( p : in Bracket_Polynomial; q : in out Bracket_Polynomial );

  -- DESCRIPTION :
  --   Copies the polynomial p to the polynomial q.
  --   Note that q := p leads to data sharing and side effects.

-- SELECTORS :

  function Dimension ( t : Bracket_Term ) return natural32;

  -- DESCRIPTION :
  --   Returns the size of the first bracket in the bracket monomial of t.

  function Dimension ( p : Bracket_Polynomial ) return natural32;

  -- DESCRIPTION :
  --   Returns the dimension of the first bracket term of p.

-- COMPARISON OPERATIONS :

  function Is_Equal ( t1,t2 : Bracket_Term ) return boolean;
  function Is_Equal ( p,q : Bracket_Polynomial ) return boolean;

  -- DESCRIPTION :
  --   Returns true if both bracket terms/polynomials are the same.

  function ">" ( t1,t2 : Bracket_Term ) return boolean;
  function "<" ( t1,t2 : Bracket_Term ) return boolean;

  -- DESCRIPTION :
  --   The order on terms is induced by the order on bracket monomials.

-- ARITHMETIC OPERATIONS :

  function "+" ( t : Bracket_Term; p : Bracket_Polynomial )
               return Bracket_Polynomial;
  function "+" ( p : Bracket_Polynomial; t : Bracket_Term )
               return Bracket_Polynomial;
  procedure Add ( p : in out Bracket_Polynomial; t : in Bracket_Term );

  -- DESCRIPTION :
  --   Adds a term t to the bracket polynomial p.

  procedure Frontal_Add ( p : in out Bracket_Polynomial; t : in Bracket_Term );

  -- DESCRIPTION :
  --   Adds the term in front of p, without respect of order or whether
  --   a term with same degrees as t already exists.
  --   Takes a copy of t, while multiplying and thus sorting the monomials.

  procedure Frontal_Construct 
               ( p : in out Bracket_Polynomial; t : in Bracket_Term );

  -- DESCRIPTION :
  --   Adds the term in front of p, without respect of order or whether
  --   a term with same degrees as t already exists.
  --   Compared to Frontal_Add, the monomials in t will not be sorted.

  function "+" ( p,q : Bracket_Polynomial ) return Bracket_Polynomial;
  procedure Add ( p : in out Bracket_Polynomial; q : in Bracket_Polynomial );

  -- DESCRIPTION : returns p+q or makes p := p+q.

  function "-" ( t : Bracket_Term ) return Bracket_Term;
  function "-" ( p : Bracket_Polynomial ) return Bracket_Polynomial;

  procedure Min ( t : in out Bracket_Term );
  procedure Min ( p : in out Bracket_Polynomial );

  -- DESCRIPTION :
  --   Returns -t, -p or changes t, p into -t, -p.

  function "-" ( t : Bracket_Term; p : Bracket_Polynomial )
               return Bracket_Polynomial;
  function "-" ( p : Bracket_Polynomial; t : Bracket_Term )
               return Bracket_Polynomial;
  procedure Min ( p : in out Bracket_Polynomial; t : in Bracket_Term );

  -- DESCRIPTION :
  --   Subtracts a term t to the bracket polynomial p.

  procedure Frontal_Min ( p : in out Bracket_Polynomial; t : in Bracket_Term );

  -- DESCRIPTION :
  --   Adds the term -t in front of p, without respect of order or whether
  --   a term with same degrees as t already exists.

  function "-" ( p,q : Bracket_Polynomial ) return Bracket_Polynomial;
  procedure Min ( p : in out Bracket_Polynomial; q : in Bracket_Polynomial );

  -- DESCRIPTION : returns p-q or makes p := p-q.

-- ITERATORS OVER MONOMIALS :

  function Number_of_Monomials ( p : Bracket_Polynomial ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of monomials in p.

  generic
    with procedure Process ( t : in Bracket_Term; continue : out boolean );
  procedure Enumerate_Terms ( p : in Bracket_Polynomial );

  -- DESCRIPTION :
  --   Visits all terms, ordered lexicographically, and applies the
  --   procedure Process to each of them.  Enumeration stops when 
  --   continue is set to false.

-- DESTRUCTOR :

  procedure Clear ( t : in out Bracket_Term );
  procedure Clear ( p : in out Bracket_Polynomial );

  -- DESCRIPTION :
  --   Deallocates the occupied memory.

private

  package Lists_of_Bracket_Terms is new Generic_Lists(Bracket_Term);
  type Bracket_Polynomial is new Lists_of_Bracket_Terms.List;

  Null_Bracket_Poly
    : constant Bracket_Polynomial
    := Bracket_Polynomial(Lists_of_Bracket_Terms.Null_List);

end Generic_Bracket_Polynomials;
