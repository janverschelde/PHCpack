with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Series_Vectors;

package Standard_Pade_Approximants is

-- DESCRIPTION :
--   Defines a data structure to represent Pade approximants
--   and vectors of Pade approximants, constructed from a series,
--   in standard double precision.

-- DATA STRUCTURES :

  type Pade is private;
  -- A Pade approximant is the quotient of two polynomials in one variable.
  -- The coefficient vectors are stored for numerator and denominator.

  type Pade_Vector is array ( integer32 range <> ) of Pade;
  -- A Pade vector is a vector of Pade approximants.

  type Link_to_Pade_Vector is access Pade_Vector;
  type Pade_VecVec is array ( integer32 range <> ) of Link_to_Pade_Vector;

-- CONSTRUCTORS :

  function Create ( num,den : Standard_Complex_Vectors.Vector )
                  return Pade;

  -- DESCRIPTION :
  --   Returns the Pade approximant defined by the coefficient
  --   vectors of the numerator and denominator,
  --   respectively given in num and den.

  -- REQUIRED : num'first = 0 = den'first.

  function Coefficients ( srv : Standard_Complex_Series_Vectors.Vector;
                          idx : integer32 )
                        return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the coefficients of series srv at the compenent
  --   with index idx.

  -- REQUIRED : idx in srv'range.

  function Create ( numdeg,dendeg : integer32;
                    srv : Standard_Complex_Series_Vectors.Vector;
                    verbose : boolean := false )
                  return Pade_Vector;

  -- DESCRIPTION :
  --   Returns a vector of Pade approximants, using the coefficients
  --   of the series vector given in srv.
  --   The degree of the numerator of each Pade approximant equals numdeg
  --   and the degree of the denominator equals dendeg.

  -- REQUIRED : each series has at least numdeg+dendeg+1 coefficients.

  function Allocate ( numdeg,dendeg : integer32 ) return Pade;

  -- DESCRIPTION :
  --   Returns an allocated Pade approximant with coefficients of
  --   numerator to degree numdeg and denominator to degree dendeg.
  --   The coefficient vectors are initialized to zero.

  function Allocate ( dim,numdeg,dendeg : integer32 ) return Pade_Vector;

  -- DESCRIPTION :
  --   Returns a vector of Pade approximants of dimension dim,
  --   allocated with coefficients of each numerator to degree numdeg
  --   and coefficients of each denominator to degree dendeg.
  --   The coefficient vectors are initialized to zero.

  procedure Create ( pv : in out Pade_Vector;
                     srv : in Standard_Complex_Series_Vectors.Vector;
                     verbose : in boolean := false );

  -- DESCRIPTION :
  --   Given an allocated Pade vector in pv,
  --   defines the coefficients of the Pade approximants using
  --   the vector of power series in srv.

  -- REQUIRED :
  --   The Pade approximants all have the same degrees 
  --   and each series has sufficiently many coefficients.

-- SELECTORS :

  function Numerator_Degree ( p : Pade ) return integer32;

  -- DESCRIPTION :
  --   Returns the degree of the numerator or -1 if p is undefined.

  function Denominator_Degree ( p : Pade ) return integer32;

  -- DESCRIPTION :
  --   Returns the degree of the denominator or -1 if p is undefined.

  function Numerator_Coefficients
             ( p : Pade ) return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 0..Numerator_Degree(p) of the
  --   complex coefficients of the Pade approximant of the numerator.

  -- REQUIRED : Numerator_Degree(p) >= 0.

  function Denominator_Coefficients
             ( p : Pade ) return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 0..Denominator_Degree(p) of the
  --   complex coefficients of the Pade approximant of the denominator.

  -- REQUIRED : Denominator_Degree(p) >= 0.

-- EVALUATORS :

  function Eval ( p : Pade; x : double_float ) return Complex_Number;
  function Eval ( p : Pade; x : Complex_Number ) return Complex_Number;

  -- DESCRIPTION :
  --   Evaluates the Pade approximant at x.

  -- REQUIRED : p is well defined.

  function Eval ( p : Pade_Vector; x : double_float )
                return Standard_Complex_Vectors.Vector;
  function Eval ( p : Pade_Vector; x : Complex_Number )
                return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Evaluates the vector of Pade approximants at x.
  --   The vector on return has the same range as p.

  -- REQUIRED : all components of p are well defined.

-- DESTRUCTORS :

  procedure Clear ( p : in out Pade );
  procedure Clear ( p : in out Pade_Vector );
  procedure Clear ( p : in out Link_to_Pade_Vector );
  procedure Clear ( p : in out Pade_VecVec );

  -- DESCRIPTION :
  --   Deallocates the memory occupied by p.

private

  type Pade_Rep ( numdeg,dendeg : integer32 ) is record
    num : Standard_Complex_Vectors.Vector(0..numdeg); -- numerator
    den : Standard_Complex_Vectors.Vector(0..dendeg); -- denominator
  end record;

  type Pade is access Pade_Rep;

end Standard_Pade_Approximants;
