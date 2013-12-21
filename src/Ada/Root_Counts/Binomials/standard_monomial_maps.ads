with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Generic_Lists;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with Standard_Integer_VecVecs;
with Standard_Integer_Matrices;

package Standard_Monomial_Maps is

-- DESCRIPTION :
--   A monomial map of the form x(i) = c(i)*tm^v(i), for i in 1..n,
--   with c(i) a complex coefficient, tm an monomial in d parameters
--   and v(i) an array of d integer exponents, represents a solution to
--   a binomial system via an explicit monomial parametrization.

-- DATA STRUCTURES :

  type Monomial_Map ( n : integer32 ) is record
    d : integer32;                              -- dimension
    c : Standard_Complex_Vectors.Vector(1..n);  -- coefficients
    v : Standard_Integer_VecVecs.VecVec(1..n);  -- exponents
  end record;

  type Link_to_Monomial_Map is access Monomial_Map;

  type Monomial_Map_Array is
    array ( integer32 range <> ) of Link_to_Monomial_Map;
  type Link_to_Monomial_Map_Array is access Monomial_Map_Array;

  package List_of_Monomial_Maps is new Generic_LIsts(Link_to_Monomial_Map);
  type Monomial_Map_List is new List_of_Monomial_Maps.List;

  type Array_of_Monomial_Map_Lists is
    array ( integer32 range <> ) of Monomial_Map_List;
  type Link_to_Array_of_Monomial_Map_Lists is
    access Array_of_Monomial_Map_Lists;

-- CREATORS :

  function Create ( n,d : integer32;
                    c : Standard_Complex_Vectors.Vector;
                    v : Standard_Integer_VecVecs.VecVec )
                  return Monomial_Map;

  -- DESCRIPTION :
  --   Returns a monomial map in n variables, with d parameters,
  --   with coefficients in c and exponents in v.

  function Create ( maps : Monomial_Map_Array ) return Monomial_Map_List;
  function Create ( maps : Monomial_Map_List ) return Monomial_Map_Array;

  -- DESCRIPTION :
  --   Allows conversion from array to list and from list to array.
  --   Deep copies are made in both conversions.

-- SELECTORS :

  function Is_Zero ( c : Complex_Number ) return boolean;

  -- DESCRIPTION : returns true if c equals zero, false otherwise.

  function Is_One ( c : Complex_Number ) return boolean;

  -- DESCRIPTION : returns true if c equals one, false otherwise.

  function Is_Equal ( m1,m2 : Monomial_Map ) return boolean;

  -- DESCRIPTION :
  --   Returns true if m1 and m2 are equal data structures.

  function Is_In ( maps : Monomial_Map_List; m : Monomial_Map ) return boolean;

  -- DESCRIPTION :
  --   Returns true if there is a map in maps equal to m.

  function Tropisms
              ( map : Monomial_Map ) return Standard_Integer_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the tropisms defined by the exponents of the map
  --   as a vector of map.d vectors where the i-th component is defined
  --   by the i-th exponent vector of the map.

  function Tropism_Configuration 
              ( map : Monomial_Map ) return Standard_Integer_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns an d-by-n matrix of n points in d-space defined by the
  --   exponents of the monomials in the map.

  function Top_Dimension ( maps : Monomial_Map_Array ) return natural32;
  function Top_Dimension ( maps : Monomial_Map_List ) return natural32;

  -- DESCRIPTION :
  --   Returns the largest dimension of the maps.

  function Degree ( map : Monomial_Map ) return natural32;

  -- DESCRIPTION :
  --   Returns the degree of the map.

  function Degrees ( maps : Monomial_Map_List )
                   return Standard_Natural_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the degrees of each map in the list.

  function Lengths ( maps : Array_of_Monomial_Map_Lists )
                   return Standard_Natural_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the lenghts of each map in the array,
  --   so on return is a vector of maps'range.
 
-- CONSTRUCTORS :

  procedure Copy ( m1 : in Monomial_Map; m2 : out Monomial_Map );

  -- DESCRIPTION :
  --   Makes a deep copy of map m1 into map m2.

  -- REQUIRED : m2.n = m1.n.

  procedure Copy ( m1 : in Link_to_Monomial_Map;
                   m2 : out Link_to_Monomial_Map );

  -- DESCRIPTION :
  --   Makes a deep copy of the content of m1 into m2.

  procedure Append ( first,last : in out Monomial_Map_List;
                     map : in Link_to_Monomial_Map );
  procedure Append ( first,last : in out Monomial_Map_List;
                     map : in Monomial_Map );

  -- DESCRIPTION :
  --   Appends the map to the list given by its first and last element.

  procedure Concatenate ( from : in Monomial_Map_List;
                          first,last : in out Monomial_Map_List );
  procedure Concatenate ( from : in Monomial_Map_Array;
                          first,last : in out Monomial_Map_List );

  -- DESCRIPTION :
  --   Append all maps in the list from to the list given by its
  --   first and last element.  If first and last are empty,
  --   then Concatenate makes a deep copy of the list from.

-- DESTRUCTORS :

  procedure Clear ( lmm : in out Link_to_Monomial_Map );

  -- DESCRIPTION :
  --   Deallocation of the space occupied by the monomial map
  --   referred to by the link lmm.  This clear is a deep cleaning
  --   operation: all exponent vectors in lmm.v are deallocated.

  procedure Clear ( amm : in out Monomial_Map_Array );
  procedure Clear ( amm : in out Link_to_Monomial_Map_Array );

  -- DESCRIPTION :
  --   Deallocation of all entries of the array of maps amm.

  procedure Clear ( lm : in out Monomial_Map_List );

  -- DESCRIPTION :
  --   Performs a deep clear of all maps in the list lm.

  procedure Clear ( maps : in out Array_of_Monomial_Map_Lists );
  procedure Clear ( maps : in out Link_to_Array_of_Monomial_Map_Lists );

  -- DESCRIPTION :
  --   Performs a deep clear on all maps.

end Standard_Monomial_Maps;
