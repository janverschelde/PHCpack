with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural64_VecVecs;       use Standard_Natural64_VecVecs;

package Monomial_Hashing is

-- DESCRIPTION :
--   This package provides enumerators for monomials and hashing functions
--   to deal efficiently with derivative operators.

  generic
    with procedure Monomial ( m : in Standard_Natural_Vectors.Vector );
  procedure Enumerate_Monomials ( d,n : in natural32 );

  -- DESCRIPTION :
  --   Enumerates all monomials of degree d in n variables,
  --   each time a new monomial m is found, Monomial(m) is called.

  generic
    with procedure Monomial ( m : in Standard_Natural_Vectors.Vector );
  procedure Enumerate_Leaves_of_Monomial_Tree ( d,n : in natural32 );

  -- DESCRIPTION :
  --   Enumerates all monomials of degree d in n variables, as they occur
  --   as leaves in the tree of enumerating all possible derivations,
  --   each time a new monomial m is found, Monomial(m) is called.

  function Monomial_Count ( d,n : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of monomials of degree d in n variables.

  function Monomial_Code
             ( d : natural32; m : Standard_Natural_Vectors.Vector )
             return natural64;

  -- DESCRIPTION :
  --   Produces a unique number to code the monomial m, using the
  --   numbers in m as the digits in a number system of basis d.

  function Monomial_Keys
             ( k,n : natural32 ) return Standard_Natural64_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns a vecvec of range 1..k with all the codes for all monomials
  --   indexing the application of derivatives to n variables up to order k.
  --   This is a hash table to index all derivative operators.

  function Search ( monkeys : Standard_Natural64_VecVecs.VecVec;
                    c : natural64; s : natural32 ) return natural32;
  function Search ( monkeys : Standard_Natural64_VecVecs.VecVec;
                    m : Standard_Natural_Vectors.Vector; s : natural32 )
                  return natural32;
  function Search ( monkeys : Standard_Natural64_VecVecs.VecVec;
                    m : Standard_Natural_Vectors.Vector ) return natural32;

  -- DESCRIPTION :
  --   Returns the index in the hash table of the monomial code,
  --   returns 0 if the code has not been found.

  -- ON ENTRY :
  --   monkeys    hash table produced by the function "Monomial_Keys";
  --   c          code of a monomial m;
  --   s          sum of the degrees in the monomial,
  --              if s > monkeys'last, then zero will be returned;
  --   m          monomial in several variables.

end Monomial_Hashing;
