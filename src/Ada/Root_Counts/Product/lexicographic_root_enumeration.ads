with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Vectors;

package Lexicographic_Root_Enumeration is

-- DESCRIPTION :
--   This package provides utilities to enumerate all (possible) roots
--   of a system whose polynomials are linear factors.

  generic

    with procedure process ( s : in Standard_Natural_Vectors.Vector;
                             continue : out boolean );

  procedure Lexicographic_Enumeration
              ( k,n : in natural32;
                d : in Standard_Natural_Vectors.Vector;
                acc : in out Standard_Natural_Vectors.Vector;
                continue : out boolean );

  -- DESCRIPTION :
  --   Lexicographic enumeration of all possible choices of n elements,
  --   the i-th element ranging between 1 and d(i).
  --   Each time a new choice is found, process is called.
  --   The enumeration is halted when process returns false for continue.

  -- ON ENTRY :
  --   k        current control value for the enumeration,
  --            call this initially with 1.
  --   n        total number of choices to make;
  --   d        upper bounds for the values to choose from;
  --   acc      accumulator, work space for the enumeration.

  -- ON RETURN :
  --   acc      updated accumulator;
  --   continue may be set to false by process.

  function Consecutive_Products
              ( d : Standard_Natural_Vectors.Vector )
              return Standard_Natural_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns an array r of consecutive products of d,
  --   in particular r(k) is the product of d(k+1) up to d(d'last).
  --   The output of this function is needed in the Root_Map function.

  function Root_Map ( n,k : natural32;
                      d,cp : Standard_Natural_Vectors.Vector )
                    return Standard_Natural_Vectors.Vector;
  function Root_Map ( n,k : natural32;
                      d : Standard_Natural_Vectors.Vector )
                    return Standard_Natural_Vectors.Vector;


  -- DESCRIPTION :
  --   Returns an n-vector with indices to map the k-th root onto
  --   a vector of degrees given by its consecutive products.

  -- ON ENTRY :
  --   n        number of equations in the system, d'range = 1..n;
  --   k        index to a solutions in the range 1..Product(d);
  --   d        degrees of the equations in a polynomial system;
  --   cp       consecutive products of the elements in d,
  --            i.e.: cp = Consecutive_Products(d),
  --            when omitted cp will be recomputed.

  procedure Add_One ( d : in Standard_Natural_Vectors.Vector;
                      x : in out Standard_Natural_Vectors.Vector;
                      ind : out integer32; fail : out boolean );

  -- DESCRIPTION :
  --   Gives the next position vector in the lexicographic order,
  --   immediately after x, bounded by the degrees in d.

  -- ON ENTRY :
  --   d        degrees of the equations in a system, d(i) > 0;
  --   x        current selection of solution vector.

  -- ON RETURN : 
  --   x        next vector in the lexicographic order, if not fail;
  --   ind      index in x of the entry which has increased, if not fail;
  --   fail     is true when x >= d on input.

end Lexicographic_Root_Enumeration;
