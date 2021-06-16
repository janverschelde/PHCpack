with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Complex_Solutions;          use Standard_Complex_Solutions;
with C_Integer_Arrays;                    use C_Integer_Arrays;
with C_Double_Arrays;                     use C_Double_Arrays;

package Coefficient_Solution_Vectors is

-- DESCRIPTION :
--   This package provides conversion routines to export solutions
--   to C programs.

-- PART I : converting to the coefficient vector representation

  function Multiplicity ( s : Solution ) return integer32;

  -- DESCRIPTION :
  --   Returns the multiplicity information of the solution.

  function Multiplicities ( s : Solution_List ) return C_Integer_Array;

  -- DESCRIPTION :
  --   Returns all multiplicities of the solutions in the list.

  function Coefficients ( s : Solution ) return C_Double_Array;

  -- DESCRIPTION :
  --   Returns the values for the continuation parameter, the solution
  --   vector and the (err,rco,res) status of Newton's method as one
  --   long vector of double floats.

  function Coefficients ( s : Solution_List ) return C_Double_Array;

  -- DESCRIPTION :
  --   Returns the coefficients representation for each solution,
  --   stacked into one long vector of double floats.

  -- REQUIRED : not Is_Null(s).

-- PART I bis : concatenating the coefficients representation

  function Concat ( n : natural32; m : C_Integer_Array; c : C_Double_Array )
                  return C_Double_Array;

  -- DESCRIPTION :
  --   Concatenates the dimension n, multiplicies m, and coefficients c
  --   into one long C double array.

  function Dimension ( x : C_Double_Array ) return natural32;

  -- DESCRIPTION :
  --   Returns the dimension of the solution vectors.

  function Multiplicities ( x : C_Double_Array ) return C_Integer_Array;

  -- DESCRIPTION :
  --   Returns the multiplicities for every solution.

  function Coefficients ( x : C_Double_Array ) return C_Double_Array;

  -- DESCRIPTION :
  --   Returns the coefficients for every solution.

-- PART II : creating solutions from the coefficient vector representation

  function Create ( n,m : natural32; c : C_Double_Array ) return Solution;

  -- DESCRIPTION :
  --   Returns an n-dimensional solution of multiplicity m from the
  --   coefficient vector representation in c.

  function Create ( n : natural32; m : C_Integer_Array; c : C_Double_Array )
                  return Solution_List;

  -- DESCRIPTION :
  --   Returns a solution list of n-dimensional solutions of multiplicities
  --   in m and with coefficients in c.

end Coefficient_Solution_Vectors;
