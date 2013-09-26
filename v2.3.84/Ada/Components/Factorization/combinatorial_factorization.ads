with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Vectors;          use Standard_Natural_Vectors;
with Standard_Natural_VecVecs;          use Standard_Natural_VecVecs;
with Sample_Point_Lists;                use Sample_Point_Lists;

package Combinatorial_Factorization is

-- DESCRIPTION :
--   Offers the tools for a combinatorial factorization of positive
--   dimensional components of solutions.

  procedure Write ( file : in file_type; f : in VecVec );

  -- DESCRIPTION :
  --   The vectors in f represent labels of witness points,
  --   grouped along their irreducible factors.

  procedure Enumerate_Factors ( n : in natural32 );

  -- DESCRIPTION :
  --   Enumerates all possible factors of a set of n witness points.
  --   The factors are written to screen.

  generic

    with procedure Output_Factors
                     ( factors : in VecVec; count,depth : in natural32 );

    -- DESCRIPTION :
    --   This routine is called every time a new factorization is found.

    -- ON ENTRY :
    --   factors    grouping of labels of witness points;
    --   count      number of this factorization;
    --   depth      depth of the enumeration tree.

  procedure Enumerate_Factorizations ( n : in natural32 );

  -- DESCRIPTION :
  --   Enumerates all possible factorizations of a set of n witness points.
  --   The factorizations are written to screen.

  generic

    with function Is_Factor ( f : Vector ) return boolean;

    -- DESCRIPTION :
    --   This routine determines whether the vector f is a factor.

  function Search_Factorization ( n : natural32 ) return VecVec;

  -- DESCRIPTION :
  --   Search for one factorization by testing enumerated factors,
  --   calling Is_Factor for each enumerated factor.
  --   On return is the factorization found.

  generic

    with function Is_Factor ( f : Vector ) return boolean;

    -- DESCRIPTION :
    --   This routine determines whether the vector f is a factor.

  function Search_Factorization_with_Output 
             ( file : file_type; n : natural32 ) return VecVec;

  -- DESCRIPTION :
  --   Search for one factorization by testing enumerated factors,
  --   calling Is_Factor for each enumerated factor.
  --   On return is the factorization found.

  function Factor ( n : natural32;
                    grid : Array_of_Standard_Sample_Lists ) return VecVec;
  function Factor ( file : file_type; n : natural32;
                    grid : Array_of_Standard_Sample_Lists ) return VecVec;

  -- DESCRIPTION :
  --   Uses the combinatorial enumeration of factors and linear traces
  --   to determine the factorization of a solution component.

  -- ON ENTRY :
  --   file       for intermediate output and diagnostics;
  --   n          #witness points, sum of degrees over all factors;
  --   grid       sample grid of n witness points, of range 0..2.

  -- ON RETURN :
  --   groups of labels of witness points for each irreducible factor.

end Combinatorial_Factorization;
