with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Integer32_Transformations;
 use Standard_Integer32_Transformations;

package Transforming_Solutions is

-- DESCRIPTION :
--   The transformations on solutions are used in solving systems
--   by monomial transformations.

  procedure Transform ( t : in Transfo; s : in out Solution );
  function  Transform ( t : Transfo; s : Solution ) return Solution;

  procedure Transform ( t : in Transfo; L : in out Solution_List );
  function  Transform ( t : Transfo; L : Solution_List ) return Solution_List;

  -- DESCRIPTION :
  --   The transformation t will be applied.

  function Insert ( c : Complex_Number; i : integer32; s : Solution )
                  return Solution;

  procedure Insert ( c : in Complex_Number; i : in integer32; 
		     L : in out Solution_List );
  function  Insert ( c : Complex_Number; i : integer32; L : Solution_List )
                   return Solution_List;

  -- DESCRIPTION :
  --   To the ith component of the solution vector, c will be inserted.

  function Insert ( cv : vector; i : integer32; s : Solution )
		  return Solution_List;

  -- DESCRIPTION :
  --   All components in the vector cv will be inserted as the
  --   the ith component in the solution vector.

end Transforming_Solutions;
