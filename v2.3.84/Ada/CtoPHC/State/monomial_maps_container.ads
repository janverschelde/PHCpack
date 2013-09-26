with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;
with Standard_Monomial_Maps;            use Standard_Monomial_Maps;

package Monomial_Maps_Container is

-- DESCRIPTION :
--   This package provides a container for monomial maps
--   designed for the interface with C.

  procedure Initialize ( sols : in Array_of_Monomial_Map_Lists ); 

  -- DESCRIPTION :
  --   Initializes the container with a list of monomial maps.

  function Retrieve return Link_to_Array_of_Monomial_Map_Lists;

  -- DESCRIPTION :
  --   Returns the lists of monomial maps stored in the container.

  function Top_Dimension return integer32;

  -- DESCRIPTION :
  --   Returns the top dimension of the monomial maps in the container,
  --   returns -1 if there are no maps stored.

  function Number_of_Maps ( d : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of maps of dimension d stored,
  --   returns -1 if there are no maps stored in the container.

  function Retrieve_Map ( d,k : integer32 ) return Link_to_Monomial_Map;

  -- DESCRIPTION :
  --   Returns the k-th map of dimension d,
  --   or the null pointer if there is no such map.

  function Degree ( d,k : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Returns the degree of the k-th map of dimension d.
  --   Returns -1 if there are no maps stored in the container,
  --   or if there are no maps of that dimension.

  function Coefficients
             ( d,k : integer32 ) return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the coefficients of the k-th map of dimension d.
  --   If there are no maps stored in the container,
  --   or if there are no maps of that dimension,
  --   then the vector on return has the range 0..0.

  function Exponents
             ( d,k : integer32 ) return Standard_Integer_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the exponents of the k-th map of dimension d.
  --   If there are no maps stored in the container,
  --   or if there are no maps of that dimension,
  --   then the vector on return has the range 0..0.

  procedure Coefficients_and_Exponents
              ( d,k : in integer32;
                cff : out Standard_Complex_Vectors.Vector;
                exp : out Standard_Integer_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Returns the coefficients and exponents of the k-th map of dimension d.
  --   If there are no maps stored in the container,
  --   or if there are no maps of that dimension,
  --   then no assignments are made to cff and exp.

  -- REQUIRED : cff'range = exp'range = 1..n,
  --   where n is the ambient dimension.

  procedure Clear;

  -- DESCRIPTION :
  --   Destroys the stored monomial maps.

end Monomial_Maps_Container;
