with Standard_Complex_Vectors;
with Standard_Random_Vectors;

package body Standard_Random_Series is

-- DESCRIPTION :
--   Exports functions that return random power series,
--   truncated to the given order.

  function Random_Series ( order : integer32 ) return Series is

    cff : Standard_Complex_Vectors.Vector(0..order)
        := Standard_Random_Vectors.Random_Vector(0,order);

  begin
    return Create(cff);
  end Random_Series;

end Standard_Random_Series;
