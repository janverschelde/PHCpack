with Standard_Complex_Vectors;
with Standard_Random_Vectors;

package body Random_Series_Generators is

  function Random_Series ( degree : integer32 ) return Series is

    cff : constant Standard_Complex_Vectors.Vector(0..degree)
        := Standard_Random_Vectors.Random_Vector(0,degree);

  begin
    return Standard_Dense_Series2.Create(cff);
  end Random_Series;

  function Random_Series ( degree : integer32 ) return Link_to_Series is

    res : constant Link_to_Series := new Series'(Random_Series(degree));

  begin
    return res;
  end Random_Series;

end Random_Series_Generators;
