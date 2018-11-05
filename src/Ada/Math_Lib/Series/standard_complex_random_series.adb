with Standard_Complex_Vectors;
with Standard_Random_Vectors;

package body Standard_Complex_Random_Series is

  function Random_Series ( degree : integer32 ) return Series is

    cff : constant Standard_Complex_Vectors.Vector(0..degree)
        := Standard_Random_Vectors.Random_Vector(0,degree);

  begin
    return Standard_Complex_Series.Create(cff);
  end Random_Series;

  function Random_Series ( degree : integer32 ) return Link_to_Series is

    res : constant Link_to_Series := new Series'(Random_Series(degree));

  begin
    return res;
  end Random_Series;

end Standard_Complex_Random_Series;
