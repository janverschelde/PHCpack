with HexaDobl_Complex_Vectors;
with HexaDobl_Random_Vectors;

package body HexaDobl_Complex_Random_Series is

  function Random_Series ( degree : integer32 ) return Series is

    cff : constant HexaDobl_Complex_Vectors.Vector(0..degree)
        := HexaDobl_Random_Vectors.Random_Vector(0,degree);

  begin
    return HexaDobl_Complex_Series.Create(cff);
  end Random_Series;

  function Random_Series ( degree : integer32 ) return Link_to_Series is

    res : constant Link_to_Series := new Series'(Random_Series(degree));

  begin
    return res;
  end Random_Series;

end HexaDobl_Complex_Random_Series;
