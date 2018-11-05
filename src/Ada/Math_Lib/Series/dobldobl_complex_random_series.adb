with DoblDobl_Complex_Vectors;
with DoblDobl_Random_Vectors;

package body DoblDobl_Complex_Random_Series is

  function Random_Series ( degree : integer32 ) return Series is

    cff : constant DoblDobl_Complex_Vectors.Vector(0..degree)
        := DoblDobl_Random_Vectors.Random_Vector(0,degree);

  begin
    return DoblDobl_Complex_Series.Create(cff);
  end Random_Series;

  function Random_Series ( degree : integer32 ) return Link_to_Series is

    res : constant Link_to_Series := new Series'(Random_Series(degree));

  begin
    return res;
  end Random_Series;

end DoblDobl_Complex_Random_Series;
