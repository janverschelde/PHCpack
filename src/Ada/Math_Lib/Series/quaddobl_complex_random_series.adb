with QuadDobl_Complex_Vectors;
with QuadDobl_Random_Vectors;

package body QuadDobl_Complex_Random_Series is

  function Random_Series ( degree : integer32 ) return Series is

    cff : constant QuadDobl_Complex_Vectors.Vector(0..degree)
        := QuadDobl_Random_Vectors.Random_Vector(0,degree);

  begin
    return QuadDobl_Complex_Series.Create(cff);
  end Random_Series;

  function Random_Series ( degree : integer32 ) return Link_to_Series is

    res : constant Link_to_Series := new Series'(Random_Series(degree));

  begin
    return res;
  end Random_Series;

end QuadDobl_Complex_Random_Series;
