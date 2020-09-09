with OctoDobl_Complex_Vectors;
with OctoDobl_Random_Vectors;

package body OctoDobl_Complex_Random_Series is

  function Random_Series ( degree : integer32 ) return Series is

    cff : constant OctoDobl_Complex_Vectors.Vector(0..degree)
        := OctoDobl_Random_Vectors.Random_Vector(0,degree);

  begin
    return OctoDobl_Complex_Series.Create(cff);
  end Random_Series;

  function Random_Series ( degree : integer32 ) return Link_to_Series is

    res : constant Link_to_Series := new Series'(Random_Series(degree));

  begin
    return res;
  end Random_Series;

end OctoDobl_Complex_Random_Series;
