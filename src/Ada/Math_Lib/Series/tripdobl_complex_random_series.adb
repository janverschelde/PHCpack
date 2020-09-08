with TripDobl_Complex_Vectors;
with TripDobl_Random_Vectors;

package body TripDobl_Complex_Random_Series is

  function Random_Series ( degree : integer32 ) return Series is

    cff : constant TripDobl_Complex_Vectors.Vector(0..degree)
        := TripDobl_Random_Vectors.Random_Vector(0,degree);

  begin
    return TripDobl_Complex_Series.Create(cff);
  end Random_Series;

  function Random_Series ( degree : integer32 ) return Link_to_Series is

    res : constant Link_to_Series := new Series'(Random_Series(degree));

  begin
    return res;
  end Random_Series;

end TripDobl_Complex_Random_Series;
