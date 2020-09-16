with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with PentDobl_Complex_Numbers_Polar;

package body PentDobl_Complex_Series_Transforms is

  procedure Maximum_Coefficient_Modulus
              ( s : in Series;
                idx : out integer32; maxcff : out penta_double ) is

    rad : penta_double;
    minone : constant penta_double := create(-1.0);

  begin
    idx := 0;
    maxcff := minone;
    for k in 1..s.deg loop
      rad := PentDobl_Complex_Numbers_Polar.Radius(s.cff(k));
      if rad > maxcff
       then maxcff := rad; idx := k;
      end if;
    end loop;
  end Maximum_Coefficient_Modulus;

  procedure Coefficient_Modulus_Transform
              ( s : in out Series;
                idx : in integer32; maxcff : in penta_double ) is

    ddidx : constant penta_double := create(idx);
    epn : constant penta_double := 1.0/ddidx;
    fac : constant penta_double := maxcff**epn;
    divfac : penta_double := fac;

  begin
    for k in 1..s.deg loop
      s.cff(k) := s.cff(k)/divfac;
      divfac := fac*divfac;
    end loop;
  end Coefficient_Modulus_Transform;


  function Scale ( s : PentDobl_Complex_Series.Series;
                   c : penta_double )
                 return PentDobl_Complex_Series.Series is

    nc : constant Complex_Number := create(c);

  begin
    return Scale(s,nc);
  end Scale;

  function Scale ( s : PentDobl_Complex_Series.Series;
                   c : Complex_Number )
                 return PentDobl_Complex_Series.Series is

    res : PentDobl_Complex_Series.Series(s.deg);
    fac : Complex_Number := c;

  begin
    res.cff(0) := s.cff(0);
    for k in 1..s.deg loop
      res.cff(k) := fac*s.cff(k);
      fac := fac*c;
    end loop;
    return res;
  end Scale;

  procedure Transform ( s : in out Series ) is

    idx : integer32;
    maxcff : penta_double;

  begin
    Maximum_Coefficient_Modulus(s,idx,maxcff);
    if idx > 0
     then Coefficient_Modulus_Transform(s,idx,maxcff);
    end if;
  end Transform;

end PentDobl_Complex_Series_Transforms;
