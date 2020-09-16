with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with OctoDobl_Complex_Numbers_Polar;

package body OctoDobl_Complex_Series_Transforms is

  procedure Maximum_Coefficient_Modulus
              ( s : in Series;
                idx : out integer32; maxcff : out octo_double ) is

    rad : octo_double;
    minone : constant octo_double := create(-1.0);

  begin
    idx := 0;
    maxcff := minone;
    for k in 1..s.deg loop
      rad := OctoDobl_Complex_Numbers_Polar.Radius(s.cff(k));
      if rad > maxcff
       then maxcff := rad; idx := k;
      end if;
    end loop;
  end Maximum_Coefficient_Modulus;

  procedure Coefficient_Modulus_Transform
              ( s : in out Series;
                idx : in integer32; maxcff : in octo_double ) is

    ddidx : constant octo_double := create(idx);
    epn : constant octo_double := 1.0/ddidx;
    fac : constant octo_double := maxcff**epn;
    divfac : octo_double := fac;

  begin
    for k in 1..s.deg loop
      s.cff(k) := s.cff(k)/divfac;
      divfac := fac*divfac;
    end loop;
  end Coefficient_Modulus_Transform;


  function Scale ( s : OctoDobl_Complex_Series.Series;
                   c : octo_double )
                 return OctoDobl_Complex_Series.Series is

    nc : constant Complex_Number := create(c);

  begin
    return Scale(s,nc);
  end Scale;

  function Scale ( s : OctoDobl_Complex_Series.Series;
                   c : Complex_Number )
                 return OctoDobl_Complex_Series.Series is

    res : OctoDobl_Complex_Series.Series(s.deg);
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
    maxcff : octo_double;

  begin
    Maximum_Coefficient_Modulus(s,idx,maxcff);
    if idx > 0
     then Coefficient_Modulus_Transform(s,idx,maxcff);
    end if;
  end Transform;

end OctoDobl_Complex_Series_Transforms;
