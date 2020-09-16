with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with DecaDobl_Complex_Numbers_Polar;

package body DecaDobl_Complex_Series_Transforms is

  procedure Maximum_Coefficient_Modulus
              ( s : in Series;
                idx : out integer32; maxcff : out deca_double ) is

    rad : deca_double;
    minone : constant deca_double := create(-1.0);

  begin
    idx := 0;
    maxcff := minone;
    for k in 1..s.deg loop
      rad := DecaDobl_Complex_Numbers_Polar.Radius(s.cff(k));
      if rad > maxcff
       then maxcff := rad; idx := k;
      end if;
    end loop;
  end Maximum_Coefficient_Modulus;

  procedure Coefficient_Modulus_Transform
              ( s : in out Series;
                idx : in integer32; maxcff : in deca_double ) is

    ddidx : constant deca_double := create(idx);
    epn : constant deca_double := 1.0/ddidx;
    fac : constant deca_double := maxcff**epn;
    divfac : deca_double := fac;

  begin
    for k in 1..s.deg loop
      s.cff(k) := s.cff(k)/divfac;
      divfac := fac*divfac;
    end loop;
  end Coefficient_Modulus_Transform;


  function Scale ( s : DecaDobl_Complex_Series.Series;
                   c : deca_double )
                 return DecaDobl_Complex_Series.Series is

    nc : constant Complex_Number := create(c);

  begin
    return Scale(s,nc);
  end Scale;

  function Scale ( s : DecaDobl_Complex_Series.Series;
                   c : Complex_Number )
                 return DecaDobl_Complex_Series.Series is

    res : DecaDobl_Complex_Series.Series(s.deg);
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
    maxcff : deca_double;

  begin
    Maximum_Coefficient_Modulus(s,idx,maxcff);
    if idx > 0
     then Coefficient_Modulus_Transform(s,idx,maxcff);
    end if;
  end Transform;

end DecaDobl_Complex_Series_Transforms;
