with Standard_Mathematical_Functions;
with Standard_Complex_Numbers_Polar;

package body Standard_Complex_Series_Transforms is

  procedure Maximum_Coefficient_Modulus
              ( s : in Series;
                idx : out integer32; maxcff : out double_float ) is

    rad : double_float;

  begin
    idx := 0;
    maxcff := -1.0;
    for k in 1..s.deg loop
      rad := Standard_Complex_Numbers_Polar.Radius(s.cff(k));
      if rad > maxcff
       then maxcff := rad; idx := k;
      end if;
    end loop;
  end Maximum_Coefficient_Modulus;

  procedure Coefficient_Modulus_Transform
              ( s : in out Series;
                idx : in integer32; maxcff : in double_float ) is

    use Standard_Mathematical_Functions;

    epn : constant double_float := 1.0/double_float(idx);
    fac : constant double_float := maxcff**epn;
    divfac : double_float := fac;

  begin
    for k in 1..s.deg loop
      s.cff(k) := s.cff(k)/divfac;
      divfac := fac*divfac;
    end loop;
  end Coefficient_Modulus_Transform;


  function Scale ( s : Standard_Complex_Series.Series;
                   c : double_float )
                 return Standard_Complex_Series.Series is

    nc : constant Complex_Number := create(c);

  begin
    return Scale(s,nc);
  end Scale;

  function Scale ( s : Standard_Complex_Series.Series;
                   c : Complex_Number )
                 return Standard_Complex_Series.Series is

    res : Standard_Complex_Series.Series(s.deg);
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
    maxcff : double_float;

  begin
    Maximum_Coefficient_Modulus(s,idx,maxcff);
    if idx > 0
     then Coefficient_Modulus_Transform(s,idx,maxcff);
    end if;
  end Transform;

end Standard_Complex_Series_Transforms;
