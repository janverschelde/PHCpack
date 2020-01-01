with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with QuadDobl_Complex_Numbers_Polar;

package body QuadDobl_Complex_Series_Transforms is

  procedure Maximum_Coefficient_Modulus
              ( s : in Series;
                idx : out integer32; maxcff : out quad_double ) is

    rad : quad_double;
    minone : constant quad_double := create(-1.0);

  begin
    idx := 0;
    maxcff := minone;
    for k in 1..s.deg loop
      rad := QuadDobl_Complex_Numbers_Polar.Radius(s.cff(k));
      if rad > maxcff
       then maxcff := rad; idx := k;
      end if;
    end loop;
  end Maximum_Coefficient_Modulus;

  procedure Coefficient_Modulus_Transform
              ( s : in out Series;
                idx : in integer32; maxcff : in quad_double ) is

    ddidx : constant quad_double := create(idx);
    epn : constant quad_double := 1.0/ddidx;
    fac : constant quad_double := maxcff**epn;
    divfac : quad_double := fac;

  begin
    for k in 1..s.deg loop
      s.cff(k) := s.cff(k)/divfac;
      divfac := fac*divfac;
    end loop;
  end Coefficient_Modulus_Transform;


  function Scale ( s : QuadDobl_Complex_Series.Series;
                   c : quad_double )
                 return QuadDobl_Complex_Series.Series is

    nc : constant Complex_Number := create(c);

  begin
    return Scale(s,nc);
  end Scale;

  function Scale ( s : QuadDobl_Complex_Series.Series;
                   c : Complex_Number )
                 return QuadDobl_Complex_Series.Series is

    res : QuadDobl_Complex_Series.Series(s.deg);
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
    maxcff : quad_double;

  begin
    Maximum_Coefficient_Modulus(s,idx,maxcff);
    if idx > 0
     then Coefficient_Modulus_Transform(s,idx,maxcff);
    end if;
  end Transform;

end QuadDobl_Complex_Series_Transforms;
