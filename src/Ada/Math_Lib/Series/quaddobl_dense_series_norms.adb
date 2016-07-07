with QuadDobl_Complex_Numbers;
with QuadDobl_Mathematical_Functions;
with QuadDobl_Complex_Numbers_Polar;
with QuadDobl_Algebraic_Series;

package body QuadDobl_Dense_Series_Norms is

  function Norm ( s : Series ) return Series is

    c : constant Series := Conjugate(s);
    p : constant Series := c*s;
    res : Series;

  begin
    res := QuadDobl_Algebraic_Series.sqrt(p,0);
    return res;
  end Norm;

  procedure Normalize ( s : in out Series ) is

    ns : constant Series := Norm(s);

  begin
    Div(s,ns);
  end Normalize;

  function Normalize ( s : Series ) return Series is

    ns : constant Series := Norm(s);
    res : constant Series := s/ns;

  begin
    return res;
  end Normalize;

  function Max_Norm ( s : Series ) return quad_double is

    res : quad_double := QuadDobl_Complex_Numbers_Polar.Radius(s.cff(0));
    rad : quad_double;

  begin
    for i in 1..s.deg loop
      rad := QuadDobl_Complex_Numbers_Polar.Radius(s.cff(i));
      if rad > res
       then res := rad;
      end if;
    end loop;
    return res;
  end Max_Norm;

  function Two_Norm ( s : Series ) return quad_double is

    use QuadDobl_Complex_Numbers;

    c : constant Series := Conjugate(s);
    p : constant Series := c*s;
    res,cff : quad_double := create(0.0);

  begin
    for i in 0..p.deg loop
      cff := REAL_PART(p.cff(i));
      res := res + cff*cff;
    end loop;
    return QuadDobl_Mathematical_Functions.SQRT(res);
  end Two_Norm;

end QuadDobl_Dense_Series_Norms;
