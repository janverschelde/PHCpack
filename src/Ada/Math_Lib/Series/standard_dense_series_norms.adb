with Standard_Complex_Numbers_Polar;
with Standard_Algebraic_Series;

package body Standard_Dense_Series_Norms is

  function Norm ( s : Series ) return Series is

    c : constant Series := Conjugate(s);
    p : constant Series := c*s;
    res : Series;

  begin
    res := Standard_Algebraic_Series.sqrt(p,0);
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

  function Max_Norm ( s : Series ) return double_float is

    res : double_float := Standard_Complex_Numbers_Polar.Radius(s.cff(0));
    rad : double_float;

  begin
    for i in 1..s.order loop
      rad := Standard_Complex_Numbers_Polar.Radius(s.cff(1));
      if rad > res
       then res := rad;
      end if;
    end loop;
    return res;
  end Max_Norm;

end Standard_Dense_Series_Norms;
