package body Geometric_Sums is

  function Double_Sum
             ( dim : integer32; ratio : double_float ) return double_float is

    gsum : double_float := 1.0;
    accu : double_float := ratio;

  begin
    for i in 1..dim loop
      gsum := gsum + accu;
      accu := ratio*accu;
    end loop;
    return gsum;
  end Double_Sum;

  function Double_Double_Sum
             ( dim : integer32; ratio : double_double ) return double_double is

    gsum : double_double := create(1.0);
    accu : double_double := ratio;

  begin
    for i in 1..dim loop
      gsum := gsum + accu;
      accu := ratio*accu;
    end loop;
    return gsum;
  end Double_Double_Sum;

  function Quad_Double_Sum
             ( dim : integer32; ratio : quad_double ) return quad_double is

    gsum : quad_double := create(1.0);
    accu : quad_double := ratio;

  begin
    for i in 1..dim loop
      gsum := gsum + accu;
      accu := ratio*accu;
    end loop;
    return gsum;
  end Quad_Double_Sum;

  function Octo_Double_Sum
             ( dim : integer32; ratio : octo_double ) return octo_double is

    gsum : octo_double := create(1.0);
    accu : octo_double := ratio;

  begin
    for i in 1..dim loop
      gsum := gsum + accu;
      accu := ratio*accu;
    end loop;
    return gsum;
  end Octo_Double_Sum;

  function Hexa_Double_Sum
             ( dim : integer32; ratio : hexa_double ) return hexa_double is

    gsum : hexa_double := create(1.0);
    accu : hexa_double := ratio;

  begin
    for i in 1..dim loop
      gsum := gsum + accu;
      accu := ratio*accu;
    end loop;
    return gsum;
  end Hexa_Double_Sum;

end Geometric_Sums;
