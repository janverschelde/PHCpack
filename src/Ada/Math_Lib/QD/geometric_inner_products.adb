package body Geometric_Inner_Products is

  function Inner_Product ( dim : integer32; rtx,rty : double_float ) 
                         return double_float is

    result : double_float := 0.0;
    x : double_float := 1.0;
    y : double_float := 1.0;

  begin
    for i in 0..(dim-1) loop
      result := result + x*y;
      x := x*rtx;
      y := y*rty;
    end loop;
    return result;
  end Inner_Product;

  function Inner_Product ( dim : integer32; rtx,rty : double_double ) 
                         return double_double is

    result : double_double := create(0.0);
    x : double_double := create(1.0);
    y : double_double := create(1.0);

  begin
    for i in 0..(dim-1) loop
      result := result + x*y;
      x := x*rtx;
      y := y*rty;
    end loop;
    return result;
  end Inner_Product;

  function Inner_Product ( dim : integer32; rtx,rty : quad_double ) 
                         return quad_double is

    result : quad_double := create(0.0);
    x : quad_double := create(1.0);
    y : quad_double := create(1.0);

  begin
    for i in 0..(dim-1) loop
      result := result + x*y;
      x := x*rtx;
      y := y*rty;
    end loop;
    return result;
  end Inner_Product;

  function Inner_Product ( dim : integer32; rtx,rty : octo_double ) 
                         return octo_double is

    result : octo_double := create(0.0);
    x : octo_double := create(1.0);
    y : octo_double := create(1.0);

  begin
    for i in 0..(dim-1) loop
      result := result + x*y;
      x := x*rtx;
      y := y*rty;
    end loop;
    return result;
  end Inner_Product;

  function Inner_Product ( dim : integer32; rtx,rty : hexa_double ) 
                         return hexa_double is

    result : hexa_double := create(0.0);
    x : hexa_double := create(1.0);
    y : hexa_double := create(1.0);

  begin
    for i in 0..(dim-1) loop
      result := result + x*y;
      x := x*rtx;
      y := y*rty;
    end loop;
    return result;
  end Inner_Product;

end Geometric_Inner_Products;
