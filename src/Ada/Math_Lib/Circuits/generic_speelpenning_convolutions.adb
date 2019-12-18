package body Generic_Speelpenning_Convolutions is

  function Allocate_Coefficients
             ( dim,deg : integer32 ) return VecVecs.VecVec is

    res : VecVecs.VecVec(1..dim);

  begin
    for k in 1..dim loop
      declare
        cff : constant Vectors.Vector(0..deg) := (0..deg => Ring.zero);
      begin
        res(k) := new Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end Allocate_Coefficients;

  procedure Multiply ( first,second,product : in Vectors.Link_to_Vector ) is

    deg : constant integer32 := first'last;

    use Ring;

  begin
    product(0) := first(0)*second(0);
    for k in 1..deg loop
      product(k) := first(0)*second(k);
      for i in 1..k loop
        product(k) := product(k) + first(i)*second(k-i);
      end loop;
    end loop;
  end Multiply;

  procedure Speel ( x : in VecVecs.VecVec;
                    forward,backward,cross : in out VecVecs.VecVec ) is
  begin
    Multiply(x(1),x(2),forward(1));
    for k in 3..x'last loop
      Multiply(forward(k-2),x(k),forward(k-1));
    end loop;
    if x'last > 2 then
      Multiply(x(x'last),x(x'last-1),backward(1));
      for k in 2..x'last-2 loop
        Multiply(backward(k-1),x(x'last-k),backward(k));
      end loop;
      Multiply(x(1),backward(x'last-3),cross(1));
      for k in 2..x'last-3 loop
        Multiply(forward(k-1),backward(x'last-2-k),cross(k));
      end loop;
      Multiply(forward(x'last-3),x(x'last),cross(x'last-2));
    end if;
  end Speel;

end Generic_Speelpenning_Convolutions;
