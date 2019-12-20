package body Generic_Speelpenning_Convolutions is

  function Allocate_Coefficients
             ( deg : integer32 ) return Vectors.Link_to_Vector is

    cff : constant Vectors.Vector(0..deg) := (0..deg => Ring.zero);
    res : constant Vectors.Link_to_Vector := new Vectors.Vector'(cff);

  begin
    return res;
  end Allocate_Coefficients;

  function Allocate_Coefficients
             ( dim,deg : integer32 ) return VecVecs.VecVec is

    res : VecVecs.VecVec(1..dim);

  begin
    for k in 1..dim loop
      res(k) := Allocate_Coefficients(deg);
    end loop;
    return res;
  end Allocate_Coefficients;

  procedure Update ( values : in Vectors.Link_to_Vector;
                     inc : in Vectors.Link_to_Vector ) is

    use Ring;

  begin
    for i in values'range loop
      values(i) := values(i) + inc(i);
    end loop;
  end Update;

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
      if x'last = 3 then
        Multiply(x(1),x(3),cross(1));
      else
        Multiply(x(1),backward(x'last-3),cross(1));
        for k in 2..x'last-3 loop
          Multiply(forward(k-1),backward(x'last-2-k),cross(k));
        end loop;
        Multiply(forward(x'last-3),x(x'last),cross(x'last-2));
      end if;
    end if;
  end Speel;

  procedure Speel ( x : in VecVecs.VecVec;
                    idx : in Standard_Integer_Vectors.Vector;
                    forward,backward,cross : in out VecVecs.VecVec ) is
  begin
    Multiply(x(idx(1)),x(idx(2)),forward(1));
    for k in 3..idx'last loop
      Multiply(forward(k-2),x(idx(k)),forward(k-1));
    end loop;
    if idx'last > 2 then
      Multiply(x(idx(idx'last)),x(idx(idx'last-1)),backward(1));
      for k in 2..idx'last-2 loop
        Multiply(backward(k-1),x(idx(idx'last-k)),backward(k));
      end loop;
      if idx'last = 3 then
        Multiply(x(idx(1)),x(idx(3)),cross(1));
      else
        Multiply(x(idx(1)),backward(idx'last-3),cross(1));
        for k in 2..idx'last-3 loop
          Multiply(forward(k-1),backward(idx'last-2-k),cross(k));
        end loop;
        Multiply(forward(idx'last-3),x(idx(idx'last)),cross(idx'last-2));
      end if;
    end if;
  end Speel;

  procedure Speel ( idx : in Standard_Integer_VecVecs.VecVec;
                    x : in VecVecs.VecVec;
                    forward,backward,cross,yd : in out VecVecs.VecVec ) is

    use Standard_Integer_Vectors;
    use Ring;

    idk : Standard_Integer_Vectors.Link_to_Vector;
    yptr : constant Vectors.Link_to_Vector := yd(yd'last);

  begin
    for k in idx'range loop
      idk := idx(k);
      if idk /= null then
        if idk'last = 1 then
          Update(yptr,x(idk(1)));
          yd(idk(1))(0) := yd(idk(1))(0) + Ring.one;
        else
          Speel(x,idk.all,forward,backward,cross);
          Update(yptr,forward(idk'last-1));
          if idk'last = 2 then
            Update(yd(idk(2)),x(idk(1)));
            Update(yd(idk(1)),x(idk(2)));
          else -- idk'last > 2 
            Update(yd(idk(1)),backward(idk'last-2));
            for j in idk'first+1..idk'last-1 loop
              Update(yd(idk(j)),cross(j-1));
            end loop;
            Update(yd(idk(idk'last)),forward(idk'last-2));
          end if;
        end if;
      end if;
    end loop;
  end Speel;

  procedure Speel ( idx : in Standard_Integer_VecVecs.VecVec;
                    cff : in VecVecs.VecVec; x : in VecVecs.VecVec;
                    forward,backward,cross,yd : in out VecVecs.VecVec;
                    wrk : Vectors.Link_to_Vector ) is

    use Standard_Integer_Vectors;

    idk : Standard_Integer_Vectors.Link_to_Vector;
    yptr : constant Vectors.Link_to_Vector := yd(yd'last);
    pcff : Vectors.Link_to_Vector;

  begin
    for k in idx'range loop
      idk := idx(k);
      if idk /= null then
        pcff := cff(k);
        if idk'last = 1 then
          Multiply(pcff,x(idk(1)),wrk);
          Update(yptr,wrk);
          Update(yd(idk(1)),pcff);
        else
          Speel(x,idk.all,forward,backward,cross);
          Multiply(pcff,forward(idk'last-1),wrk);
          Update(yptr,wrk);
          if idk'last = 2 then
            Multiply(pcff,x(idk(1)),wrk);
            Update(yd(idk(2)),wrk);
            Multiply(pcff,x(idk(2)),wrk);
            Update(yd(idk(1)),wrk);
          else -- idk'last > 2 
            Multiply(pcff,backward(idk'last-2),wrk);
            Update(yd(idk(1)),wrk);
            for j in idk'first+1..idk'last-1 loop
              Multiply(pcff,cross(j-1),wrk);
              Update(yd(idk(j)),wrk);
            end loop;
            Multiply(pcff,forward(idk'last-2),wrk);
            Update(yd(idk(idk'last)),wrk);
          end if;
        end if;
      end if;
    end loop;
  end Speel;

end Generic_Speelpenning_Convolutions;
