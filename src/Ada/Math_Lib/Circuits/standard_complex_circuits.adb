with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Complex_Numbers;

package body Standard_Complex_Circuits is

  procedure Forward ( x : in Standard_Complex_Vectors.Link_to_Vector;
                      f : in Standard_Complex_Vectors.Link_to_Vector ) is

    use Standard_Complex_Numbers;

  begin
    f(f'first) := x(x'first)*x(x'first+1);
    for k in 2..x'last-1 loop
      f(k) := f(k-1)*x(k+1);
    end loop;
  end Forward;

  procedure Forward_Backward
              ( x,f,b : in Standard_Complex_Vectors.Link_to_Vector ) is

    use Standard_Complex_Numbers;

  begin
    f(f'first) := x(x'first)*x(x'first+1);
    for k in 2..x'last-1 loop
      f(k) := f(k-1)*x(k+1);
    end loop;
    b(b'first) := x(x'last)*x(x'last-1);
    for k in 2..x'last-2 loop
      b(k) := b(k-1)*x(x'last-k);
    end loop;
  end Forward_Backward;

  procedure Fused_Forward_Backward
              ( x,f,b : in Standard_Complex_Vectors.Link_to_Vector ) is

    use Standard_Complex_Numbers;

  begin
    f(f'first) := x(x'first)*x(x'first+1);
    b(b'first) := x(x'last)*x(x'last-1);
    for k in 2..x'last-2 loop
      f(k) := f(k-1)*x(k+1);
      b(k) := b(k-1)*x(x'last-k);
    end loop;
    if f'last > 1
     then f(f'last) := f(f'last-1)*x(x'last);
    end if;
  end Fused_Forward_Backward;

  procedure Forward_Backward_Cross
              ( x,f,b,c : in Standard_Complex_Vectors.Link_to_Vector ) is

    use Standard_Complex_Numbers;

  begin
    f(f'first) := x(x'first)*x(x'first+1);
    for k in 2..x'last-1 loop
      f(k) := f(k-1)*x(k+1);
    end loop;
    if x'last > 2 then
      b(b'first) := x(x'last)*x(x'last-1);
      for k in 2..x'last-2 loop
        b(k) := b(k-1)*x(x'last-k);
      end loop;
      if x'last = 3 then
        c(1) := x(1)*x(3);
      else
        c(1) := x(1)*b(x'last-3);
        for k in 2..x'last-3 loop
          c(k) := f(k-1)*b(x'last-2-k);
        end loop;
        c(x'last-2) := x(x'last)*f(x'last-3);
      end if;
    end if;
  end Forward_Backward_Cross;

  procedure Forward_Backward_Cross
              ( idx : in Standard_Integer_Vectors.Vector;
                x,f,b,c : in Standard_Complex_Vectors.Link_to_Vector ) is

    use Standard_Complex_Numbers;

  begin
    f(1) := x(idx(1))*x(idx(2));
    for k in 2..idx'last-1 loop
      f(k) := f(k-1)*x(idx(k+1));
    end loop;
    if idx'last > 2 then
      b(b'first) := x(idx(idx'last))*x(idx(idx'last-1));
      for k in 2..idx'last-2 loop
        b(k) := b(k-1)*x(idx(idx'last-k));
      end loop;
      if idx'last = 3 then
        c(1) := x(idx(1))*x(idx(3));
      else
        c(1) := x(idx(1))*b(idx'last-3);
        for k in 2..idx'last-3 loop
          c(k) := f(k-1)*b(idx'last-2-k);
        end loop;
        c(idx'last-2) := x(idx(idx'last))*f(idx'last-3);
      end if;
    end if;
  end Forward_Backward_Cross;

  procedure Fused_Forward_Backward_Cross
              ( x,f,b,c : in Standard_Complex_Vectors.Link_to_Vector ) is

    use Standard_Complex_Numbers;

    firstend,lastend,plusidx,minidx : integer32;

  begin
    if x'last >= 8 then
      if x'last mod 2 = 0 then
        lastend := x'last-4;
        firstend := lastend/2;
        f(f'first) := x(x'first)*x(x'first+1);
        b(b'first) := x(x'last)*x(x'last-1);
        for k in 2..firstend loop
          f(k) := f(k-1)*x(k+1);
          b(k) := b(k-1)*x(x'last-k);
        end loop;
        minidx := firstend+1; plusidx := minidx+1;
        for k in firstend+1..lastend loop
          f(k) := f(k-1)*x(k+1);
          b(k) := b(k-1)*x(x'last-k);
          c(plusidx) := f(plusidx-1)*b(x'last-2-plusidx);
          plusidx := plusidx + 1;
          c(minidx) := f(minidx-1)*b(x'last-2-minidx);
          minidx := minidx - 1;
        end loop;
      else
        lastend := x'last-4;
        firstend := (x'last-3)/2;
        f(f'first) := x(x'first)*x(x'first+1);
        b(b'first) := x(x'last)*x(x'last-1);
        for k in 2..firstend loop
          f(k) := f(k-1)*x(k+1);
          b(k) := b(k-1)*x(x'last-k);
        end loop;
        plusidx := firstend+1;
        c(plusidx) := f(plusidx-1)*b(x'last-2-plusidx);
        minidx := plusidx;
        for k in firstend+1..lastend loop
          f(k) := f(k-1)*x(k+1);
          b(k) := b(k-1)*x(x'last-k);
          plusidx := plusidx + 1;
          c(plusidx) := f(plusidx-1)*b(x'last-2-plusidx);
          minidx := minidx - 1;
          c(minidx) := f(minidx-1)*b(x'last-2-minidx);
        end loop;
      end if;
      plusidx := lastend+1;
      f(plusidx) := f(plusidx-1)*x(plusidx+1);
      b(plusidx) := b(plusidx-1)*x(x'last-plusidx);
      plusidx := plusidx+1;
      f(plusidx) := f(plusidx-1)*x(plusidx+1);
      b(plusidx) := b(plusidx-1)*x(x'last-plusidx);
      f(f'last) := f(f'last-1)*x(x'last);
      c(1) := x(1)*b(x'last-3);
      c(x'last-2) := x(x'last)*f(x'last-3);
    else
      f(f'first) := x(x'first)*x(x'first+1);
      b(b'first) := x(x'last)*x(x'last-1);
      for k in 2..x'last-2 loop
        f(k) := f(k-1)*x(k+1);
        b(k) := b(k-1)*x(x'last-k);
      end loop;
      if f'last > 1
       then f(f'last) := f(f'last-1)*x(x'last);
      end if;
      if x'last > 3 then
        c(1) := x(1)*b(x'last-3);
        for k in 2..x'last-3 loop
          c(k) := f(k-1)*b(x'last-2-k);
        end loop;
        c(x'last-2) := x(x'last)*f(x'last-3);
      elsif x'last = 3 then
        c(1) := x(1)*x(3);
      end if;
    end if;
  end Fused_Forward_Backward_Cross;

  function Allocate
             ( mxe : Standard_Integer_Vectors.Vector )
             return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(mxe'range);
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);

  begin
    for k in mxe'range loop
      if mxe(k) > 1 then
        res(k) := new Standard_Complex_Vectors.Vector'(1..mxe(k)-1 => zero);
      end if;
    end loop;
    return res;
  end Allocate;

  procedure Power_Table
              ( mxe : in Standard_Integer_Vectors.Vector;
                x : in Standard_Complex_Vectors.Link_to_Vector;
                pwt : in Standard_Complex_VecVecs.VecVec ) is

    lnk : Standard_Complex_Vectors.Link_to_Vector;

    use Standard_Complex_Numbers;

  begin
    for k in x'range loop
      if mxe(k) > 1 then
        lnk := pwt(k);
        lnk(1) := x(k)*x(k);
        for i in 2..mxe(k)-1 loop
          lnk(i) := lnk(i-1)*x(k);
        end loop;
      end if;
    end loop;
  end Power_Table;

end Standard_Complex_Circuits;
