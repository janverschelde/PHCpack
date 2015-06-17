with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Multprec_Complex_Numbers;            use Multprec_Complex_Numbers;
with Standard_Speelpenning_Products;

package body Multprec_Speelpenning_Products is

  function Straight_Speel
             ( x : Multprec_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector is

    n : constant integer32 := x'last;
    res : Multprec_Complex_Vectors.Vector(0..n);

  begin
    Copy(x(1),res(0));
    Copy(x(1),res(n));
    for i in 2..(n-1) loop   -- does 2*(n-1) multiplications
      Mul(res(0),x(i));
      Mul(res(n),x(i));
    end loop;
    Mul(res(0),x(n));        -- 2*(n-1) + 1 = 2*n - 1
    for i in 1..(n-1) loop
      Copy(x(n),res(i));
      for j in 1..n-1 loop   -- does n*(n-2) multiplications
        if i /= j
         then Mul(res(i),x(j));
        end if;
      end loop;
    end loop;
    return res;
  end Straight_Speel;

  function Straight_Speel
             ( e : Standard_Natural_Vectors.Vector;
               x : Multprec_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector is

    n : constant integer32 := x'last;
    res : Multprec_Complex_Vectors.Vector(0..n);

  begin
    if e(1) = 0
     then res(0) := Create(integer(1));
     else Copy(x(1),res(0));
    end if;
    for i in 2..(n-1) loop   -- does |e|-1 multiplications
      if e(i) > 0
       then Mul(res(0),x(i));
      end if;
    end loop;
    if e(n) = 0 then
      res(n) := Create(integer(0));
    else
      res(0) := res(0)*x(n);
      if e(1) = 0
       then res(n) := Create(integer(1));
       else Copy(x(1),res(n));
      end if;
      for i in 2..(n-1) loop
        if e(i) > 0
         then Mul(res(n),x(i));
        end if;
      end loop;
    end if;
    for i in 1..(n-1) loop
      if e(i) = 0 then
        res(i) := Create(integer(0));
      else
        if e(n) = 0
         then res(i) := Create(integer(1));
         else Copy(x(n),res(i));
        end if;
        for j in 1..n-1 loop   -- does n*(n-2) multiplications
          if i /= j then
            if e(j) > 0
             then Mul(res(i),x(j));
            end if;
          end if;
        end loop;
      end if;
    end loop;
    return res;
  end Straight_Speel;

  function Reverse_Speel
             ( x : Multprec_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector is

    n : constant integer32 := x'last;
    res : Multprec_Complex_Vectors.Vector(0..n);
    fwd : Multprec_Complex_Vectors.Vector(1..n);
    back : Complex_Number;

  begin
    Copy(x(1),fwd(1));           -- fwd (forward) accumulates products
    for i in 2..n loop                     -- of consecutive variables
      fwd(i) := fwd(i-1)*x(i);             -- does n-1 multiplications
    end loop;
    Copy(fwd(n),res(0));
    Copy(fwd(n-1),res(n));
    Copy(x(n),back);              -- back accumulates backward products
    for i in reverse 2..n-1 loop  -- loop does 2*(n-1) multiplications
      res(i) := fwd(i-1)*back;
      Mul(back,x(i));
    end loop;
    Copy(back,res(1));
    Multprec_Complex_Vectors.Clear(fwd);
    Multprec_Complex_Numbers.Clear(back);
    return res;
  end Reverse_Speel;

  procedure Nonzeroes
             ( e : in Standard_Natural_Vectors.Vector;
               x : in Multprec_Complex_Vectors.Vector;
               idx : out Standard_Integer_Vectors.Vector;
               enz : out Standard_Natural_Vectors.Vector;
               xnz : out Multprec_Complex_Vectors.Vector ) is

    ind : integer32 := e'first-1;

  begin
    for i in e'range loop
      if e(i) /= 0 then
        ind := ind + 1;
        idx(ind) := i;
        enz(ind) := e(i);
        Copy(x(i),xnz(ind));
      end if;
    end loop;
  end Nonzeroes;

  function Indexed_Speel
             ( nx,nz : integer32;
               indnz : Standard_Integer_Vectors.Vector;
               x : Multprec_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(0..nx);
    nzx : Multprec_Complex_Vectors.Vector(1..nz);
    eva : Multprec_Complex_Vectors.Vector(0..nz);

  begin
    for k in 1..nz loop
      Copy(x(indnz(k)),nzx(k));
    end loop;
    eva := Reverse_Speel(nzx);
    Copy(eva(0),res(0));
    res(1..nx) := (1..nx => Create(integer(0)));
    for i in 1..nz loop
      Copy(eva(i),res(indnz(i)));
    end loop;
    Multprec_Complex_Vectors.Clear(nzx);
    Multprec_Complex_Vectors.Clear(eva);
    return res;
  end Indexed_Speel;

  function Reverse_Speel
             ( e : Standard_Natural_Vectors.Vector;
               x : Multprec_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector is

    n : constant integer32 := x'last;
    cntnz : integer32 := 0;         -- counts nonzero exponents in e
    res : Multprec_Complex_Vectors.Vector(0..n);
    fwd : Multprec_Complex_Vectors.Vector(1..n);
    idx : Standard_Integer_Vectors.Vector(1..n);
    ind : integer32 := 0;
    back : Complex_Number;

  begin
    for i in e'range loop             -- initialize forward products
      if e(i) /= 0
       then Copy(x(i),fwd(1)); ind := i+1; idx(1) := i; exit; 
      end if;
    end loop;
    res(1..n) := (1..n => Create(integer(0)));
    if ind = 0 then              -- case of a constant: all e(i) = 0
      res(0) := Create(integer(1));
    else
      cntnz := 1;
      for i in ind..e'last loop             -- make forward products
        if e(i) /= 0 then             -- skipping the zero exponents
          cntnz := cntnz + 1; idx(cntnz) := i;
          fwd(cntnz) := fwd(cntnz-1)*x(i);
        end if;
      end loop;
      if cntnz = 1 then
        Copy(x(ind-1),res(0));
        for i in 1..n loop
          if i = ind-1 
           then res(i) := Create(integer(1));
           else res(i) := Create(integer(0));
          end if;
        end loop;
      else
        Copy(fwd(cntnz),res(0));                -- value of monomial
        Copy(fwd(cntnz-1),res(idx(cntnz)));       -- last derivative
        Copy(x(idx(cntnz)),back);   -- accumulates backward products
        for i in reverse 2..cntnz-1 loop   
          res(idx(i)) := fwd(i-1)*back;
          Mul(back,x(idx(i)));
        end loop;
        Copy(back,res(idx(1)));
      end if;
    end if;
    Multprec_Complex_Vectors.Clear(fwd);
    Multprec_Complex_Numbers.Clear(back);
    return res;
  end Reverse_Speel;

  function Indexed_Reverse_Speel
             ( idx : Standard_Integer_Vectors.Vector;
               x : Multprec_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector is

    n : constant integer32 := x'last;
    res : Multprec_Complex_Vectors.Vector(0..n);
    fwd : Multprec_Complex_Vectors.Vector(idx'range);
    back : Complex_Number;

  begin
    res(1..n) := (1..n => Create(integer(0)));
    if idx'last < idx'first then               -- case of a constant
      res(0) := Create(integer(1));
    elsif idx'last = idx'first then          -- case of one variable
      Copy(x(idx(idx'first)),res(0));
      res(idx(idx'first)) := Create(integer(1));
    else
      Copy(x(idx(idx'first)),fwd(1));
      for i in idx'first+1..idx'last loop   -- make forward products
        fwd(i) := fwd(i-1)*x(idx(i));
      end loop;
      Copy(fwd(idx'last),res(0));               -- value of monomial
      Copy(fwd(idx'last-1),res(idx(idx'last)));   -- last derivative
      Copy(x(idx(idx'last)),back);  -- accumulates backward products
      for i in reverse 2..idx'last-1 loop   
        res(idx(i)) := fwd(i-1)*back;              -- cross products
        Mul(back,x(idx(i)));
      end loop;
      Copy(back,res(idx(1)));
    end if;
    Multprec_Complex_Vectors.Clear(fwd);
    Multprec_Complex_Numbers.Clear(back);
    return res;
  end Indexed_Reverse_Speel;

end Multprec_Speelpenning_Products;
