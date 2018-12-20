with unchecked_deallocation;

package body Generic_Monomials is

-- CONSTRUCTORS :

  function Create ( c : Ring.number; dim : natural32 ) return Monomial is

    res : Monomial;

  begin
    res.cff := c;
    res.dim := dim;
    res.nvr := 0;
    res.n_base := 0;
    return res;
  end Create;

  function Create ( c : Ring.number;
                    e : Standard_Natural_Vectors.Vector ) return Monomial is

    res : Monomial;
    invr,in_base,idxpos,idxbase : integer32;

  begin
    res.cff := c;
    res.dim := natural32(e'last);
    res.nvr := 0;
    res.n_base := 0;
    for i in e'range loop
      if e(i) > 0 then
        res.nvr := res.nvr + 1;
        if e(i) > 1
         then res.n_base := res.n_base + 1;
        end if;
      end if;
    end loop;
    if res.n_base > 0 then
      in_base := integer32(res.n_base);
      res.pos_base := new Standard_Natural_Vectors.Vector(1..in_base);
      res.exp_base := new Standard_Natural_Vectors.Vector(1..in_base);
      res.exp_tbl_base := new Standard_Natural_Vectors.Vector(1..in_base);
    end if;
    if res.nvr > 0 then
      invr := integer32(res.nvr);
      res.pos := new Standard_Natural_Vectors.Vector(1..invr);
      res.exp := new Standard_Natural_Vectors.Vector(1..invr);
      idxpos := 1; idxbase := 1;
      for i in e'range loop
        if e(i) > 0 then
          res.pos(idxpos) := natural32(i);
          res.exp(idxpos) := e(i);
          idxpos := idxpos + 1;
          if e(i) > 1 then
            res.pos_base(idxbase) := natural32(i);
            res.exp_base(idxbase) := e(i);
            res.exp_tbl_base(idxbase) := e(i)-2;
            idxbase := idxbase + 1;
          end if;
        end if;
      end loop;
    end if;
    return res;
  end Create;

  function Create ( c : Ring.number;
                    e : Standard_Natural_Vectors.Vector )
                  return Link_to_Monomial is

    res : Link_to_Monomial;
    mon : constant Monomial := Create(c,e);

  begin
    res := new Monomial'(mon);
    return res;
  end Create;

-- SELECTORS :

  function Degree ( m : Monomial ) return integer32 is

    res : integer32 := 0;

  begin
    if m.n_base = 0 then -- no degrees higher than one
      if m.nvr = 0 then
        if Ring.Equal(m.cff,Ring.zero)
         then res := -1;
         else res := 0;
        end if;
      else
        res := integer32(m.nvr);
      end if;
    else -- return the sum of the exponents of the variables in m
      for i in 1..integer32(m.nvr) loop
        res := res + integer32(m.exp(i));
      end loop;
    end if;
    return res;
  end Degree;

  function Degree ( m : Link_to_Monomial ) return integer32 is
  begin
    if m = null
     then return -1;
     else return Degree(m.all);
    end if;
  end Degree;

  function Largest_Exponent ( m : Monomial ) return natural32 is

    res : natural32 := 0;

  begin
    for i in 1..integer32(m.nvr) loop
      if m.exp(i) > res
       then res := m.exp(i);
      end if;
    end loop;
    return res;
  end Largest_Exponent;

  function Largest_Exponent ( m : Link_to_Monomial ) return natural32 is
  begin
    if m = null
     then return 0;
     else return Largest_Exponent(m.all);
    end if;
  end Largest_Exponent;

-- EVALUATION and DIFFERENTIATION :

  function Eval ( m : Monomial; x : Vectors.Vector ) return Ring.number is

    res : Ring.number;

    use Ring;

  begin
    Copy(m.cff,res);
    for i in 1..integer32(m.nvr) loop
      for j in 1..m.exp(i) loop
        Mul(res,x(integer32(m.pos(i)))); -- works in place
      end loop;
    end loop;
    return res;
  end Eval;

  function Eval ( m : Link_to_Monomial;
                  x : Vectors.Vector ) return Ring.number is

    res : Ring.number;

    use Ring;

  begin
    Copy(m.cff,res);
    for i in 1..integer32(m.nvr) loop
      for j in 1..m.exp(i) loop
        Mul(res,x(integer32(m.pos(i)))); -- works in place
      end loop;
    end loop;
    return res;
  end Eval;

  procedure Diff ( m : in Monomial; x : in Vectors.Vector;
                   yd : in out Vectors.Vector ) is

    idx : integer32;
    intfac : Ring.number;

    use Ring;

  begin
    for i in 1..integer32(m.nvr) loop
      idx := integer32(m.pos(i));
      Copy(m.cff,yd(i));
      if m.exp(i) > 1 then
        for k in 1..integer32(m.exp(i))-1 loop
          Mul(yd(i),x(idx));
        end loop;
        intfac := Ring.create(integer(m.exp(i)));
        Mul(yd(i),intfac);
      end if;
      for j in 1..integer32(m.nvr) loop
        if i /= j then
          idx := integer32(m.pos(j));
          for k in 1..integer32(m.exp(j)) loop
            Mul(yd(i),x(idx));
          end loop;
        end if;
      end loop;
    end loop;
  end Diff;

  procedure Diff ( m : in Link_to_Monomial; x : in Vectors.Vector;
                   yd : in out Vectors.Vector ) is
  begin
    if m /= null
     then Diff(m.all,x,yd);
    end if;
  end Diff;

  procedure Speel ( m : in Monomial; x : in Vectors.Vector;
                    y : in out Ring.number; yd : in out Vectors.Vector ) is
  begin
    if m.n_base = 0 then
      Speel_on_Product(m,x,y,yd);
    else
      declare
        powtab : VecVecs.VecVec(x'range);
      begin
        Power_Table(m,x,powtab);
        Speel(m,x,powtab,y,yd);
        VecVecs.Clear(powtab);
      end;
    end if;
  end Speel;

  procedure Speel ( m : in Link_to_Monomial; x : in Vectors.Vector;
                    y : in out Ring.number; yd : in out Vectors.Vector ) is
  begin
    if m.n_base = 0 then
      Speel_on_Product(m,x,y,yd);
    else
      declare
        powtab : VecVecs.VecVec(x'range);
      begin
        Power_Table(m,x,powtab);
        Speel(m,x,powtab,y,yd);
        VecVecs.Clear(powtab);
      end;
    end if;
  end Speel;

  procedure Speel_on_Product
              ( m : in Monomial; x : in Vectors.Vector;
                y : in out Ring.number; yd : in out Vectors.Vector ) is

    use Ring;

  begin
    Copy(x(integer32(m.pos(1))),yd(2));
    for i in 2..integer32(m.nvr-1) loop
      yd(i+1) := yd(i)*x(integer32(m.pos(i)));
    end loop;
    Copy(m.cff,y);
    for i in reverse 2..integer32(m.nvr) loop
      Mul(yd(i),y);
      Mul(y,x(integer32(m.pos(i))));
    end loop;
    Copy(y,yd(1));
    Mul(y,x(integer32(m.pos(1))));
  end Speel_on_Product;

  procedure Speel_on_Product
              ( m : in Link_to_Monomial; x : in Vectors.Vector;
                y : in out Ring.number; yd : in out Vectors.Vector ) is

    use Ring;

  begin
    Copy(x(integer32(m.pos(1))),yd(2));
    for i in 2..integer32(m.nvr-1) loop
      yd(i+1) := yd(i)*x(integer32(m.pos(i)));
    end loop;
    Copy(m.cff,y);
    for i in reverse 2..integer32(m.nvr) loop
      Mul(yd(i),y);
      Mul(y,x(integer32(m.pos(i))));
    end loop;
    Copy(y,yd(1));
    Mul(y,x(integer32(m.pos(1))));
  end Speel_on_Product;

  procedure Power_Table
              ( m : in monomial; x : Vectors.Vector;
                powtab : out VecVecs.VecVec ) is

    row,size : integer32;
    lv : Vectors.Link_to_Vector;

    use Ring;

  begin
    if m.n_base /= 0 then
      for i in 1..integer32(m.nvr) loop
        row := integer32(m.pos(i));
        size := integer32(m.exp(i));
        powtab(row) := new Vectors.Vector(0..size);
        lv := powtab(row);
        Ring.Copy(x(row),lv(0));
        for j in 1..size loop
          lv(j) := lv(j-1)*x(row);
        end loop;
      end loop;
    end if;
  end Power_Table;

  procedure Power_Table
              ( m : in Link_to_monomial; x : Vectors.Vector;
                powtab : out VecVecs.VecVec ) is
  begin
    if m /= null
     then Power_Table(m.all,x,powtab);
    end if;
  end Power_Table;

  procedure Common_Factor
              ( m : in Monomial; powtab : in VecVecs.VecVec;
                y : in out Ring.number ) is

    row,col : integer32;

  begin
    Ring.Copy(m.cff,y);
    for i in 1..integer32(m.n_base) loop
      row := integer32(m.pos_base(i));
      col := integer32(m.exp_tbl_base(i));
      Ring.Mul(y,powtab(row)(col));
    end loop;
  end Common_Factor;

  procedure Common_Factor
              ( m : in Monomial; powtab : in VecVecs.Link_to_VecVec;
                y : in out Ring.number ) is

    row,col : integer32;

  begin
    Ring.Copy(m.cff,y);
    for i in 1..integer32(m.n_base) loop
      row := integer32(m.pos_base(i));
      col := integer32(m.exp_tbl_base(i));
      Ring.Mul(y,powtab(row)(col));
    end loop;
  end Common_Factor;

  procedure Common_Factor
              ( m : in Link_to_Monomial; powtab : in VecVecs.VecVec;
                y : in out Ring.number ) is

    row,col : integer32;

  begin
    Ring.Copy(m.cff,y);
    for i in 1..integer32(m.n_base) loop
      row := integer32(m.pos_base(i));
      col := integer32(m.exp_tbl_base(i));
      Ring.Mul(y,powtab(row)(col));
    end loop;
  end Common_Factor;

  procedure Common_Factor
              ( m : in Link_to_Monomial;
                powtab : in VecVecs.Link_to_VecVec;
                y : in out Ring.number ) is

    row,col : integer32;

  begin
    Ring.Copy(m.cff,y);
    for i in 1..integer32(m.n_base) loop
      row := integer32(m.pos_base(i));
      col := integer32(m.exp_tbl_base(i));
      Ring.Mul(y,powtab(row)(col));
    end loop;
  end Common_Factor;

  procedure Speel ( m : in Monomial;
                    x : in Vectors.Vector; b : in Ring.number;
                    y : in out Ring.number; yd : in out Vectors.Vector ) is

    intfac : Ring.number;

    use Ring;

  begin
    if m.nvr > 1
     then Copy(x(integer32(m.pos(1))),yd(2));
    end if;
    for i in 2..integer32(m.nvr-1) loop
      yd(i+1) := yd(i)*x(integer32(m.pos(i)));
    end loop;
    Copy(b,y);
    for i in reverse 2..integer32(m.nvr) loop
      if m.exp(i) > 1 then
        intfac := Ring.create(integer(m.exp(i)));
        Mul(yd(i),intfac);
        Ring.clear(intfac);
      end if;
      Mul(yd(i),y);
      Mul(y,x(integer32(m.pos(i))));
    end loop;
    Copy(y,yd(1));
    if m.exp(1) > 1 then
      intfac := Ring.create(integer(m.exp(1)));
      Mul(yd(1),intfac);
      Ring.clear(intfac);
    end if;
    Mul(y,x(integer32(m.pos(1))));
  end Speel;

  procedure Speel ( m : in Link_to_Monomial;
                    x : in Vectors.Vector; b : in Ring.number;
                    y : in out Ring.number; yd : in out Vectors.Vector ) is

    intfac : Ring.number;

    use Ring;

  begin
    if m.nvr > 1
     then Copy(x(integer32(m.pos(1))),yd(2));
    end if;
    for i in 2..integer32(m.nvr-1) loop
      yd(i+1) := yd(i)*x(integer32(m.pos(i)));
    end loop;
    Copy(b,y);
    for i in reverse 2..integer32(m.nvr) loop
      if m.exp(i) > 1 then
        intfac := Ring.create(integer(m.exp(i)));
        Mul(yd(i),intfac);
        Ring.clear(intfac);
      end if;
      Mul(yd(i),y);
      Mul(y,x(integer32(m.pos(i))));
    end loop;
    Copy(y,yd(1));
    if m.exp(1) > 1 then
      intfac := Ring.create(integer(m.exp(1)));
      Mul(yd(1),intfac);
      Ring.clear(intfac);
    end if;
    Mul(y,x(integer32(m.pos(1))));
  end Speel;

  procedure Speel ( m : in Monomial;
                    x : in Vectors.Vector; powtab : in VecVecs.VecVec;
                    y : in out Ring.number; yd : in out Vectors.Vector ) is

    b : Ring.number;

  begin
    Common_Factor(m,powtab,b);
    Speel(m,x,b,y,yd);
  end Speel;

  procedure Speel ( m : in Monomial; x : in Vectors.Vector;
                    powtab : in VecVecs.Link_to_VecVec;
                    y : in out Ring.number; yd : in out Vectors.Vector ) is

    b : Ring.number;

  begin
    Common_Factor(m,powtab,b);
    Speel(m,x,b,y,yd);
  end Speel;

  procedure Speel ( m : in Link_to_Monomial; x : in Vectors.Vector;
                    powtab : in VecVecs.VecVec;
                    y : in out Ring.number; yd : in out Vectors.Vector ) is

    b : Ring.number;

  begin
    Common_Factor(m,powtab,b);
    Speel(m,x,b,y,yd);
  end Speel;

  procedure Speel ( m : in Link_to_Monomial; x : in Vectors.Vector;
                    powtab : in VecVecs.Link_to_VecVec;
                    y : in out Ring.number; yd : in out Vectors.Vector ) is

    b : Ring.number;

  begin
    Common_Factor(m,powtab,b);
    Speel(m,x,b,y,yd);
  end Speel;

-- DESTRUCTORS :
    
  procedure Clear ( m : in out Monomial ) is
  begin
    if m.nvr > 0 then
      Standard_Natural_Vectors.Clear(m.pos);
      Standard_Natural_Vectors.Clear(m.exp);
    end if;
    if m.n_base > 0 then
      Standard_Natural_Vectors.Clear(m.pos_base);
      Standard_Natural_Vectors.Clear(m.exp_base);
      Standard_Natural_Vectors.Clear(m.exp_tbl_base);
    end if;
  end Clear;

  procedure Clear ( m : in out Link_to_Monomial ) is

    procedure free is new unchecked_deallocation(Monomial,Link_to_Monomial);

  begin
    if m /= null then
      Clear(m.all);
      free(m);
    end if;
  end Clear;

end Generic_Monomials;
