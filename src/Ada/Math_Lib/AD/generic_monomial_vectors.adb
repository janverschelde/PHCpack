with unchecked_deallocation;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;

package body Generic_Monomial_Vectors is

-- CONSTRUCTORS :

  procedure Power_Update ( p : in out Polynomial ) is
  begin
    p.deg1 := Compute_Deg1(p.mons);
    Largest_Exponents(p.mons,p.maxexp);
    if not p.deg1
     then p.powtab := Allocate_Power_Table(p.maxexp);
    end if;
  end Power_Update;

  procedure Power_Update ( p : in out Link_to_Polynomial ) is
  begin
    if p /= null then
      p.deg1 := Compute_Deg1(p.mons);
      Largest_Exponents(p.mons,p.maxexp);
      if not p.deg1
       then p.powtab := Allocate_Power_Table(p.maxexp);
      end if;
    end if;
  end Power_Update;

  function Compute_Deg1 ( v : Monomial_Vector ) return boolean is
  begin
    for i in v'range loop
      if v(i).n_base > 0
       then return false;
      end if;
    end loop;
    return true;
  end Compute_Deg1;

  function Compute_Deg1 ( v : Link_to_Monomial_Vector ) return boolean is
  begin
    if v = null
     then return true;
     else return Compute_Deg1(v.all);
    end if;
  end Compute_Deg1;

  function Compute_Deg1 ( p : Polynomial ) return boolean is
  begin
    return Compute_Deg1(p.mons);
  end Compute_Deg1;

  function Compute_Deg1 ( p : Link_to_Polynomial ) return boolean is
  begin
    if p = null
     then return true;
     else return Compute_Deg1(p.mons);
    end if;
  end Compute_Deg1;

  procedure Largest_Exponents
              ( v : in Monomial_Vector;
                e : out Standard_Natural_Vectors.Vector ) is

    idx : integer32;

  begin
    e := (e'range => 0);
    for i in v'range loop
      for j in v(i).exp'range loop
        idx := integer32(v(i).pos(j));
        if v(i).exp(j) > e(idx)
         then e(idx) := v(i).exp(j);
        end if;
      end loop;
    end loop;
  end Largest_Exponents;

  procedure Largest_Exponents
              ( v : in Link_to_Monomial_Vector;
                e : out Standard_Natural_Vectors.Vector ) is
  begin
    if v /= null
     then Largest_Exponents(v.all,e);
    end if;
  end Largest_Exponents;

  procedure Largest_Exponents
              ( p : in Polynomial;
                e : out Standard_Natural_Vectors.Vector ) is
  begin
    Largest_Exponents(p.mons,e);
  end Largest_Exponents;

  procedure Largest_Exponents
              ( p : in Link_to_Polynomial;
                e : out Standard_Natural_Vectors.Vector ) is
  begin
    if p /= null
     then Largest_Exponents(p.all,e);
    end if;
  end Largest_Exponents;

  function Allocate_Power_Table
              ( maxexp : Standard_Natural_Vectors.Vector )
              return VecVecs.Link_to_VecVec is

    res : VecVecs.Link_to_VecVec;
    powtab : VecVecs.VecVec(maxexp'range);

  begin
    Allocate_Power_Table(powtab,maxexp);
    res := new VecVecs.VecVec'(powtab);
    return res;
  end Allocate_Power_Table;

  procedure Allocate_Power_Table
              ( powtab : in out VecVecs.VecVec;
                maxexp : in Standard_Natural_Vectors.Vector ) is
  begin
    for i in powtab'range loop
      powtab(i) := new Vectors.Vector(0..integer32(maxexp(i)));
    end loop;
  end Allocate_Power_Table;

  procedure Update_Power_Table
              ( powtab : in out VecVecs.VecVec; x : in Vectors.Vector ) is

    size : integer32;
    lv : Vectors.Link_to_Vector;

    use Ring;

  begin
    for i in powtab'range loop
      lv := powtab(i);
      Ring.Copy(x(i),lv(0));
      size := lv'last;
      for j in 1..size loop
        lv(j) := lv(j-1)*x(i);
      end loop;
    end loop;
  end Update_Power_Table;

  procedure Update_Power_Table
              ( powtab : in VecVecs.Link_to_VecVec;
                x : in Vectors.Vector ) is

    size : integer32;
    lv : Vectors.Link_to_Vector;

    use Ring;

  begin
    for i in powtab'range loop
      lv := powtab(i);
      Ring.Copy(x(i),lv(0));
      size := lv'last;
      for j in 1..size loop
        lv(j) := lv(j-1)*x(i);
      end loop;
    end loop;
  end Update_Power_Table;

-- SELECTORS :

  function Degree ( v : Monomial_Vector ) return integer32 is

    res : integer32 := Monomials.Degree(v(v'first));
    deg : integer32;

  begin
    for i in v'first+1..v'last loop
      deg := Monomials.Degree(v(i));
      if deg > res
       then res := deg;
      end if;
    end loop;
    return res;
  end Degree;

  function Degree ( v : Link_to_Monomial_Vector ) return integer32 is
  begin
    if v = null
     then return -1;
     else return Degree(v.all);
    end if;
  end Degree;

  function Degree ( p : Polynomial ) return integer32 is

    res : integer32 := Degree(p.mons);

  begin
    if res >= 0 then
      return res;
    elsif not Ring.Equal(p.cff0,Ring.zero) then
      return 0;
    else
      return -1;
    end if;
  end Degree;

  function Degree ( p : Link_to_Polynomial ) return integer32 is
  begin
    if p = null
     then return -1;
     else return Degree(p.all);
    end if;
  end Degree;

-- EVALUATION and DIFFERENTIATION :

  function Eval ( v : Monomial_Vector;
                  x : Vectors.Vector ) return Ring.number is

    res : Ring.number := Monomials.Eval(v(v'first),x);
    eva : Ring.number;

    use Ring;

  begin
    for i in v'first+1..v'last loop
      eva := Monomials.Eval(v(i),x);
      Add(res,eva);
      Clear(eva);
    end loop;
    return res;
  end Eval;

  function Eval ( v : Link_to_Monomial_Vector;
                  x : Vectors.Vector ) return Ring.number is

    res : Ring.number := Monomials.Eval(v(v'first),x);
    eva : Ring.number;

    use Ring;

  begin
    for i in v'first+1..v'last loop
      eva := Monomials.Eval(v(i),x);
      Add(res,eva);
      Clear(eva);
    end loop;
    return res;
  end Eval;

  function Eval ( p : Polynomial;
                  x : Vectors.Vector ) return Ring.number is

    res : Ring.number := Eval(p.mons,x);

  begin
    Ring.add(res,p.cff0);
    return res;
  end Eval;

  function Eval ( p : Link_to_Polynomial;
                  x : Vectors.Vector ) return Ring.number is

    res : Ring.number;

  begin
    if p = null
     then Ring.copy(Ring.zero,res);
     else res := Eval(p.all,x);
    end if;
    return res;
  end Eval;

  procedure Diff ( v : in Monomial_Vector; x : in Vectors.Vector;
                   yd,wrk : in out Vectors.Vector ) is
  begin
    for i in wrk'range loop
      Ring.copy(Ring.zero,yd(i));
    end loop;
    for i in v'range loop
      Monomials.Diff(v(i),x,wrk);
      for j in 1..integer32(v(i).nvr) loop
        Ring.add(yd(integer32(v(i).pos(j))),wrk(j));
      end loop;
    end loop;
  end Diff;

  procedure Diff ( v : in Monomial_Vector; x : in Vectors.Vector;
                   yd : in out Vectors.Vector ) is

    wrk : Vectors.Vector(x'range);

  begin
    Diff(v,x,yd,wrk);
  end Diff;

  procedure Diff ( v : in Link_to_Monomial_Vector; x : in Vectors.Vector;
                   yd,wrk : in out Vectors.Vector ) is
  begin
    if v /= null
     then Diff(v.all,x,yd,wrk);
    end if;
  end Diff;

  procedure Diff ( v : in Link_to_Monomial_Vector; x : in Vectors.Vector;
                   yd : in out Vectors.Vector ) is
  begin
    if v /= null
     then Diff(v.all,x,yd);
    end if;
  end Diff;

  procedure Speel ( p : in Polynomial; x : in Vectors.Vector;
                    y : in out Ring.number;
                    yd,wrk : in out Vectors.Vector ) is
  begin
    if p.deg1 then
      Speel_on_Product(p,x,y,yd,wrk);
    else
      Update_Power_Table(p.powtab,x);
      Speel(p.mons,x,p.powtab,y,yd,wrk);
      Ring.add(y,p.cff0);
    end if;
  end Speel;

  procedure Speel ( p : in Polynomial; x : in Vectors.Vector;
                    y : in out Ring.number; yd : in out Vectors.Vector ) is

    wrk : Vectors.Vector(x'range);

  begin
    Speel(p,x,y,yd,wrk);
  end Speel;

  procedure Speel ( p : in Link_to_Polynomial; x : in Vectors.Vector;
                    y : in out Ring.number;
                    yd,wrk : in out Vectors.Vector ) is
  begin
    if p /= null then
      if p.deg1 then
        Speel_on_Product(p,x,y,yd,wrk);
      else
        Update_Power_Table(p.powtab,x);
        Speel(p.mons,x,p.powtab,y,yd,wrk);
        Ring.add(y,p.cff0);
      end if;
    end if;
  end Speel;

  procedure Speel ( p : in Link_to_Polynomial; x : in Vectors.Vector;
                    y : in out Ring.number; yd : in out Vectors.Vector ) is

    wrk : Vectors.Vector(x'range);

  begin
    Speel(p,x,y,yd,wrk);
  end Speel;

  procedure Speel_on_Product
              ( v : in Monomial_Vector; x : in Vectors.Vector;
                y : in out Ring.number; yd,wrk : in out Vectors.Vector ) is

    mon : Monomials.Link_to_Monomial;
    val : Ring.number;

  begin
    for i in wrk'range loop
      Ring.copy(Ring.zero,yd(i));
    end loop;
    mon := v(v'first);
    Monomials.Speel_on_Product(mon,x,y,wrk);
    for j in 1..integer32(mon.nvr) loop
      Ring.add(yd(integer32(mon.pos(j))),wrk(j));
    end loop;
    for i in v'first+1..v'last loop
      mon := v(i);
      Monomials.Speel_on_Product(mon,x,val,wrk);
      Ring.add(y,val);
      for j in 1..integer32(mon.nvr) loop
        Ring.add(yd(integer32(mon.pos(j))),wrk(j));
      end loop;
    end loop;
  end Speel_on_Product;

  procedure Speel_on_Product
              ( v : in Monomial_Vector; x : in Vectors.Vector;
                y : in out Ring.number; yd : in out Vectors.Vector ) is

    wrk : Vectors.Vector(x'range);

  begin
    Speel_on_Product(v,x,y,yd,wrk);
  end Speel_on_Product;

  procedure Speel_on_Product
              ( v : in Link_to_Monomial_Vector; x : in Vectors.Vector;
                y : in out Ring.number; yd,wrk : in out Vectors.Vector ) is
  begin
    if v /= null
     then Speel_on_Product(v.all,x,y,yd,wrk);
    end if;
  end Speel_on_Product;

  procedure Speel_on_Product
              ( v : in Link_to_Monomial_Vector; x : in Vectors.Vector;
                y : in out Ring.number; yd : in out Vectors.Vector ) is
  begin
    if v /= null
     then Speel_on_Product(v.all,x,y,yd);
    end if;
  end Speel_on_Product;

  procedure Speel_on_Product
              ( p : in Polynomial; x : in Vectors.Vector;
                y : in out Ring.number; yd,wrk : in out Vectors.Vector ) is
  begin
    Speel_on_Product(p.mons,x,y,yd,wrk);
    Ring.add(y,p.cff0);
  end Speel_on_Product;

  procedure Speel_on_Product
              ( p : in Polynomial; x : in Vectors.Vector;
                y : in out Ring.number; yd : in out Vectors.Vector ) is
  begin
    Speel_on_Product(p.mons,x,y,yd);
    Ring.add(y,p.cff0);
  end Speel_on_Product;

  procedure Speel_on_Product
              ( p : in Link_to_Polynomial; x : in Vectors.Vector;
                y : in out Ring.number; yd,wrk : in out Vectors.Vector ) is
  begin
    if p /= null
     then Speel_on_Product(p.all,x,y,yd,wrk);
    end if;
  end Speel_on_Product;

  procedure Speel_on_Product
              ( p : in Link_to_Polynomial; x : in Vectors.Vector;
                y : in out Ring.number; yd : in out Vectors.Vector ) is
  begin
    if p /= null
     then Speel_on_Product(p.all,x,y,yd);
    end if;
  end Speel_on_Product;

  procedure Speel ( v : in Monomial_Vector;
                    x : in Vectors.Vector; powtab : in VecVecs.VecVec;
                    y : in out Ring.number;
                    yd,wrk : in out Vectors.Vector ) is

    mon : Monomials.Link_to_Monomial;
    val : Ring.number;

  begin
    for i in wrk'range loop
      Ring.copy(Ring.zero,yd(i));
    end loop;
    mon := v(v'first);
    if mon.n_base = 0
     then Monomials.Speel_on_Product(mon,x,y,wrk);
     else Monomials.Speel(mon,x,powtab,y,wrk);
    end if;
    for j in 1..integer32(mon.nvr) loop
      Ring.add(yd(integer32(mon.pos(j))),wrk(j));
    end loop;
    for i in v'first+1..v'last loop
      mon := v(i);
      if mon.n_base = 0
       then Monomials.Speel_on_Product(mon,x,val,wrk);
       else Monomials.Speel(mon,x,powtab,val,wrk);
      end if;
      Ring.add(y,val);
      for j in 1..integer32(mon.nvr) loop
        Ring.add(yd(integer32(mon.pos(j))),wrk(j));
      end loop;
    end loop;
  end Speel;

  procedure Speel ( v : in Monomial_Vector;
                    x : in Vectors.Vector;
                    powtab : in VecVecs.Link_to_VecVec;
                    y : in out Ring.number;
                    yd,wrk : in out Vectors.Vector ) is

    mon : Monomials.Link_to_Monomial;
    val : Ring.number;

  begin
    for i in wrk'range loop
      Ring.copy(Ring.zero,yd(i));
    end loop;
    mon := v(v'first);
    if mon.n_base = 0
     then Monomials.Speel_on_Product(mon,x,y,wrk);
     else Monomials.Speel(mon,x,powtab,y,wrk);
    end if;
    for j in 1..integer32(mon.nvr) loop
      Ring.add(yd(integer32(mon.pos(j))),wrk(j));
    end loop;
    for i in v'first+1..v'last loop
      mon := v(i);
      if mon.n_base = 0
       then Monomials.Speel_on_Product(mon,x,val,wrk);
       else Monomials.Speel(mon,x,powtab,val,wrk);
      end if;
      Ring.add(y,val);
      for j in 1..integer32(mon.nvr) loop
        Ring.add(yd(integer32(mon.pos(j))),wrk(j));
      end loop;
    end loop;
  end Speel;

  procedure Speel ( v : in Link_to_Monomial_Vector;
                    x : in Vectors.Vector; powtab : in VecVecs.VecVec;
                    y : in out Ring.number; yd,wrk : in out Vectors.Vector ) is
  begin
    if v /= null
     then Speel(v.all,x,powtab,y,yd,wrk);
    end if;
  end Speel;

  procedure Speel ( v : in Link_to_Monomial_Vector;
                    x : in Vectors.Vector; powtab : in VecVecs.Link_to_VecVec;
                    y : in out Ring.number; yd,wrk : in out Vectors.Vector ) is
  begin
    if v /= null
     then Speel(v.all,x,powtab,y,yd,wrk);
    end if;
  end Speel;

  procedure Speel_without_Cache
              ( v : in Monomial_Vector; x : in Vectors.Vector;
                y : in out Ring.number;
                yd,wrk : in out Vectors.Vector ) is

    mon : Monomials.Link_to_Monomial;
    val : Ring.number;

  begin
    for i in wrk'range loop
      Ring.copy(Ring.zero,yd(i));
    end loop;
    mon := v(v'first);
    Monomials.Speel(mon,x,y,wrk);
    for j in 1..integer32(mon.nvr) loop
      Ring.add(yd(integer32(mon.pos(j))),wrk(j));
    end loop;
    for i in v'first+1..v'last loop
      mon := v(i);
      Monomials.Speel(mon,x,val,wrk);
      Ring.add(y,val);
      for j in 1..integer32(mon.nvr) loop
        Ring.add(yd(integer32(mon.pos(j))),wrk(j));
      end loop;
    end loop;
  end Speel_without_Cache;

  procedure Speel_without_Cache
              ( p : in Polynomial; x : in Vectors.Vector;
                y : in out Ring.number; yd : in out Vectors.Vector ) is

    wrk : Vectors.Vector(x'range);

  begin
    Speel_without_Cache(p.mons,x,y,yd,wrk);
    Ring.add(y,p.cff0);
  end Speel_without_Cache;

  procedure Speel_without_Cache
              ( p : in Link_to_Polynomial; x : in Vectors.Vector;
                y : in out Ring.number; yd : in out Vectors.Vector ) is

    wrk : Vectors.Vector(x'range);

  begin
    if p /= null then
      Speel_without_Cache(p.mons,x,y,yd,wrk);
      Ring.add(y,p.cff0);
    end if;
  end Speel_without_Cache;

-- DESTRUCTORS :

  procedure Clear ( v : in out Monomial_Vector ) is
  begin
    for i in v'range loop
      Monomials.Clear(v(i));
    end loop;
  end Clear;

  procedure Clear ( v : in out Link_to_Monomial_Vector ) is

    procedure free is
      new unchecked_deallocation(Monomial_Vector,Link_to_Monomial_Vector);

  begin
    if v /= null then
      Clear(v.all);
      free(v);
    end if;
  end Clear;

  procedure Clear ( p : in out Polynomial ) is
  begin
    Clear(p.mons);
    if not p.deg1
     then VecVecs.Deep_Clear(p.powtab);
    end if;
  end Clear;

  procedure Clear ( p : in out Link_to_Polynomial ) is

    procedure free is
      new unchecked_deallocation(Polynomial,Link_to_Polynomial);

  begin
    if p /= null then
      Clear(p.all);
      free(p);
    end if;
  end Clear;

end Generic_Monomial_Vectors;
