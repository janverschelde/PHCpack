with unchecked_deallocation;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

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

-- EVALUATORS :

  function Eval ( m : Monomial; x : Vectors.Vector ) return Ring.number is

    res : Ring.number := m.cff;

    use Ring;

  begin
    for i in 1..integer32(m.nvr) loop
      for j in 1..m.exp(i) loop
        res := res*x(integer32(m.pos(i)));
      end loop;
    end loop;
    return res;
  end Eval;

  function Eval ( m : Link_to_Monomial;
                  x : Vectors.Vector ) return Ring.number is

    res : Ring.number := m.cff;

    use Ring;

  begin
    for i in 1..integer32(m.nvr) loop
      for j in 1..m.exp(i) loop
        res := res*x(integer32(m.pos(i)));
      end loop;
    end loop;
    return res;
  end Eval;

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
