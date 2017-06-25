with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;

package body Permute_Operations is

  function "*" ( p : Permutation; v : Standard_Natural_Vectors.Vector )
	       return Standard_Natural_Vectors.Vector is

    r : Standard_Natural_Vectors.Vector(v'range);

  begin
    for i in p'range loop
      if p(i) >= 0
       then r(i) := v(p(i));
       else r(i) := natural32(-integer32(v(-p(i))));
      end if;
    end loop;
    return r;
  end "*";

  function "*" ( p : Permutation; v : Standard_Integer_Vectors.Vector )
	       return Standard_Integer_Vectors.Vector is

    r : Standard_Integer_Vectors.Vector(v'range);

  begin
    for i in p'range loop
      if p(i) >= 0
       then r(i) := v(p(i));
       else r(i) := -v(-p(i));
      end if;
    end loop;
    return r;
  end "*";

  function "*" ( p : Permutation; v : Standard_Floating_Vectors.Vector )
               return Standard_Floating_Vectors.Vector is

    r : Standard_Floating_Vectors.Vector(v'range);

  begin
    for i in p'range loop
      if p(i) >= 0
       then r(i) := v(p(i));
       else r(i) := -v(-p(i));
      end if;
    end loop;
    return r;
  end "*";

  function "*" ( p : Permutation; v : Standard_Complex_Vectors.Vector )
               return Standard_Complex_Vectors.Vector is

    r : Standard_Complex_Vectors.Vector(v'range);
    use Standard_Complex_Numbers;

  begin
    for i in p'range loop
      if p(i) >= 0
       then r(i) := v(p(i));
       else r(i) := -v(-p(i));
      end if;
    end loop;
    return r;
  end "*";

  function "*" ( p : Permutation; v : DoblDobl_Complex_Vectors.Vector )
               return DoblDobl_Complex_Vectors.Vector is

    r : DoblDobl_Complex_Vectors.Vector(v'range);
    use DoblDobl_Complex_Numbers;

  begin
    for i in p'range loop
      if p(i) >= 0
       then r(i) := v(p(i));
       else r(i) := -v(-p(i));
      end if;
    end loop;
    return r;
  end "*";

  function "*" ( p : Permutation; v : QuadDobl_Complex_Vectors.Vector )
               return QuadDobl_Complex_Vectors.Vector is

    r : QuadDobl_Complex_Vectors.Vector(v'range);
    use QuadDobl_Complex_Numbers;

  begin
    for i in p'range loop
      if p(i) >= 0
       then r(i) := v(p(i));
       else r(i) := -v(-p(i));
      end if;
    end loop;
    return r;
  end "*";

  function Permutable ( v1,v2 : Standard_Natural_Vectors.Vector )
                      return boolean is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last then
      return false;  -- the dimensions must correspond !
    else
      declare
        p : Permutation(v1'first..v1'last);
      begin
        for k in p'range loop
          p(k) := 0;
          for l in v2'range loop
            if v2(l) = v1(k) then
              p(k) := l;
              for j in 1..(k-1) loop
                if p(j) = l
                 then p(k) := 0;
                end if;
              end loop;
            end if;
            exit when p(k) /= 0;
          end loop;
          if p(k) = 0
           then return false;
          end if;
        end loop;
      end;
      return true;
    end if;
  end Permutable;

  function Permutable ( v1,v2 : Standard_Integer_Vectors.Vector )
                      return boolean is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last then
      return false;  -- the dimensions must correspond !
    else
      declare
        p : Permutation(v1'first..v1'last);
      begin
        for k in p'range loop
          p(k) := 0;
          for l in v2'range loop
            if v2(l) = v1(k) then
              p(k) := l;
              for j in 1..(k-1) loop
                if p(j) = l
                 then p(k) := 0;
                end if;
              end loop;
            end if;
            exit when p(k) /= 0;
          end loop;
          if p(k) = 0
           then return false;
          end if;
        end loop;
      end;
      return true;
    end if;
  end Permutable;

  function Permutable ( v1,v2 : Standard_Floating_Vectors.Vector )
                      return boolean is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last then
      return false;  -- the dimensions must correspond !
    else
      declare
        p : Permutation(v1'first..v1'last);
      begin
        for k in p'range loop
          p(k) := 0;
          for l in v2'range loop
            if v2(l) = v1(k) then
              p(k) := l;
              for j in 1..(k-1) loop
                if p(j) = l
                 then p(k) := 0;
                end if;
              end loop;
            end if;
            exit when p(k) /= 0;
          end loop;
          if p(k) = 0
           then return false;
          end if;
        end loop;
      end;
      return true;
    end if;
  end Permutable;

  function Permutable ( v1,v2 : Standard_Complex_Vectors.Vector )
                      return boolean is

    use Standard_Complex_Numbers;

  begin
    if v1'first /= v2'first or else v1'last /= v2'last then
      return false;  -- the dimensions must correspond !
    else
      declare
        p : Permutation(v1'first..v1'last);
      begin
        for k in p'range loop
          p(k) := 0;
          for l in v2'range loop
            if v2(l) = v1(k) then
              p(k) := l;
              for j in 1..(k-1) loop
                if p(j) = l
                 then p(k) := 0;
                end if;
              end loop;
            end if;
            exit when p(k) /= 0;
          end loop;
          if p(k) = 0
           then return false;
          end if;
        end loop;
      end;
      return true;
    end if;
  end Permutable;

  function Permutable ( v1,v2 : Standard_Floating_Vectors.Vector;
                        tol : double_float ) return boolean is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last then
      return false;  -- the dimensions must correspond !
    else
      declare
        p : Permutation(v1'first..v1'last);
      begin
        for k in p'range loop
          p(k) := 0;
          for l in v2'range loop
            if ABS(v2(l) - v1(k)) <= tol then
              p(k) := l;
              for j in 1..(k-1) loop
                if p(j) = l
                 then p(k) := 0;
                end if;
              end loop;
            end if;
            exit when p(k) /= 0;
          end loop;
          if p(k) = 0
           then return false;
          end if;
        end loop;
      end;
      return true;
    end if;
  end Permutable;

  function Permutable ( v1,v2 : Standard_Complex_Vectors.Vector;
                        tol : double_float ) return boolean is

    use Standard_Complex_Numbers;

  begin
    if v1'first /= v2'first or else v1'last /= v2'last then
      return false;  -- the dimensions must correspond !
    else
      declare
        p : Permutation(v1'first..v1'last);
      begin
        for k in p'range loop
          p(k) := 0;
          for l in v2'range loop
            if (ABS(REAL_PART(v2(l)) - REAL_PART(v1(k))) <= tol)
             and then (ABS(IMAG_PART(v2(l)) - IMAG_PART(v1(k))) <= tol) then
              p(k) := l;
              for j in 1..(k-1) loop
                if p(j) = l
                 then p(k) := 0;
                end if;
              end loop;
            end if;
            exit when p(k) /= 0;
          end loop;
          if p(k) = 0
           then return false;
          end if;
        end loop;
      end;
      return true;
    end if;
  end Permutable;

  function Sign_Permutable ( v1,v2 : Standard_Natural_Vectors.Vector )
                           return boolean is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last then
      return false;  -- the dimensions must correspond !
    else
      declare
        p : Permutation(v1'first..v1'last);
      begin
        for k in p'range loop
          p(k) := 0;
          for l in v2'range loop
            if v2(l) = v1(k) or else v2(l) = -v1(k) then
              p(k) := l;
              for j in 1..(k-1) loop
                if p(j) = l
                 then p(k) := 0;
                end if;
              end loop;
            end if;
            exit when p(k) /= 0;
          end loop;
          if p(k) = 0
           then return false;
          end if;
        end loop;
      end;
      return true;
    end if;
  end Sign_Permutable;

  function Sign_Permutable ( v1,v2 : Standard_Integer_Vectors.Vector )
                           return boolean is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last then
      return false;  -- the dimensions must correspond !
    else
      declare
        p : Permutation(v1'first..v1'last);
      begin
        for k in p'range loop
          p(k) := 0;
          for l in v2'range loop
            if v2(l) = v1(k) or else v2(l) = -v1(k) then
              p(k) := l;
              for j in 1..(k-1) loop
                if p(j) = l
                 then p(k) := 0;
                end if;
              end loop;
            end if;
            exit when p(k) /= 0;
          end loop;
          if p(k) = 0
           then return false;
          end if;
        end loop;
      end;
      return true;
    end if;
  end Sign_Permutable;

  function Sign_Permutable ( v1,v2 : Standard_Floating_Vectors.Vector )
                           return boolean is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last then
      return false;  -- the dimensions must correspond !
    else
      declare
        p : Permutation(v1'first..v1'last);
      begin
        for k in p'range loop
          p(k) := 0;
          for l in v2'range loop
            if v2(l) = v1(k) or else v2(l) = -v1(k) then
              p(k) := l;
              for j in 1..(k-1) loop
                if p(j) = l
                 then p(k) := 0;
                end if;
              end loop;
            end if;
            exit when p(k) /= 0;
          end loop;
          if p(k) = 0
           then return false;
          end if;
        end loop;
      end;
      return true;
    end if;
  end Sign_Permutable;

  function Sign_Permutable ( v1,v2 : Standard_Complex_Vectors.Vector )
                           return boolean is

    use Standard_Complex_Numbers;

  begin
    if v1'first /= v2'first or else v1'last /= v2'last then
      return false;  -- the dimensions must correspond !
    else
      declare
        p : Permutation(v1'first..v1'last);
      begin
        for k in p'range loop
          p(k) := 0;
          for l in v2'range loop
            if v2(l) = v1(k) or else v2(l) = -v1(k) then
              p(k) := l;
              for j in 1..(k-1) loop
                if p(j) = l
                 then p(k) := 0;
                end if;
              end loop;
            end if;
            exit when p(k) /= 0;
          end loop;
          if p(k) = 0
           then return false;
          end if;
        end loop;
      end;
      return true;
    end if;
  end Sign_Permutable;

  function Sign_Permutable ( v1,v2 : Standard_Floating_Vectors.Vector;
                             tol : double_float ) return boolean is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last then
      return false;  -- the dimensions must correspond !
    else
      declare
        p : Permutation(v1'first..v1'last);
      begin
        for k in p'range loop
          p(k) := 0;
          for l in v2'range loop
            if (ABS(v2(l) - v1(k)) <= tol)
       or else (ABS(v2(l) + v1(k)) <= tol)
            then
              p(k) := l;
              for j in 1..(k-1) loop
                if p(j) = l
                 then p(k) := 0;
                end if;
              end loop;
            end if;
            exit when p(k) /= 0;
          end loop;
          if p(k) = 0
           then return false;
          end if;
        end loop;
      end;
      return true;
    end if;
  end Sign_Permutable;

  function Sign_Permutable ( v1,v2 : Standard_Complex_Vectors.Vector;
                             tol : double_float ) return boolean is

    use Standard_Complex_Numbers;

  begin
    if v1'first /= v2'first or else v1'last /= v2'last then
      return false;  -- the dimensions must correspond !
    else
      declare
        p : Permutation(v1'first..v1'last);
      begin
        for k in p'range loop
          p(k) := 0;
          for l in v2'range loop
            if ((ABS(REAL_PART(v2(l)) - REAL_PART(v1(k))) <= tol)
       and then (ABS(IMAG_PART(v2(l)) - IMAG_PART(v1(k))) <= tol))
       or else ((ABS(REAL_PART(v2(l)) + REAL_PART(v1(k))) <= tol)
       and then (ABS(IMAG_PART(v2(l)) + IMAG_PART(v1(k))) <= tol))
            then
              p(k) := l;
              for j in 1..(k-1) loop
                if p(j) = l 
                 then p(k) := 0;
                end if;
              end loop;
            end if;
            exit when p(k) /= 0;
          end loop;
          if p(k) = 0
           then return false;
          end if;
        end loop;
      end;
      return true;
    end if;
  end Sign_Permutable;

  function "*" ( p : Permutation; t : Standard_Complex_Polynomials.Term )
               return Standard_Complex_Polynomials.Term is

    res : Standard_Complex_Polynomials.Term;
    use Standard_Complex_Numbers;

  begin
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector(t.dg'range);
    for i in p'range loop
      if p(i) >= 0 then
        res.dg(i) := t.dg(p(i));
      else
        res.dg(i) := t.dg(-p(i));
        res.cf := -res.cf;
      end if;
    end loop;
    return res;
  end "*";

  function "*" ( p : Permutation; t : DoblDobl_Complex_Polynomials.Term )
               return DoblDobl_Complex_Polynomials.Term is

    res : DoblDobl_Complex_Polynomials.Term;
    use DoblDobl_Complex_Numbers;

  begin
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector(t.dg'range);
    for i in p'range loop
      if p(i) >= 0 then
        res.dg(i) := t.dg(p(i));
      else
        res.dg(i) := t.dg(-p(i));
        res.cf := -res.cf;
      end if;
    end loop;
    return res;
  end "*";

  function "*" ( p : Permutation; t : QuadDobl_Complex_Polynomials.Term )
               return QuadDobl_Complex_Polynomials.Term is

    res : QuadDobl_Complex_Polynomials.Term;
    use QuadDobl_Complex_Numbers;

  begin
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector(t.dg'range);
    for i in p'range loop
      if p(i) >= 0 then
        res.dg(i) := t.dg(p(i));
      else
        res.dg(i) := t.dg(-p(i));
        res.cf := -res.cf;
      end if;
    end loop;
    return res;
  end "*";

  function "*" ( p : Permutation; s : Standard_Complex_Polynomials.Poly )
               return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;
    res : Poly := Null_Poly;

    procedure Permute_Term ( t : in Term; continue : out boolean ) is
      tt : Term := p*t;
    begin
      Add(res,tt);
      Clear(tt);
      continue := true;
    end Permute_Term;
    procedure Permute_Terms is new Visiting_Iterator(Permute_Term);

  begin
    Permute_Terms(s);
    return res;
  end "*";

  function "*" ( p : Permutation; s : DoblDobl_Complex_Polynomials.Poly )
               return DoblDobl_Complex_Polynomials.Poly is

    use DoblDobl_Complex_Polynomials;
    res : Poly := Null_Poly;

    procedure Permute_Term ( t : in Term; continue : out boolean ) is
      tt : Term := p*t;
    begin
      Add(res,tt);
      Clear(tt);
      continue := true;
    end Permute_Term;
    procedure Permute_Terms is new Visiting_Iterator(Permute_Term);

  begin
    Permute_Terms(s);
    return res;
  end "*";

  function "*" ( p : Permutation; s : QuadDobl_Complex_Polynomials.Poly )
               return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Polynomials;
    res : Poly := Null_Poly;

    procedure Permute_Term ( t : in Term; continue : out boolean ) is
      tt : Term := p*t;
    begin
      Add(res,tt);
      Clear(tt);
      continue := true;
    end Permute_Term;
    procedure Permute_Terms is new Visiting_Iterator(Permute_Term);

  begin
    Permute_Terms(s);
    return res;
  end "*";

  function "*" ( p : Permutation; t : Standard_Complex_Laurentials.Term )
               return Standard_Complex_Laurentials.Term is

    res : Standard_Complex_Laurentials.Term;
    use Standard_Complex_Numbers;

  begin
    res.cf := t.cf;
    res.dg := new Standard_Integer_Vectors.Vector(t.dg'range);
    for i in p'range loop
      if p(i) >= 0 then
        res.dg(i) := t.dg(p(i));
      else
        res.dg(i) := t.dg(-p(i));
        res.cf := -res.cf;
      end if;
    end loop;
    return res;
  end "*";

  function "*" ( p : Permutation; t : DoblDobl_Complex_Laurentials.Term )
               return DoblDobl_Complex_Laurentials.Term is

    res : DoblDobl_Complex_Laurentials.Term;
    use DoblDobl_Complex_Numbers;

  begin
    res.cf := t.cf;
    res.dg := new Standard_Integer_Vectors.Vector(t.dg'range);
    for i in p'range loop
      if p(i) >= 0 then
        res.dg(i) := t.dg(p(i));
      else
        res.dg(i) := t.dg(-p(i));
        res.cf := -res.cf;
      end if;
    end loop;
    return res;
  end "*";

  function "*" ( p : Permutation; t : QuadDobl_Complex_Laurentials.Term )
               return QuadDobl_Complex_Laurentials.Term is

    res : QuadDobl_Complex_Laurentials.Term;
    use QuadDobl_Complex_Numbers;

  begin
    res.cf := t.cf;
    res.dg := new Standard_Integer_Vectors.Vector(t.dg'range);
    for i in p'range loop
      if p(i) >= 0 then
        res.dg(i) := t.dg(p(i));
      else
        res.dg(i) := t.dg(-p(i));
        res.cf := -res.cf;
      end if;
    end loop;
    return res;
  end "*";

  function "*" ( p : Permutation; s : Standard_Complex_Laurentials.Poly )
               return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Laurentials;
    res : Poly := Null_Poly;

    procedure Permute_Term ( t : in Term; continue : out boolean ) is

      tt : Term := p*t;

    begin
      Add(res,tt);
      Clear(tt);
      continue := true;
    end Permute_Term;
    procedure Permute_Terms is new Visiting_Iterator(Permute_Term);

  begin
    Permute_Terms(s);
    return res;
  end "*";

  function "*" ( p : Permutation; s : DoblDobl_Complex_Laurentials.Poly )
               return DoblDobl_Complex_Laurentials.Poly is

    use DoblDobl_Complex_Laurentials;
    res : Poly := Null_Poly;

    procedure Permute_Term ( t : in Term; continue : out boolean ) is

      tt : Term := p*t;

    begin
      Add(res,tt);
      Clear(tt);
      continue := true;
    end Permute_Term;
    procedure Permute_Terms is new Visiting_Iterator(Permute_Term);

  begin
    Permute_Terms(s);
    return res;
  end "*";

  function "*" ( p : Permutation; s : QuadDobl_Complex_Laurentials.Poly )
               return QuadDobl_Complex_Laurentials.Poly is

    use QuadDobl_Complex_Laurentials;
    res : Poly := Null_Poly;

    procedure Permute_Term ( t : in Term; continue : out boolean ) is

      tt : Term := p*t;

    begin
      Add(res,tt);
      Clear(tt);
      continue := true;
    end Permute_Term;
    procedure Permute_Terms is new Visiting_Iterator(Permute_Term);

  begin
    Permute_Terms(s);
    return res;
  end "*";

  function "*" ( s : Standard_Complex_Poly_Systems.Poly_Sys;
                 p : Permutation )
               return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(s'range);

  begin
    for k in res'range loop
      res(k) := p*s(k);
    end loop;
    return res;
  end "*";

  function "*" ( s : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 p : Permutation )
               return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(s'range);

  begin
    for k in res'range loop
      res(k) := p*s(k);
    end loop;
    return res;
  end "*";

  function "*" ( s : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 p : Permutation )
               return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(s'range);

  begin
    for k in res'range loop
      res(k) := p*s(k);
    end loop;
    return res;
  end "*";

  function "*" ( s : Standard_Complex_Laur_Systems.Laur_Sys;
                 p : Permutation )
               return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(s'range);

  begin
    for k in res'range loop
      res(k) := p*s(k);
    end loop;
    return res;
  end "*";

  function "*" ( s : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 p : Permutation )
               return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    res : DoblDobl_Complex_Laur_Systems.Laur_Sys(s'range);

  begin
    for k in res'range loop
      res(k) := p*s(k);
    end loop;
    return res;
  end "*";

  function "*" ( s : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 p : Permutation )
               return QuadDobl_Complex_Laur_Systems.Laur_Sys is

    res : QuadDobl_Complex_Laur_Systems.Laur_Sys(s'range);

  begin
    for k in res'range loop
      res(k) := p*s(k);
    end loop;
    return res;
  end "*";

  function "*" ( p : Permutation;
                 s : Standard_Complex_Poly_Systems.Poly_Sys )
               return Standard_Complex_Poly_Systems.Poly_Sys is

    r : Standard_Complex_Poly_Systems.Poly_Sys(s'range);
    use Standard_Complex_Polynomials;

  begin
    for i in p'range loop
      if p(i) >= 0
       then Copy(s(p(i)),r(i));
       else r(i) := -s(-p(i));
      end if;
    end loop;
    return r;
  end "*";

  function "*" ( p : Permutation;
                 s : DoblDobl_Complex_Poly_Systems.Poly_Sys )
               return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    r : DoblDobl_Complex_Poly_Systems.Poly_Sys(s'range);
    use DoblDobl_Complex_Polynomials;

  begin
    for i in p'range loop
      if p(i) >= 0
       then Copy(s(p(i)),r(i));
       else r(i) := -s(-p(i));
      end if;
    end loop;
    return r;
  end "*";

  function "*" ( p : Permutation;
                 s : QuadDobl_Complex_Poly_Systems.Poly_Sys )
               return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    r : QuadDobl_Complex_Poly_Systems.Poly_Sys(s'range);
    use QuadDobl_Complex_Polynomials;

  begin
    for i in p'range loop
      if p(i) >= 0
       then Copy(s(p(i)),r(i));
       else r(i) := -s(-p(i));
      end if;
    end loop;
    return r;
  end "*";

  function "*" ( p : Permutation;
                 s : Standard_Complex_Laur_Systems.Laur_Sys )
               return Standard_Complex_Laur_Systems.Laur_Sys is

    r : Standard_Complex_Laur_Systems.Laur_Sys(s'range);
    use Standard_Complex_Laurentials;

  begin
    for i in p'range loop
      if p(i) >= 0
       then Copy(s(p(i)),r(i));
       else r(i) := -s(-p(i));
      end if;
    end loop;
    return r;
  end "*";

  function "*" ( p : Permutation;
                 s : DoblDobl_Complex_Laur_Systems.Laur_Sys )
               return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    r : DoblDobl_Complex_Laur_Systems.Laur_Sys(s'range);
    use DoblDobl_Complex_Laurentials;

  begin
    for i in p'range loop
      if p(i) >= 0
       then Copy(s(p(i)),r(i));
       else r(i) := -s(-p(i));
      end if;
    end loop;
    return r;
  end "*";

  function "*" ( p : Permutation;
                 s : QuadDobl_Complex_Laur_Systems.Laur_Sys )
               return QuadDobl_Complex_Laur_Systems.Laur_Sys is

    r : QuadDobl_Complex_Laur_Systems.Laur_Sys(s'range);
    use QuadDobl_Complex_Laurentials;

  begin
    for i in p'range loop
      if p(i) >= 0
       then Copy(s(p(i)),r(i));
       else r(i) := -s(-p(i));
      end if;
    end loop;
    return r;
  end "*";

end Permute_Operations;
