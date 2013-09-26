with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;

package body Curves_into_Grassmannian is

-- CREATOR :

  function Symbolic_Create ( m,p,q : natural32; top,bottom : Bracket )
                           return Standard_Complex_Poly_Matrices.Matrix is

    res : Standard_Complex_Poly_Matrices.Matrix
            (1..integer32(m+p),1..integer32(p));
    rws : constant natural32 := (m+p)*(q+1);
    n : constant integer32 := integer32(Number_of_Variables(top,bottom) + 2);
    row,ind : integer32;
    s_deg,t_deg : natural32;
    t : Term;

  begin
    for i in res'range(1) loop                    -- initialization
      for j in res'range(2) loop
        res(i,j) := Null_Poly;
      end loop;
    end loop;
    t.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    t.cf := Create(1.0);
    ind := 0;                                     -- ind counts #variables
    for j in 1..integer32(p) loop                 -- assign columnwise
      t_deg := (bottom(j)-1)/(m+p);               -- degree in t to homogenize
      row := 0; s_deg := 0;
      for i in 1..rws loop
        row := row + 1;
        if i >= top(j) and i <= bottom(j) then
          ind := ind+1;
          t.dg(n-1) := s_deg; t.dg(n) := t_deg;
          t.dg(ind) := 1;
          Add(res(row,j),t);
          t.dg(n-1) := 0; t.dg(n) := 0;
          t.dg(ind) := 0;
        end if;
        if i mod (m+p) = 0 then
          row := 0; s_deg := s_deg+1;
          if t_deg > 0
           then t_deg := t_deg-1;
          end if;
        end if;
      end loop;
    end loop;
    Clear(t);
    return res;
  end Symbolic_Create;

-- SELECTORS :

  function Number_of_Variables ( top,bottom : Bracket ) return natural32 is

    cnt : natural32 := 0;

  begin
    for j in top'range loop
      cnt := cnt + (bottom(j) - top(j) + 1);
    end loop;
    return cnt;
  end Number_of_Variables;

  function Standard_Coordinate_Frame
             ( m,p,q : natural32; top,bottom : Bracket;
               coeff : Standard_Complex_Matrices.Matrix )
             return Standard_Natural_Matrices.Matrix is

    rws : constant integer32 := integer32((m+p)*(q+1));
    res : Standard_Natural_Matrices.Matrix(1..rws,1..integer32(p));
    tol : constant double_float := 10.0**(-10);
    first : boolean;

  begin
    for j in 1..integer32(p) loop
      first := true;
      for i in 1..rws loop
        if i < integer32(top(j)) or i > integer32(bottom(j)) then
          res(i,j) := 0;
        elsif (first and (AbsVal(coeff(i,j)) > tol)) then
          res(i,j) := 1; first := false;
        else
          res(i,j) := 2;
        end if;
      end loop;
    end loop;
    return res;
  end Standard_Coordinate_Frame;

  function Eval ( c : Term; s,t : Complex_Number ) return Term is

    res : Term;

  begin
    Copy(c,res);
    for i in 1..res.dg(res.dg'last-1) loop        -- evaluate s
      res.cf := res.cf*s;
    end loop;
    res.dg(res.dg'last-1) := 0;
    for i in 1..res.dg(res.dg'last) loop          -- evaluate t
      res.cf := res.cf*t;
    end loop;
    res.dg(res.dg'last) := 0;
    return res;
  end Eval;

  function Eval ( p : Poly; s,t : Complex_Number ) return Poly is

    res : Poly := Null_Poly;

    procedure Eval_Term ( ct : in Term; continue : out boolean ) is

      et : Term := Eval(ct,s,t);

    begin
      Add(res,et);
      Clear(et);
      continue := true;
    end Eval_Term;
    procedure Eval_Terms is new Visiting_Iterator(Eval_Term);

  begin
    Eval_Terms(p);
    return res;
  end Eval;

  function Eval ( c : Standard_Complex_Poly_Matrices.Matrix;
                  s,t : Complex_Number )
                return Standard_Complex_Poly_Matrices.Matrix is

    res : Standard_Complex_Poly_Matrices.Matrix(c'range(1),c'range(2));

  begin
    for i in c'range(1) loop
      for j in c'range(2) loop
        if c(i,j) = Null_Poly
         then res(i,j) := Null_Poly;
         else res(i,j) := Eval(c(i,j),s,t);
        end if;
      end loop;
    end loop;
    return res;
  end Eval;

  function Elim ( c : Term; s,t : Complex_Number ) return Term is

    res : Term;

  begin
    res.dg := new Standard_Natural_Vectors.Vector'
                    (c.dg(c.dg'first..c.dg'last-2));
    res.cf := c.cf;
    for i in 1..c.dg(c.dg'last-1) loop        -- evaluate s
      res.cf := res.cf*s;
    end loop;
    for i in 1..c.dg(c.dg'last) loop          -- evaluate t
      res.cf := res.cf*t;
    end loop;
    return res;
  end Elim;

  function Elim ( p : Poly; s,t : Complex_Number ) return Poly is

    res : Poly := Null_Poly;

    procedure Elim_Term ( ct : in Term; continue : out boolean ) is

      et : Term := Elim(ct,s,t);

    begin
      Add(res,et);
      Clear(et);
      continue := true;
    end Elim_Term;
    procedure Elim_Terms is new Visiting_Iterator(Elim_Term);

  begin
    Elim_Terms(p);
    return res;
  end Elim;

  function Elim ( c : Standard_Complex_Poly_matrices.Matrix;
                  s,t : Complex_Number )
                return Standard_Complex_Poly_Matrices.Matrix is

    res : Standard_Complex_Poly_Matrices.Matrix(c'range(1),c'range(2));

  begin
    for i in c'range(1) loop
      for j in c'range(2) loop
        if c(i,j) = Null_Poly
         then res(i,j) := Null_Poly;
         else res(i,j) := Elim(c(i,j),s,t);
        end if;
      end loop;
    end loop;
    return res;
  end Elim;

  function Substitute ( t : Term; v : Standard_Complex_Vectors.Vector )
                      return Term is

    res : Term;
    n : constant integer32 := v'length;

  begin
    if n > t.dg'last then
      return t;
    else
      res.dg := new Standard_Natural_Vectors.Vector(1..t.dg'last-n);
      for i in res.dg'range loop
        res.dg(i) := t.dg(i+n);
      end loop;
      res.cf := t.cf;
      for i in v'range loop
        for j in 1..t.dg(i) loop
          res.cf := res.cf*v(i);
        end loop;
      end loop;
    end if;
    return res;
  end Substitute;

  function Substitute ( p : Poly; v : Standard_Complex_Vectors.Vector )
                      return Poly is

    res : Poly := Null_Poly;

    procedure Substitute_Term ( t : in Term; continue : out boolean ) is

      st : Term := Substitute(t,v);

    begin
      Add(res,st);
      Clear(st);
      continue := true;
    end Substitute_Term;
    procedure Substitute_Terms is new Visiting_Iterator(Substitute_Term);

  begin
    Substitute_Terms(p);
    return res;
  end Substitute;

  function Substitute ( c : Standard_Complex_Poly_Matrices.Matrix;
                        v : Standard_Complex_Vectors.Vector )
                      return Standard_Complex_Poly_Matrices.Matrix is

    res : Standard_Complex_Poly_Matrices.Matrix(c'range(1),c'range(2));

  begin
    for i in c'range(1) loop
      for j in c'range(2) loop
        if c(i,j) = Null_Poly
         then res(i,j) := Null_Poly;
         else res(i,j) := Substitute(c(i,j),v);
        end if;
      end loop;
    end loop;
    return res;
  end Substitute;

  function Convert ( p : Poly ) return Complex_Number is

    res : Complex_Number := Create(0.0);

    procedure Convert_Term ( t : in Term; continue : out boolean ) is
    begin
      res := t.cf;
      continue := false;
    end Convert_Term;
    procedure Convert_Terms is new Visiting_Iterator(Convert_Term);

  begin
    Convert_Terms(p);
    return res;
  end Convert;

  function Convert ( c : Standard_Complex_Poly_Matrices.Matrix )
                   return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(c'range(1),c'range(2));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        if c(i,j) = Null_Poly
         then res(i,j) := Create(0.0);
         else res(i,j) := Convert(c(i,j));
        end if;
      end loop;
    end loop;
    return res;
  end Convert;

  function Column_Localize ( top,bottom : Bracket;
                             locmap : Standard_Natural_Matrices.Matrix;
                             t : Term ) return Term is

  -- DESCRIPTION :
  --   Applies the localization map to the term, eliminating those xij's
  --   xij for which the corresponding entry in locmap is either 0 or 1.

  -- NOTE : 
  --   This localization assumes that t.dg(k) = 0 with k for which the
  --   corresponding (i,j) with locmap(i,j) = 0.
  --   The localization pattern is traversed columnwise.

    res : Term;
    ndg : Standard_Natural_Vectors.Vector(t.dg'range);
    cnt : integer32 := t.dg'first-1;
    ind : integer32 := cnt;

  begin
    for j in locmap'range(2) loop       -- columnwise order of the variables
      for i in top(j)..bottom(j) loop   -- restricted range skips the zeros
        ind := ind + 1;
        if locmap(integer32(i),j) = 2 then         -- skip the ones
          cnt := cnt + 1;
          ndg(cnt) := t.dg(ind);
        end if;
      end loop;
    end loop;
    for i in ind+1..t.dg'last loop      -- leave the lifting !
      cnt := cnt + 1;
      ndg(cnt) := t.dg(i);
    end loop;
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector'(ndg(1..cnt));
    return res;
  end Column_Localize;

  function Column_Localize ( top,bottom : Bracket;
                             locmap : Standard_Natural_Matrices.Matrix;
                             p : Poly ) return Poly is

  -- DESCRIPTION :
  --   Applies the localization map to the polynomial, eliminating
  --   those xij's for which locmap(i,j) is either 0 or 1.

    res : Poly := Null_Poly;

    procedure Column_Localize_Term ( t : in Term; continue : out boolean ) is

      lt : Term := Column_Localize(top,bottom,locmap,t);

    begin
      Add(res,lt);
      Clear(lt.dg);
      continue := true;
    end Column_Localize_Term;
    procedure Column_Localize_Terms is
      new Visiting_Iterator(Column_Localize_Term);

  begin
    Column_Localize_Terms(p);
    return res;
  end Column_Localize;

  function Column_Localize ( top,bottom : Bracket;
                             locmap : Standard_Natural_Matrices.Matrix;
                             pm : Standard_Complex_Poly_Matrices.Matrix )
                           return Standard_Complex_Poly_Matrices.Matrix is
  
    res : Standard_Complex_Poly_Matrices.Matrix(pm'range(1),pm'range(2));

  begin
    for i in pm'range(1) loop
      for j in pm'range(2) loop
        res(i,j) := Column_Localize(top,bottom,locmap,pm(i,j));
      end loop;
    end loop;
    return res;
  end Column_Localize;

  function Column_Localize ( top,bottom : Bracket;
                             locmap : Standard_Natural_Matrices.Matrix;
                             p : Poly_Sys ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Column_Localize(top,bottom,locmap,p(i));
    end loop;
    return res;
  end Column_Localize;

  function Column_Vector_Rep
             ( top,bottom : Bracket;
               cffmat : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Vectors.Vector is

    dim : constant natural32 := Number_of_Variables(top,bottom);
    res : Standard_Complex_Vectors.Vector(1..integer32(dim));
    cnt : integer32 := 0;

  begin
    for j in cffmat'range(2) loop
      for i in integer32(top(j))..integer32(bottom(j)) loop
        cnt := cnt + 1;
        res(cnt) := cffmat(i,j);
      end loop;
    end loop;
    return res;
  end Column_Vector_Rep;

  function Column_Vector_Rep ( locmap : Standard_Natural_Matrices.Matrix;
                               cffmat : Standard_Complex_Matrices.Matrix )
                             return Standard_Complex_Vectors.Vector is

    dim : constant integer32 := cffmat'length(1)*cffmat'length(2);
    res : Standard_Complex_Vectors.Vector(1..dim);
    cnt : integer32 := 0;

  begin
    for j in cffmat'range(2) loop
      for i in cffmat'range(1) loop
        if locmap(i,j) = 2 then
          cnt := cnt + 1;
          res(cnt) := cffmat(i,j);
        end if;
      end loop;
    end loop;
    return res(1..cnt);
  end Column_Vector_Rep;

  function Column_Matrix_Rep
              ( locmap : Standard_Natural_Matrices.Matrix;
                cffvec : Standard_Complex_Vectors.Vector )
              return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(locmap'range(1),locmap'range(2));
    cnt : integer32 := 0;

  begin
    for j in locmap'range(2) loop
      for i in locmap'range(1) loop
        if locmap(i,j) = 0 then
          res(i,j) := Create(0.0);
        elsif locmap(i,j) = 1 then
          res(i,j) := Create(1.0);
        else
          cnt := cnt + 1;
          res(i,j) := cffvec(cnt);
        end if;
      end loop;
    end loop;
    return res;
  end Column_Matrix_Rep;

  procedure Swap ( p : in out Poly; k,l : in integer32 ) is

    procedure Swap_Term ( t : in out Term; continue : out boolean ) is

      lval : constant natural32 := t.dg(l);

    begin
      t.dg(l) := t.dg(k);
      t.dg(k) := lval;
      continue := true;
    end Swap_Term;
    procedure Swap_Terms is new Changing_Iterator(Swap_Term);

  begin
    Swap_Terms(p);
  end Swap;

  procedure Swap ( c : in out Standard_Complex_Poly_Matrices.Matrix;
                   k,l : in integer32 ) is
  begin
    for i in c'range(1) loop
      for j in c'range(2) loop
        if c(i,j) /= Null_Poly
         then Swap(c(i,j),k,l);
        end if;
      end loop;
    end loop;
  end Swap;

  function Insert ( p : Poly; k : integer32 ) return Poly is

    res : Poly := Null_Poly;

    procedure Insert_Term ( t : in Term; continue : out boolean ) is

      rt : Term;

    begin
      rt.cf := t.cf;
      rt.dg := new Standard_Natural_Vectors.Vector(t.dg'first..t.dg'last+1);
      rt.dg(t.dg'first..k-1) := t.dg(t.dg'first..k-1);
      rt.dg(k) := 0;
      rt.dg(k+1..rt.dg'last) := t.dg(k..t.dg'last);
      Add(res,rt);
      Clear(rt);
      continue := true;
    end Insert_Term;
    procedure Insert_Terms is new Visiting_Iterator(Insert_Term);

  begin
    Insert_Terms(p);
    return res;
  end Insert;

  function Insert ( c : Standard_Complex_Poly_Matrices.Matrix; k : integer32 )
                  return Standard_Complex_Poly_Matrices.Matrix is

    res : Standard_Complex_Poly_Matrices.Matrix(c'range(1),c'range(2));

  begin
    for i in c'range(1) loop
      for j in c'range(2) loop
        if c(i,j) = Null_Poly
         then res(i,j) := Null_Poly;
         else res(i,j) := Insert(c(i,j),k);
        end if;
      end loop;
    end loop;
    return res;
  end Insert;

end Curves_into_Grassmannian;
