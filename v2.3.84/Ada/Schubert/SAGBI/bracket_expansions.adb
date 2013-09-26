with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;

package body Bracket_Expansions is

-- AUXILIARIES :

  function Create_Term ( n,d,i,j : natural32 ) return Term is

  -- DESCRIPTION :
  --   Returns the term that corresponds with xij.

    res : Term;

  begin
    res.cf := Create(1.0);
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(d*n) => 0);
    res.dg(integer32((i-1)*d + j)) := 1;
    return res;
  end Create_Term;

  function Create_Localized_Term ( n,d,i,j,k : natural32 ) return Term is

  -- DESCRIPTION :
  --   Creates a term in localized variables:
  --     if i >= k, then xij equals either 1 if i=j, or 0 if i/=j.
  --   If i >= k, then the `dg' of the term on return is null.

    res : Term;

  begin
    if i < k then
      res.cf := Create(1.0);
      res.dg := new Standard_Natural_Vectors.Vector'
                      (1..integer32(d*(n-d)) => 0);
      res.dg(integer32((i-1)*d + j)) := 1;
    else
      if i = j+k-1
       then res.cf := Create(1.0);
       else res.cf := Create(0.0);
      end if;
      return res;
    end if;
    return res;
  end Create_Localized_Term;

  function Create_Term ( loc,n,d,i,j : natural32 ) return Term is

  -- DESCRIPTION :
  --   Returns the term that corresponds with xij, with respect to
  --   the localization information in loc, which is 0,1, or 2.

    res : Term;

  begin
    if loc = 0
     then res.cf := Create(0.0);
     else res.cf := Create(1.0);
    end if;
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(d*n) => 0);
    if loc = 2
     then res.dg(integer32((i-1)*d + j)) := 1;
    end if;
    return res;
  end Create_Term;

-- TARGET ROUTINES :

  function Expand2 ( n,d : natural32; b : Bracket ) return Poly is

  -- DESCRIPTION :
  --   Does the expansion in the two-dimensional case.

    res,subtmp : Poly;
    b11 : Term := Create_Term(n,d,b(1),1);
    b12 : Term := Create_Term(n,d,b(1),2);
    b21 : Term := Create_Term(n,d,b(2),1);
    b22 : Term := Create_Term(n,d,b(2),2);
    d1 : Poly := Create(b11);
    d2 : Poly := Create(b21);

  begin                      -- res := b11*b22 - b21*b12;
    res := b22*d1;
    subtmp := b12*d2;
    Sub(res,subtmp);
    Clear(subtmp);
    Clear(b22); Clear(b12);
    Clear(b11); Clear(b21);
    Clear(d1);  Clear(d2);
    return res;
  end Expand2;

  function Expand ( n,d : natural32; b : Bracket ) return Poly is

    res : Poly;

  begin
    if b'last = 2 then
      res := Expand2(n,d,b);
    else
      res := Null_Poly;
      declare
        sig : integer := +1;
        b1 : Bracket(1..b'last-1);
      begin
        b1(1..b'last-1) := b(1..b'last-1);
        if b'last mod 2 = 0
         then sig := -1;
        end if;
        for i in b'range loop
          declare
            xt : Term := Create_Term(n,d,b(i),natural32(b'last));
            expb1,xtexpb1 : Poly;
          begin
            b1(1..i-1) := b(1..i-1);
            b1(i..b1'last) := b(i+1..b'last);
            expb1 := Expand(n,d,b1);
            xtexpb1 := xt*expb1;
            if sig > 0
             then Add(res,xtexpb1);
             else Sub(res,xtexpb1);
            end if;
            sig := -sig;
            Clear(expb1); Clear(xt); Clear(xtexpb1);
          end;
        end loop;
      end;
    end if;
    return res;
  end Expand;

  function Expand ( n,d : natural32; b : Bracket_Monomial ) return Poly is

    res : Poly := Null_Poly;
    fst : boolean := true;

    procedure Expand_Bracket ( bb : in Bracket; continue : out boolean ) is

      expbb : Poly;

    begin
      if fst then
        res := Expand(n,d,bb);
        fst := false;
      else
        expbb := Expand(n,d,bb);
        Mul(res,expbb);
        Clear(expbb);
      end if;
      continue := true;
    end Expand_Bracket;
    procedure Expand_Brackets is new Enumerate_Brackets(Expand_Bracket);

  begin
    Expand_Brackets(b);
    return res;
  end Expand;

  function Expand ( n,d : natural32; b : Bracket_Term ) return Poly is

    res : Poly := Expand(n,d,b.monom);

  begin
    Mul(res,b.coeff);
    return res;
  end Expand;

  function Expand ( n,d : natural32; b : Bracket_Polynomial ) return Poly is

    res : Poly := Null_Poly;

    procedure Expand_Term ( t : in Bracket_Term; continue : out boolean ) is

      expt : Poly := Expand(n,d,t);

    begin
      Add(res,expt);
      Clear(expt);
      continue := true;
    end Expand_Term;
    procedure Expand_Terms is new Enumerate_Terms(Expand_Term);

  begin
    Expand_Terms(b);
    return res;
  end Expand;

  function Localized_Expand2
             ( n,d : natural32; b : Bracket; k : natural32 ) return Poly is

  -- DESCRIPTION :
  --   Computes the localized expansion for two dimensions.
  --   Exploits the fact that b(1) < b(2).

    res : Poly;
    b11 : Term := Create_Localized_Term(n,d,b(1),1,k);
    b12 : Term := Create_Localized_Term(n,d,b(1),2,k);
    b21 : Term := Create_Localized_Term(n,d,b(2),1,k);
    b22 : Term := Create_Localized_Term(n,d,b(2),2,k);

  begin
    if b(2) < k then
      declare
        subtmp : Poly;
        d1 : Poly := Create(b11);
        d2 : Poly := Create(b21);
      begin                      -- res := b11*b22 - b21*b12;
        res := b22*d1;
        subtmp := b12*d2;
        Sub(res,subtmp);
        Clear(subtmp);
        Clear(b22); Clear(b12);
        Clear(b11); Clear(b21);
        Clear(d1);  Clear(d2);
      end;
    else
      if b(1) < k then
        b11.cf := b11.cf*b22.cf;
        b12.cf := b12.cf*b21.cf;
        res := Create(b11);
        Sub(res,b12);
        Clear(b11); Clear(b12);
      else
        declare
          rt : Term;
          dd : constant natural32 := d*(n-d);
        begin
          rt.cf := b11.cf*b22.cf - b21.cf*b12.cf;
          rt.dg := new Standard_Natural_Vectors.Vector'(1..integer32(dd) => 0);
          res := Create(rt);
        end;
      end if;
    end if;
    return res;
  end Localized_Expand2;

  function Localized_Expand ( n,d : natural32; b : Bracket ) return Poly is

    res : Poly;
    k : constant natural32 := n-d+1;

  begin
    if b'last = 2 then
      res := Localized_Expand2(n,d,b,k);
    else
      res := Null_Poly;
      declare
        sig : integer32 := +1;
        b1 : Bracket(1..b'last-1);
      begin
        b1(1..b'last-1) := b(1..b'last-1);
        if b'last mod 2 = 0
         then sig := -1;
        end if;
        for i in b'range loop
          declare
            xt : Term := Create_Localized_Term(n,d,b(i),natural32(b'last),k);
            expb1,xtexpb1 : Poly;
          begin
            if xt.cf /= Create(0.0) then
              b1(1..i-1) := b(1..i-1);
              b1(i..b1'last) := b(i+1..b'last);
              expb1 := Localized_Expand(n,d,b1);
              xtexpb1 := xt*expb1;
              if sig > 0
               then Add(res,xtexpb1);
               else Sub(res,xtexpb1);
              end if;
              Clear(expb1); Clear(xt); Clear(xtexpb1);
            end if;
            sig := -sig;
          end;
        end loop;
      end;
    end if;
    return res;
  end Localized_Expand;

  function Localization_Map ( n,d : natural32 ) return Matrix is

    res : Matrix(1..integer32(n),1..integer32(d));
    row,col : integer32;

  begin
    for i in res'range(1) loop            -- initialization
      for j in res'range(2) loop
        res(i,j) := 0;                    -- means empty space
      end loop;
    end loop;
    col := 0;
    for i in res'range(1) loop            -- one free element in every row
      col := col+1;
      if col > integer32(d)
       then col := 1;
      end if;
      res(i,col) := 2;
    end loop;
    row := 0;
    for j in res'range(2) loop            -- one 1 in every column
      loop
        row := row+1;
        if row > integer32(n)
         then row := 1;
        end if;
        exit when (res(row,j) = -1);      -- empty space found
      end loop;
      res(row,j) := 1;
    end loop;
    for j in res'range(2) loop            -- fill in d-1 zeros in every column
      for i in 1..integer32(d-1) loop
        row := 0;
        loop
          row := row+1;
          if row > integer32(n)
           then row := 1;
          end if;
          exit when (res(row,j) = -1);    -- empty space found
        end loop;
        res(row,j) := 0;
      end loop;
    end loop;
    for i in res'range(1) loop            -- fill rest with free elements
      for j in res'range(2) loop
        if res(i,j) = -1
         then res(i,j) := 2;
        end if;
      end loop;
    end loop;
    return res;
  end Localization_Map;

  function Expand2 ( locmap : Matrix; b : Bracket ) return Poly is

  -- DESCRIPTION :
  --   Expands a 2-by-2 minor of the matrix, selecting the rows with
  --   entries in the 2-bracket b, respecting the localization map.

    res,subtmp : Poly;
    n : constant natural32 := natural32(locmap'length(1));
    d : constant natural32 := natural32(locmap'length(2));
    b11 : Term := Create_Term(locmap(integer32(b(1)),1),n,d,b(1),1);
    b12 : Term := Create_Term(locmap(integer32(b(1)),2),n,d,b(1),2);
    b21 : Term := Create_Term(locmap(integer32(b(2)),1),n,d,b(2),1);
    b22 : Term := Create_Term(locmap(integer32(b(2)),2),n,d,b(2),2);
    d1 : Poly := Create(b11);
    d2 : Poly := Create(b21);

  begin                      -- res := b11*b22 - b21*b12;
    res := b22*d1;
    subtmp := b12*d2;
    Sub(res,subtmp);
    Clear(subtmp);
    Clear(b22); Clear(b12);
    Clear(b11); Clear(b21);
    Clear(d1);  Clear(d2);
    return res;
  end Expand2;

  function Expand ( locmap : Matrix; b : Bracket ) return Poly is

  -- DESCRIPTION :
  --   Expands a d-by-d minor of the matrix, selecting the rows with
  --   entries in b, respecting the localization map.
  --   The expansion starts at the last column and proceeds recursively.

    res : Poly;
    n : constant natural32 := natural32(locmap'length(1));
    d : constant natural32 := natural32(locmap'length(2));

  begin
    if b'last = 2 then
      res := Expand2(locmap,b);
    else
      res := Null_Poly;
      declare
        sig : integer := +1;
        b1 : Bracket(1..b'last-1);
      begin
        b1(1..b'last-1) := b(1..b'last-1);
        if b'last mod 2 = 0
         then sig := -1;
        end if;
        for i in b'range loop
          declare
            xt : Term
               := Create_Term(locmap(integer32(b(i)),b'last),n,d,
                              b(i),natural32(b'last));
            expb1,xtexpb1 : Poly;
          begin
            if xt.cf /= Create(0.0) then
              b1(1..i-1) := b(1..i-1);
              b1(i..b1'last) := b(i+1..b'last);
              expb1 := Expand(locmap,b1);
              xtexpb1 := xt*expb1;
              if sig > 0
               then Add(res,xtexpb1);
               else Sub(res,xtexpb1);
              end if;
              Clear(expb1); Clear(xtexpb1);
            end if;
            Clear(xt);
          end;
          sig := -sig;
        end loop;
      end;
    end if;
    return res;
  end Expand;

  function Reduce_Variables
             ( locmap : Matrix; dg : Standard_Natural_Vectors.Vector )
             return Standard_Natural_Vectors.Vector is

  -- DESCRIPTION :
  --   Removes the #variables in the degrees vectors, removing all variables
  --   that correspond to zeros in the localization map.

    res : Standard_Natural_Vectors.Vector(dg'range) := dg;
    cnt : integer32 := res'last;
    d : constant integer32 := locmap'length(2);

  begin
    for i in reverse locmap'range(1) loop
      for j in reverse locmap'range(2) loop
        if locmap(i,j) /= 2 then
          cnt := cnt - 1;
          for k in ((i-1)*d+j)..cnt loop
            res(k) := res(k+1);
          end loop;
        end if;
      end loop;
    end loop;
    return res(res'first..cnt);
  end Reduce_Variables;

  function Reduce_Variables ( locmap : Matrix; t : Term ) return Term is

  -- DESCRIPTION :
  --   Reduces the #variables in the term, removing all variables that
  --   correspond to zeros in the localization map.

    res : Term;

  begin
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector'
                    (Reduce_Variables(locmap,t.dg.all));
    return res;
  end Reduce_Variables;

  procedure Reduce_Variables ( locmap : in Matrix; p : in out Poly ) is

    rp : Poly := Null_Poly;

    procedure Reduce_Term ( t : in Term; continue : out boolean ) is

     rt : constant Term := Reduce_Variables(locmap,t);

    begin
      Add(rp,rt);
      continue := true;
    end Reduce_Term;
    procedure Reduce_Terms is new Visiting_Iterator(Reduce_Term);

  begin
    Reduce_Terms(p);
    Clear(p);
    p := rp;
  end Reduce_Variables;

end Bracket_Expansions;
