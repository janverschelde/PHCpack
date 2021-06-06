with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Characters_and_Numbers;
with Standard_Complex_Numbers;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Random_Vectors;
with Symbol_Table;
with Double_Laurent_Series;
with Random_Laurent_Series;             use Random_Laurent_Series;

package body Double_Lseries_Polynomials is

  procedure Write ( plead : in Standard_Integer_Vectors.Vector;
                    pcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                    pmons : in Standard_Integer_VecVecs.VecVec;
                    s : in string := "p" ) is
  begin
    for k in plead'range loop
      put(s & "("); put(k,1); put(") :"); put(pmons(k)); new_line;
      Double_Laurent_Series.Write(plead(k),pcffs(k).all);
    end loop;
  end Write;

  procedure Write ( tab : in Table; s : in string := "p" ) is
  begin
    Write(tab.lead.all,tab.cffs,tab.mons.all,s);
  end Write;

  procedure Write ( tab : in Table_Vector; s : in string := "p" ) is
  begin
    for i in 1..tab.nbt loop
      declare
        stri : constant string := Characters_and_Numbers.Convert(i);
      begin
        Write(tab.lead(i).all,tab.cffs(i),tab.mons(i).all,s & stri);
      end;
    end loop;
  end Write;

  procedure Write ( tva : in Table_Vector_Array; s : in string := "p" ) is

    ltv : Link_to_Table_Vector;

  begin
    for i in tva'range loop
      ltv := tva(i); 
      for j in 1..ltv.nbt loop
        declare
          stri : constant string := Characters_and_Numbers.Convert(i);
          strj : constant string := Characters_and_Numbers.Convert(j);
          tab : Table;
        begin
          tab.nvr := ltv.nvr;
          tab.lead := ltv.lead(j);
          tab.cffs := ltv.cffs(j);
          tab.mons := ltv.mons(j);
          tab.nbt := tab.lead'last;
          Write(tab,s & stri & "," & strj);
        end;
      end loop;
    end loop;
  end Write;

-- BASIC OPERATIONS :

  procedure Make_Random_Polynomial
              ( dim,nbr,deg,pwr,low,upp : in integer32;
                lead : out Standard_Integer_Vectors.Vector;
                cffs : out Standard_Complex_VecVecs.Link_to_VecVec;
                mons : out Standard_Integer_VecVecs.VecVec ) is
  begin
    for k in 1..nbr loop
      declare
        mon : constant Standard_Integer_Vectors.Vector(1..dim)
            := Standard_Random_Vectors.Random_Vector(1,dim,0,pwr);
      begin
        mons(k) := new Standard_Integer_Vectors.Vector'(mon);
      end;
    end loop;
    Random_Vector(nbr,deg,low,upp,lead,cffs);
  end Make_Random_Polynomial;

  function Random_Table ( dim,nbr,deg,pwr,low,upp : integer32 ) return Table is

    res : Table;
    lead : Standard_Integer_Vectors.Vector(1..nbr);
    mons : Standard_Integer_VecVecs.VecVec(1..nbr);

  begin
    res.nbt := nbr;
    res.nvr := dim;
    Make_Random_Polynomial(dim,nbr,deg,pwr,low,upp,lead,res.cffs,mons);
    res.lead := new Standard_Integer_Vectors.Vector'(lead);
    res.mons := new Standard_Integer_VecVecs.VecVec'(mons);
    return res;
  end Random_Table;

-- EVALUATORS :

  procedure Eval ( deg,mlead : in integer32;
                   cff : in Standard_Complex_Vectors.Link_to_Vector;
                   mon : in Standard_Integer_Vectors.Link_to_Vector;
                   xlead : in Standard_Integer_Vectors.Vector;
                   xcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                   ye : out integer32;
                   yc : out Standard_Complex_Vectors.Vector ) is

    ze : integer32;
    zc : Standard_Complex_Vectors.Vector(0..deg);

  begin
    ye := mlead;
    for k in 0..deg loop -- initialize result with monomial coefficient
      yc(k) := cff(k);
    end loop;
    for i in mon'range loop -- mon(i) is the power of the i-th variable
      if mon(i) > 0 then
        for j in 1..mon(i) loop
          Double_Laurent_Series.Multiply
            (deg,ye,xlead(i),yc,xcffs(i).all,ze,zc);
          ye := ze;
          for k in 0..deg loop
            yc(k) := zc(k);
          end loop;
        end loop;
      end if;
    end loop;
  end Eval;

  procedure Eval ( deg : in integer32;
                   plead : in Standard_Integer_Vectors.Vector;
                   pcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                   pmons : in Standard_Integer_VecVecs.VecVec;
                   xlead : in Standard_Integer_Vectors.Vector;
                   xcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                   ye : out integer32;
                   yc : out Standard_Complex_Vectors.Vector;
                   verbose : in boolean := true ) is

    ze,ewrk : integer32;
    zc,cwrk : Standard_Complex_Vectors.Vector(0..deg);

  begin
    if plead'last < plead'first then
      ye := 0;
      for k in 0..deg loop
        yc(k) := Standard_Complex_Numbers.Create(0.0);
      end loop;
    else
      Eval(deg,plead(1),pcffs(1),pmons(1),xlead,xcffs,ye,yc);
      if verbose then
        put_line("After evaluating the first monomial :");
        Double_Laurent_Series.Write(ye,yc);
      end if;
      for i in 2..plead'last loop
        Eval(deg,plead(i),pcffs(i),pmons(i),xlead,xcffs,ze,zc);
        if verbose then
          put("The value of monomial : "); put(i,1); put_line(" :");
          Double_Laurent_Series.Write(ze,zc);
        end if;
        Double_Laurent_Series.Add(deg,ye,ze,yc,zc,ewrk,cwrk);
        ye := ewrk;
        for k in 0..deg loop
          yc(k) := cwrk(k);
        end loop;
        if verbose then
          put("After update with the value of monomial ");
          put(i,1); put_line(" :");
          Double_Laurent_Series.Write(ye,yc);
        end if;
      end loop;
    end if;
  end Eval;

  procedure Eval ( deg : in integer32; tab : in Table;
                   xlead : in Standard_Integer_Vectors.Vector;
                   xcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                   ye : out integer32;
                   yc : out Standard_Complex_Vectors.Vector;
                   verbose : in boolean := true ) is

  begin
    Eval(deg,tab.lead.all,tab.cffs,tab.mons.all,xlead,xcffs,ye,yc,verbose);
  end Eval;

  procedure Eval ( deg : in integer32; tab : in Table_Vector;
                   xlead : in Standard_Integer_Vectors.Vector;
                   xcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                   ylead : out Standard_Integer_Vectors.Vector;
                   ycffs : out Standard_Complex_VecVecs.VecVec;
                   verbose : in boolean := true ) is

     use Standard_Complex_Vectors;
     cff : Standard_Complex_Vectors.Link_to_Vector;
     zero : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Complex_Numbers.Create(0.0);

  begin
    for i in 1..tab.nbt loop
      if ycffs(i) /= null then
        cff := ycffs(i);
      else
        declare
          v : constant Standard_Complex_Vectors.Vector(0..deg)
            := (0..deg => zero);
        begin
          ycffs(i) := new Standard_Complex_Vectors.Vector'(v);
          cff := ycffs(i);
        end;
      end if;
      Eval(deg,tab.lead(i).all,tab.cffs(i),tab.mons(i).all,
           xlead,xcffs,ylead(i),cff.all,verbose);
    end loop;
  end Eval;

  procedure Eval ( deg : in integer32; tab : in Table_Vector_Array;
                   xlead : in Standard_Integer_Vectors.Vector;
                   xcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                Alead : in out Standard_Integer_Matrices.Matrix;
                Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                verbose : in boolean := true ) is

    ltv : Link_to_Table_Vector;
    row : Standard_Complex_VecVecs.Link_to_VecVec;
    cff : Standard_Complex_Vectors.Link_to_Vector;

    use Standard_Complex_Vectors;

  begin
    for i in tab'range loop
      ltv := tab(i); 
      row := Acffs(i);
      for j in 1..ltv.nbt loop
        cff := row(j);
        Eval(deg,ltv.lead(j).all,ltv.cffs(j),ltv.mons(j).all,
             xlead,xcffs,Alead(i,j),cff.all,verbose);
      end loop;
    end loop;
  end Eval;

  function tsymbol_index return integer32 is

    res : integer32 := 0;
    nbr : constant natural32 := Symbol_Table.Number;

  begin
    for k in 1..nbr loop
      declare
        sb : constant Symbol_Table.Symbol := Symbol_Table.get(k);
      begin
        if sb(1) = 't'
         then res := integer32(k); exit;
        end if;
      end;
    end loop;
    return res;
  end tsymbol_index;

  function Index_of_Degrees
             ( mons : Standard_Integer_VecVecs.VecVec;
               idx : integer32;
               degs : Standard_Integer_Vectors.Vector ) return integer32 is

    monexp : Standard_Integer_Vectors.Link_to_Vector;
    allequal : boolean;

  begin
    for k in mons'first..(idx-1) loop
      monexp := mons(k);
      allequal := true;
      for i in monexp'range loop
        allequal := (monexp(i) = degs(i));
        exit when not allequal;
      end loop;
      if allequal
       then return k;
      end if;
    end loop;
    return idx;
  end Index_of_Degrees;

  function Minimum_Laurent_Series_Degree
             ( p : Poly; tdx : integer32 ) return integer32 is

    maxdegt : constant integer32 := Maximal_Degree(p,tdx);
    mindegt : constant integer32 := Minimal_Degree(p,tdx);

  begin
    return (maxdegt - mindegt - 1);
  end Minimum_Laurent_Series_Degree;

  function Minimum_Laurent_Series_Degree
             ( p : Laur_Sys; tdx : integer32 ) return integer32 is

    res : integer32 := Minimum_Laurent_Series_Degree(p(p'first),tdx);
    deg : integer32;

  begin
    for k in p'first+1..p'last loop
      deg := Minimum_Laurent_Series_Degree(p(k),tdx);
      if deg > res
       then res := deg;
      end if;
    end loop;
    return res;
  end Minimum_Laurent_Series_Degree;

  procedure Make_Series_Polynomial
              ( p : in Poly; dim,nvr,tdx,deg : in integer32;
                lead : out Standard_Integer_Vectors.Link_to_Vector;
                cffs : out Standard_Complex_VecVecs.Link_to_VecVec;
                mons : out Standard_Integer_VecVecs.Link_to_VecVec ) is

    nbr : constant integer32 := integer32(Number_of_Terms(p));
    plead : Standard_Integer_Vectors.Vector(1..nbr);
    wcffs : Standard_Complex_VecVecs.VecVec(1..nbr);
    pmons : Standard_Integer_VecVecs.VecVec(1..nbr);
    cnt : integer32 := 0;

    procedure Visit_Term ( t : in Term; continue : out boolean ) is

    -- DESCRIPTION :
    --   If tdx is zero, then the monomial is just copied,
    --   otherwise the terms in the Laurent series need to be
    --   collected from the expanded form.

      mon : Standard_Integer_Vectors.Vector(1..nvr);
      cff : Standard_Complex_Vectors.Vector(0..deg);
      located,newlead,gap : integer32;
      cfflocated : Standard_Complex_Vectors.Link_to_Vector;

    begin
      cnt := cnt + 1; 
      if tdx = 0 then          -- constant Laurent series
        plead(cnt) := 0;
        for k in 1..dim loop
          mon(k) := integer32(t.dg(k));
        end loop;
        pmons(cnt) := new Standard_Integer_Vectors.Vector'(mon);
        cff(0) := t.cf;
        for k in 1..deg loop
          cff(k) := Standard_Complex_Numbers.Create(0.0);
        end loop;
        wcffs(cnt) := new Standard_Complex_Vectors.Vector'(cff);
      else
        for k in 1..(tdx-1) loop
          mon(k) := integer32(t.dg(k));
        end loop;
        for k in (tdx+1)..dim loop
          mon(k-1) := integer32(t.dg(k));
        end loop;
        located := Index_of_Degrees(pmons,cnt,mon);
        if located = cnt then                        --  a new monomial
          plead(cnt) := integer32(t.dg(tdx));
          pmons(cnt) := new Standard_Integer_Vectors.Vector'(mon);
          cff(0) := t.cf;
          for k in 1..deg loop
            cff(k) := Standard_Complex_Numbers.Create(0.0);
          end loop;
          wcffs(cnt) := new Standard_Complex_Vectors.Vector'(cff);
        else
          cnt := cnt-1;                              -- no new monomial
          newlead := integer32(t.dg(tdx));
          cfflocated := wcffs(located);            -- fit in new t-term
          if newlead = plead(located) then
            Standard_Complex_Numbers.Add(cfflocated(0),t.cf);
          elsif newlead > plead(located) then
            gap := newlead - plead(located);   -- keep leading exponent
            if gap <= deg then
              Standard_Complex_Numbers.Add(cfflocated(gap),t.cf);
            -- else the new coefficient will be ignored ...
            end if;
          else -- newlead < plead(located) 
            gap := plead(located) - newlead; -- leading exponent changes!
            for k in reverse 0..deg-gap loop        -- shift coefficients
              cfflocated(gap+k) := cfflocated(k);
            end loop;
            for k in 1..(gap-1) loop
              exit when (k > deg);
              cfflocated(k) := Standard_Complex_Numbers.Create(0.0);
            end loop;
            cfflocated(0) := t.cf;
            plead(located) := newlead;
          end if;
        end if;
      end if;
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    lead := new Standard_Integer_Vectors.Vector'(plead(1..cnt));
    cffs := new Standard_Complex_VecVecs.VecVec'(wcffs(1..cnt));
    mons := new Standard_Integer_VecVecs.VecVec'(pmons(1..cnt));
  end Make_Series_Polynomial;

  function Make_Table ( p : Poly; dim,nvr,tdx,deg : integer32 ) return Table is

    res : Table;

  begin
    res.nbt := dim;
    res.nvr := nvr;
    Make_Series_Polynomial(p,dim,nvr,tdx,deg,res.lead,res.cffs,res.mons);
    return res;
  end Make_Table;

  procedure Make_Series_System
              ( p : in Laur_Sys; dim,nvr,tdx,deg : in integer32;
                lead : out Standard_Integer_VecVecs.VecVec;
                cffs : out Standard_Complex_VecVecs.Array_of_VecVecs;
                mons : out Standard_Integer_VecVecs.Array_of_VecVecs ) is
  begin
    for k in p'range loop
      Make_Series_Polynomial(p(k),dim,nvr,tdx,deg,lead(k),cffs(k),mons(k));
    end loop;
  end Make_Series_System;

  function Make_Table_Vector
             ( p : Laur_Sys; dim,nvr,tdx,deg : integer32 )
             return Table_Vector is

    neq : constant integer32 := p'last;
    res : Table_Vector(neq);

  begin
    res.nvr := nvr;
    Make_Series_System(p,dim,nvr,tdx,deg,res.lead,res.cffs,res.mons);
    return res;
  end Make_Table_Vector;

  function Make_Table_Vector_Array
             ( jp : Jaco_Mat; tdx,deg : integer32 )
             return Table_Vector_Array is

    res : Table_Vector_Array(jp'range(1));
    dim : constant integer32 := jp'length(2);
    nbt,nvr : integer32;

  begin
    if tdx = 0
     then nbt := dim;
     else nbt := dim-1;
    end if;
    nvr := nbt;
    for i in jp'range(1) loop
      declare
        tv : Table_Vector(nbt);
      begin
        tv.nvr := nvr;
        if tdx = 0 then
          for j in jp'range(2) loop
            declare
              tab : constant Table := Make_Table(jp(i,j),dim,nvr,tdx,deg);
            begin
              tv.lead(j) := new Standard_Integer_Vectors.Vector'(tab.lead.all);
              tv.cffs(j) := tab.cffs;
              tv.mons(j) := tab.mons;
            end;
          end loop;
        else
          for j in jp'first(2)..(tdx-1) loop
            declare
              tab : constant Table := Make_Table(jp(i,j),dim,nvr,tdx,deg);
            begin
              tv.lead(j) := new Standard_Integer_Vectors.Vector'(tab.lead.all);
              tv.cffs(j) := tab.cffs;
              tv.mons(j) := tab.mons;
            end;
          end loop;
          for j in (tdx+1)..jp'last(2) loop
            declare
              tab : constant Table := Make_Table(jp(i,j),dim,nvr,tdx,deg);
              k : constant integer32 := j-1;
            begin
              tv.lead(k) := new Standard_Integer_Vectors.Vector'(tab.lead.all);
              tv.cffs(k) := tab.cffs;
              tv.mons(k) := tab.mons;
            end;
          end loop;
        end if;
        res(i) := new Table_Vector'(tv);
      end;
    end loop;
    return res;
  end Make_Table_Vector_Array;

end Double_Lseries_Polynomials;
