with text_io;                           use text_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Random_Vectors;
with Standard_Laurent_Series;
with Random_Laurent_Series;             use Random_Laurent_Series;

package body Standard_Lseries_Polynomials is

  procedure Write ( plead : in Standard_Integer_Vectors.Vector;
                    pcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                    pmons : in Standard_Integer_VecVecs.VecVec;
                    s : in string := "p" ) is
  begin
    for k in plead'range loop
      put(s & "("); put(k,1); put(") :"); put(pmons(k)); new_line;
      Standard_Laurent_Series.Write(plead(k),pcffs(k).all);
    end loop;
  end Write;

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
        ye := ye + xlead(i)*mon(i);
        for j in 1..mon(i) loop
          Standard_Laurent_Series.Multiply
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
                   yc : out Standard_Complex_Vectors.Vector ) is

    ze,ewrk : integer32;
    zc,cwrk : Standard_Complex_Vectors.Vector(0..deg);

  begin
    Eval(deg,plead(1),pcffs(1),pmons(1),xlead,xcffs,ye,yc);
    for i in 2..plead'last loop
      Eval(deg,plead(i),pcffs(i),pmons(i),xlead,xcffs,ze,zc);
      Standard_Laurent_Series.Add(deg,ye,ze,yc,zc,ewrk,cwrk);
      ye := ewrk;
      for k in 0..deg loop
        yc(k) := cwrk(k);
      end loop;
    end loop;
  end Eval;

  procedure Make_Series_Polynomial
              ( p : in Poly; dim,nvr,tdx,deg : in integer32 ) is

    nbr : constant integer32 := integer32(Number_of_Terms(p));
    plead : Standard_Integer_Vectors.Vector(1..nbr);
    cffs : Standard_Complex_VecVecs.VecVec(1..nbr);
    pcffs : constant Standard_Complex_VecVecs.Link_to_VecVec
          := new Standard_Complex_VecVecs.VecVec'(cffs);
    pmons : Standard_Integer_VecVecs.VecVec(1..nbr);
    cnt : integer32 := 0;

    procedure Visit_Term ( t : in Term; continue : out boolean ) is

      mon : Standard_Integer_Vectors.Vector(1..nvr);
      cff : Standard_Complex_Vectors.Vector(0..deg);

    begin
      cnt := cnt + 1; 
      if tdx = 0 then
        plead(cnt) := 0;
        for k in 1..dim loop
          mon(k) := integer32(t.dg(k));
        end loop;
      else
        for k in 1..(tdx-1) loop
          mon(k) := integer32(t.dg(k));
        end loop;
        plead(cnt) := integer32(t.dg(tdx));
        for k in (tdx+1)..dim loop
          mon(k-1) := integer32(t.dg(k));
        end loop;
      end if;
      pmons(cnt) := new Standard_Integer_Vectors.Vector'(mon);
      cff(0) := t.cf;
      for k in 1..deg loop
        cff(k) := Standard_Complex_Numbers.Create(0.0);
      end loop;
      pcffs(cnt) := new Standard_Complex_Vectors.Vector'(cff);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    put_line("The polynomial with Laurent series coefficients :");
    Write(plead,pcffs,pmons);
  end Make_Series_Polynomial;

end Standard_Lseries_Polynomials;
