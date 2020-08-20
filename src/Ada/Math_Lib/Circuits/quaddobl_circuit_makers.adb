with text_io;                             use text_io;
with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Quad_Double_Numbers;                 use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers_io;         use QuadDobl_Complex_Numbers_io;
with Standard_Random_Numbers;
with QuadDobl_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;         use Standard_Natural_Vectors_io;
with Standard_Integer_Vectors_io;         use Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
with Standard_Random_Vectors;
with QuadDobl_Random_Vectors;
with QuadDobl_Complex_Poly_Functions;
with Exponent_Indices;

package body QuadDobl_Circuit_Makers is

  function Random_Indices
             ( dim : integer32 )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..dim);
    cnt : integer32 := 0;
    rnd : integer32;

  begin
    loop
      for k in 1..dim loop
        rnd := Standard_Random_Numbers.Random(0,1);
        if rnd = 1 then
          cnt := cnt + 1;
          res(cnt) := k;
        end if;
      end loop;
      exit when (cnt > 1); -- at least two indices
      cnt := 0;            -- reset the counter
    end loop;
    return res(1..cnt);
  end Random_Indices;

  function Random_Indices
             ( dim,size : integer32 )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..size) := (1..size => 0);
    rnd,minval,minidx : integer32;
    redo : boolean; -- check for duplicate indices

  begin
    for k in 1..size loop
      redo := true;
      while redo loop
        rnd := Standard_Random_Numbers.Random(1,dim);
        redo := false;
        for i in 1..k-1 loop
          if res(i) = rnd
           then redo := true; exit;
          end if;
        end loop;
      end loop;
      res(k) := rnd;
    end loop;
   -- sort the indices
    for k in 1..size-1 loop
      minval := res(k); minidx := k;
      for i in k+1..size loop
        if minval > res(i)
         then minidx := i; minval := res(i);
        end if;
      end loop;
      if minidx /= k
       then res(minidx) := res(k); res(k) := minval;
      end if;
    end loop;
    return res;
  end Random_Indices;

  function Random_Complex_Circuit
             ( nbr,dim : integer32 )
             return QuadDobl_Complex_Circuits.Circuit is

    res : QuadDobl_Complex_Circuits.Circuit(nbr)
        := QuadDobl_Complex_Circuits.Allocate(nbr,dim);

  begin
    for k in 1..nbr loop
      res.xps(k) := new Standard_Integer_Vectors.Vector'(Random_Indices(dim));
    end loop;
    res.cff := QuadDobl_Random_Vectors.Random_Vector(1,nbr);
    res.cst := QuadDobl_Random_Numbers.Random1;
    res.pdg := Exponent_Indices.Polynomial_Degree(res.xps);
    return res;
  end Random_Complex_Circuit;

  function Random_Complex_Circuit
             ( nbr,dim,pwr : integer32 )
             return QuadDobl_Complex_Circuits.Circuit is

    res : QuadDobl_Complex_Circuits.Circuit(nbr)
        := QuadDobl_Complex_Circuits.Allocate(nbr,dim);
    xpk : Standard_Integer_Vectors.Vector(1..dim);

  begin
    for k in 1..nbr loop
      xpk := Standard_Random_Vectors.Random_Vector(1,dim,0,pwr);
      res.xps(k) := new Standard_Integer_Vectors.Vector'(xpk);
      res.idx(k) := Exponent_Indices.Exponent_Index(res.xps(k));
      res.fac(k) := Exponent_Indices.Factor_Index(res.xps(k));
    end loop;
    res.cff := QuadDobl_Random_Vectors.Random_Vector(1,nbr);
    res.cst := QuadDobl_Random_Numbers.Random1;
    res.pdg := Exponent_Indices.Polynomial_Degree(res.xps);
    return res;
  end Random_Complex_Circuit;

  function Random_Complex_Circuit
             ( nbr,dim,pwr : integer32 )
             return QuadDobl_Complex_Circuits.Link_to_Circuit is

    crc : constant QuadDobl_Complex_Circuits.Circuit(nbr)
        := Random_Complex_Circuit(nbr,dim,pwr);
    res : constant QuadDobl_Complex_Circuits.Link_to_Circuit
        := new QuadDobl_Complex_Circuits.Circuit'(crc);

  begin
    return res;
  end Random_Complex_Circuit;

  function Random_Complex_Circuits
             ( neq,nbr,dim,pwr : integer32 )
             return QuadDobl_Complex_Circuits.Circuits is

    res : QuadDobl_Complex_Circuits.Circuits(1..neq);

  begin
    for k in 1..neq loop
      res(k) := Random_Complex_Circuit(nbr,dim,pwr);
    end loop;
    return res;
  end Random_Complex_Circuits;

  function Random_Complex_System
             ( neq,nbr,dim,pwr : integer32 )
             return QuadDobl_Complex_Circuits.System is

    crc : constant QuadDobl_Complex_Circuits.Circuits(1..neq)
        := Random_Complex_Circuits(neq,nbr,dim,pwr);
    res : constant QuadDobl_Complex_Circuits.System(neq,dim)
        := QuadDobl_Complex_Circuits.Create(crc,dim);

  begin
    return res;
  end Random_Complex_System;

  function Random_Complex_System
             ( neq,nbr,dim,pwr : integer32 )
             return QuadDobl_Complex_Circuits.Link_to_System is

    sys : constant QuadDobl_Complex_Circuits.System(neq,dim)
        := Random_Complex_System(neq,nbr,dim,pwr);
    res : constant QuadDobl_Complex_Circuits.Link_to_System
        := new QuadDobl_Complex_Circuits.System'(sys);

  begin
    return res;
  end Random_Complex_System;

  function Make_Polynomial
             ( c : QuadDobl_Complex_Circuits.Circuit;
               index : boolean := false )
             return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Polynomials;

    res : Poly;
    t : Term;
    lnk : Standard_Integer_Vectors.Link_to_Vector;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..c.dim => 0);
    t.cf := c.cst;
    res := Create(t);
    for k in 1..c.nbr loop
      t.cf := c.cff(k);
      lnk := c.xps(k);
      if index then
        t.dg.all := (1..c.dim => 0);
        for i in lnk'range loop
          t.dg(lnk(i)) := 1;
        end loop;
      else
        for i in 1..c.dim loop
          t.dg(i) := natural32(lnk(i));
        end loop;
      end if;
      Add(res,t);
    end loop;
    Clear(t);
    return res;
  end Make_Polynomial;

  function Gradient ( p : QuadDobl_Complex_Polynomials.Poly;
                      x : QuadDobl_Complex_Vectors.Vector )
                    return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(x'range);
    dp : QuadDobl_Complex_Polynomials.Poly;

  begin
    for k in x'range loop
      dp := QuadDobl_Complex_Polynomials.Diff(p,k);
      res(k) := QuadDobl_Complex_Poly_Functions.Eval(dp,x);
      QuadDobl_Complex_Polynomials.Clear(dp);
    end loop;
    return res;
  end Gradient;

  function Hessian ( p : QuadDobl_Complex_Polynomials.Poly;
                     x : QuadDobl_Complex_Vectors.Vector )
                   return QuadDobl_Complex_Matrices.Matrix is

    res : QuadDobl_Complex_Matrices.Matrix(x'range,x'range);

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Functions;

    dp,dp2 : Poly;

  begin
    for i in res'range(1) loop
      dp := Diff(p,i);
      for j in res'range(2) loop
        dp2 := Diff(dp,j);
        res(i,j) := Eval(dp2,x);
        Clear(dp2);
      end loop;
      Clear(dp);
    end loop;
    return res;
  end Hessian;

  function Constant_Coefficient
             ( p : QuadDobl_Complex_Polynomials.Poly )
             return QuadDobl_Complex_Numbers.Complex_Number is

    use QuadDobl_Complex_Polynomials;

    res : QuadDobl_Complex_Numbers.Complex_Number;
    dim : constant integer32 := integer32(Number_of_Unknowns(p));
    zdg : constant Standard_Natural_Vectors.Vector(1..dim) := (1..dim => 0);
    deg : Degrees;

  begin
   -- put_line("inside constant coefficient ...");
    deg := new Standard_Natural_Vectors.Vector'(zdg);
   -- put("deg : "); put(deg.all); new_line;
    res := Coeff(p,deg);
    Clear(deg);
    return res;
  end Constant_Coefficient;

  function Is_NonZero ( c : QuadDobl_Complex_Numbers.Complex_Number )
                   return integer32 is

    zero : constant quad_double := create(0.0);

  begin
    if QuadDobl_Complex_Numbers.REAL_PART(c) /= zero then
      return 1;
    elsif QuadDobl_Complex_Numbers.IMAG_PART(c) /= zero then
      return 1;
    else
      return 0;
    end if;
  end Is_NonZero;

  function Make_Complex_Circuit
             ( p : QuadDobl_Complex_Polynomials.Poly;
               verbose : boolean := true )
             return QuadDobl_Complex_Circuits.Circuit is

    use QuadDobl_Complex_Polynomials;

    nbr : constant integer32 := integer32(Number_of_Terms(p));
    dim : constant integer32 := integer32(Number_of_Unknowns(p));
    cst : constant QuadDobl_Complex_Numbers.Complex_Number
        := Constant_Coefficient(p);
    isz : constant integer32 := Is_NonZero(cst);
    res : QuadDobl_Complex_Circuits.Circuit(nbr-isz)
        := QuadDobl_Complex_Circuits.Allocate(nbr-isz,dim);
    cnt : integer32 := 0;

    function Is_Zero ( d : in Degrees ) return boolean is
    
    -- DESCRIPTION :
    --   Returns true if all entries of d are zero.

    begin
      for k in d'range loop
        if d(k) /= 0
         then return false;
        end if;
      end loop;
      return true;
    end Is_Zero;

    procedure Visit_Term ( t : in Term; c : out boolean ) is

      xp : Standard_Integer_Vectors.Vector(t.dg'range);

    begin
      if not Is_Zero(t.dg) then
        cnt := cnt + 1;
        if verbose then
          put("visiting term "); put(cnt,1);
          put(" with exponents : "); put(t.dg.all); new_line;
        end if;
        res.cff(cnt) := t.cf;
        for i in xp'range loop
          xp(i) := integer32(t.dg(i));
        end loop;
        res.xps(cnt) := new Standard_Integer_Vectors.Vector'(xp);
        res.idx(cnt) := Exponent_Indices.Exponent_Index(res.xps(cnt));
        res.fac(cnt) := Exponent_Indices.Factor_Index(res.xps(cnt));
      end if;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    res.dim := dim;
    res.pdg := Degree(p);
    res.cst := cst;
    Visit_Terms(p);
    return res;
  end Make_Complex_Circuit;

  function Make_Complex_System
             ( p : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
               verbose : in boolean := true )
             return QuadDobl_Complex_Circuits.Link_to_System is

    use QuadDobl_Complex_Circuits;

    res : Link_to_System;
    c : Circuits(p'range);
    d : integer32;

  begin
    for k in c'range loop
      if verbose then
        put("The circuit for polynomial "); put(k,1);
        put_line(" :");
      end if;
      c(k) := new Circuit'(Make_Complex_Circuit(p(k)));
      if verbose then
        for i in 1..c(k).nbr loop
          put(c(k).cff(i)); put(c(k).xps(i)); new_line;
        end loop;
        put(c(k).cst); new_line;
        put("polynomial degree : "); put(c(k).pdg); new_line;
      end if;
    end loop;
    d := c(c'first).dim;
    res := new System'(Create(c,d));
    return res;
  end Make_Complex_System;

  procedure Write_Matrix ( A : in QuadDobl_Complex_Matrices.Matrix ) is
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("A["); put(i,1); put(","); put(j,1); put("] : ");
        put(A(i,j)); new_line;
      end loop;
    end loop;
  end Write_Matrix;

-- FROM CONVOLUTION CIRCUITS TO COMPLEX CIRCUITS :

  function Make_Complex_Circuit
             ( c : QuadDobl_Speelpenning_Convolutions.Circuit )
             return QuadDobl_Complex_Circuits.Circuit is

    res : QuadDobl_Complex_Circuits.Circuit(c.nbr)
        := QuadDobl_Complex_Circuits.Allocate(c.nbr,c.dim);
    lnk : QuadDobl_Complex_Vectors.Link_to_Vector;

    use QuadDobl_Complex_Vectors;

  begin
    Standard_Integer_VecVecs.Copy(c.xps,res.xps);
    Standard_Integer_VecVecs.Copy(c.idx,res.idx);
    Standard_Integer_VecVecs.Copy(c.fac,res.fac);
    res.pdg := Exponent_Indices.Polynomial_Degree(res.xps);
    for k in 1..c.nbr loop 
      lnk := c.cff(k);         -- coefficient power series
      res.cff(k) := lnk(0);    -- take leading coefficient
    end loop;
    if c.cst = null
     then res.cst := QuadDobl_Complex_Numbers.Create(integer(0));
     else res.cst := c.cst(0);
    end if;
    return res;
  end Make_Complex_Circuit;

  function Make_Complex_Circuit
             ( c : QuadDobl_Speelpenning_Convolutions.Link_to_Circuit )
             return QuadDobl_Complex_Circuits.Link_to_Circuit is

    res : QuadDobl_Complex_Circuits.Link_to_Circuit;

    use QuadDobl_Speelpenning_Convolutions;

  begin
    if c /= null then
      declare
        crc : constant QuadDobl_Complex_Circuits.Circuit(c.nbr)
            := Make_Complex_Circuit(c.all);
      begin
        res := new QuadDobl_Complex_Circuits.Circuit'(crc);
      end;
    end if;
    return res;
  end Make_Complex_Circuit;

  function Make_Complex_System
             ( s : QuadDobl_Speelpenning_Convolutions.System )
             return QuadDobl_Complex_Circuits.System is

    res : QuadDobl_Complex_Circuits.System(s.neq,s.dim);
    crc : QuadDobl_Complex_Circuits.Circuits(s.crc'range);

  begin
    for k in crc'range loop
      crc(k) := Make_Complex_Circuit(s.crc(k));
    end loop;
    res := QuadDobl_Complex_Circuits.Create(crc,s.dim);
    return res;
  end Make_Complex_System;

  function Make_Complex_System
             ( s : QuadDobl_Speelpenning_Convolutions.Link_to_System )
             return QuadDobl_Complex_Circuits.Link_to_System is

    res : QuadDobl_Complex_Circuits.Link_to_System;

    use QuadDobl_Speelpenning_Convolutions;

  begin
    if s /= null then
      declare
        cfs : constant QuadDobl_Complex_Circuits.System(s.neq,s.dim)
            := Make_Complex_System(s.all);
      begin
        res := new QuadDobl_Complex_Circuits.System'(cfs);
      end;
    end if;
    return res;
  end Make_Complex_System;

end QuadDobl_Circuit_Makers;
