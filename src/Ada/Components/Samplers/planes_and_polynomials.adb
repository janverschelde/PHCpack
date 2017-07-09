with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;

package body Planes_and_Polynomials is

-- AUXILIARIES :

  function Is_In ( v : Standard_Integer_Vectors.Vector; k : integer32 )
                 return boolean is

  -- DESCRIPTION :
  --   Returns true if there is an index i such that v(i) = k;
  --   returns false otherwise.

  begin
    for i in v'range loop
      if v(i) = k
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

  function Filter ( c : Standard_Complex_Numbers.Complex_Number;
                    tol : double_float )
                  return Standard_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   If the real or imaginary part of c is less than tol,
  --   then it will be equal to zero in the number on return.

    res : Standard_Complex_Numbers.Complex_Number;
    rec,imc : double_float;

  begin
    rec := REAL_PART(c);
    imc := IMAG_PART(c);
    if AbsVal(rec) > tol then
      if AbsVal(imc) > tol
       then res := c;
       else res := Create(rec);
      end if;
    else
      if AbsVal(imc) > tol
       then res := Create(0.0,imc);
       else res := Create(0.0);
      end if;
    end if;
    return res;
  end Filter;

  function Filter ( c : Multprec_Complex_Numbers.Complex_Number;
                    tol : double_float )
                  return Multprec_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   If the real or imaginary part of c is less than tol,
  --   then it will be equal to zero in the number on return.

    res : Multprec_Complex_Numbers.Complex_Number;
    rec,imc,absrec,absimc : Floating_Number;
    zero : Floating_Number;

  begin
    rec := REAL_PART(c); absrec := AbsVal(rec);
    imc := IMAG_PART(c); absimc := AbsVal(imc);
    if absrec > tol then
      if absimc > tol
       then Copy(c,res);
       else res := Create(rec);
      end if;
    else
      if absimc > tol then
        zero := Create(integer(0));
        res := Create(zero,imc);
        Clear(zero);
      else
        res := Create(integer(0));
      end if;
    end if;
    Clear(rec); Clear(absrec);
    Clear(imc); Clear(absimc);
    return res;
  end Filter;

  procedure Filter ( p : in out Standard_Complex_Polynomials.Poly;
                     tol : in double_float ) is

  -- DESCRIPTION :
  --   Removes terms in p with coefficients less than tol in magnitude.

    use Standard_Complex_Polynomials;
    fp : Standard_Complex_Polynomials.Poly := Null_Poly;

    procedure Filter_Term ( t : in Term; cont : out boolean ) is
    begin
      if AbsVal(t.cf) > tol
       then Add(fp,t);
      end if;
      cont := true;
    end Filter_Term;
    procedure Filter_Terms is new Visiting_Iterator(Filter_Term);

  begin
    Filter_Terms(p);
    Clear(p);
    p := fp;
  end Filter;

  procedure Filter ( p : in out Multprec_Complex_Polynomials.Poly;
                     tol : in double_float ) is

  -- DESCRIPTION :
  --   Removes terms in p with coefficients less than tol in magnitude.

    use Multprec_Complex_Polynomials;
    fp : Multprec_Complex_Polynomials.Poly := Null_Poly;

    procedure Filter_Term ( t : in Term; cont : out boolean ) is

      abscff : Floating_Number := AbsVal(t.cf);

    begin
      if abscff > tol
       then Add(fp,t);
      end if;
      Clear(abscff);
      cont := true;
    end Filter_Term;
    procedure Filter_Terms is new Visiting_Iterator(Filter_Term);

  begin
    Filter_Terms(p);
    Clear(p);
    p := fp;
  end Filter;

-- CONVERTORS :

  function Hyperplane ( cff : Standard_Complex_Vectors.Vector )
                      return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;
    res : Poly := Null_Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..cff'last => 0);
    if cff(0) /= Create(0.0) then
      t.cf := cff(0);
      Add(res,t);
    end if;
    for i in 1..cff'last loop
      if cff(i) /= Create(0.0) then
        t.dg(i) := 1;
        t.cf := cff(i);
        Add(res,t);
        t.dg(i) := 0;
      end if;
    end loop;
    Clear(t);
    return res;
  end Hyperplane;

  function Hyperplane ( cff : Standard_Complex_Vectors.Vector )
                      return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Laurentials;
    res : Poly := Null_Poly;
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..cff'last => 0);
    if cff(0) /= Create(0.0) then
      t.cf := cff(0);
      Add(res,t);
    end if;
    for i in 1..cff'last loop
      if cff(i) /= Create(0.0) then
        t.dg(i) := 1;
        t.cf := cff(i);
        Add(res,t);
        t.dg(i) := 0;
      end if;
    end loop;
    Clear(t);
    return res;
  end Hyperplane;

  function Hyperplane ( cff : Standard_Complex_Vectors.Vector;
                        tol : double_float )
                      return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;
    res : Poly := Null_Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..cff'last => 0);
    if AbsVal(cff(0)) > tol then
      t.cf := Filter(cff(0),tol);
      Add(res,t);
    end if;
    for i in 1..cff'last loop
      if AbsVal(cff(i)) > tol then
        t.dg(i) := 1;
        t.cf := Filter(cff(i),tol);
        Add(res,t);
        t.dg(i) := 0;
      end if;
    end loop;
    Clear(t);
    return res;
  end Hyperplane;

  function Hyperplane ( cff : DoblDobl_Complex_Vectors.Vector )
                      return DoblDobl_Complex_Polynomials.Poly is

    use DoblDobl_Complex_Numbers,DoblDobl_Complex_Polynomials;
    ddzero : constant double_double := Double_Double_Numbers.Create(0.0);
    res : Poly := Null_Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..cff'last => 0);
    if cff(0) /= Create(ddzero) then
      t.cf := cff(0);
      Add(res,t);
    end if;
    for i in 1..cff'last loop
      if cff(i) /= Create(ddzero) then
        t.dg(i) := 1;
        t.cf := cff(i);
        Add(res,t);
        t.dg(i) := 0;
      end if;
    end loop;
    Clear(t);
    return res;
  end Hyperplane;

  function Hyperplane ( cff : DoblDobl_Complex_Vectors.Vector )
                      return DoblDobl_Complex_Laurentials.Poly is

    use DoblDobl_Complex_Numbers,DoblDobl_Complex_Laurentials;
    ddzero : constant double_double := Double_Double_Numbers.Create(0.0);
    res : Poly := Null_Poly;
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..cff'last => 0);
    if cff(0) /= Create(ddzero) then
      t.cf := cff(0);
      Add(res,t);
    end if;
    for i in 1..cff'last loop
      if cff(i) /= Create(ddzero) then
        t.dg(i) := 1;
        t.cf := cff(i);
        Add(res,t);
        t.dg(i) := 0;
      end if;
    end loop;
    Clear(t);
    return res;
  end Hyperplane;

  function Hyperplane ( cff : QuadDobl_Complex_Vectors.Vector )
                      return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Numbers,QuadDobl_Complex_Polynomials;
    ddzero : constant quad_double
           := Quad_Double_Numbers.Create(0.0);
    res : Poly := Null_Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..cff'last => 0);
    if cff(0) /= Create(ddzero) then
      t.cf := cff(0);
      Add(res,t);
    end if;
    for i in 1..cff'last loop
      if cff(i) /= Create(ddzero) then
        t.dg(i) := 1;
        t.cf := cff(i);
        Add(res,t);
        t.dg(i) := 0;
      end if;
    end loop;
    Clear(t);
    return res;
  end Hyperplane;

  function Hyperplane ( cff : QuadDobl_Complex_Vectors.Vector )
                      return QuadDobl_Complex_Laurentials.Poly is

    use QuadDobl_Complex_Numbers,QuadDobl_Complex_Laurentials;
    ddzero : constant quad_double
           := Quad_Double_Numbers.Create(0.0);
    res : Poly := Null_Poly;
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..cff'last => 0);
    if cff(0) /= Create(ddzero) then
      t.cf := cff(0);
      Add(res,t);
    end if;
    for i in 1..cff'last loop
      if cff(i) /= Create(ddzero) then
        t.dg(i) := 1;
        t.cf := cff(i);
        Add(res,t);
        t.dg(i) := 0;
      end if;
    end loop;
    Clear(t);
    return res;
  end Hyperplane;

  function Hyperplane ( cff : Multprec_Complex_Vectors.Vector )
                      return Multprec_Complex_Polynomials.Poly is

    use Multprec_Complex_Polynomials;
    res : Poly := Null_Poly;
    t : Term;
    zero : Multprec_Complex_Numbers.Complex_Number := Create(integer(0));

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..cff'last => 0);
    if not Equal(cff(0),zero) then
      Copy(cff(0),t.cf);
      Add(res,t);
    end if;
    for i in 1..cff'last loop
      if not Equal(cff(i),zero) then
        t.dg(i) := 1;
        Copy(cff(i),t.cf);
        Add(res,t);
        t.dg(i) := 0;
      end if;
    end loop;
    Clear(t);
    Clear(zero);
    return res;
  end Hyperplane;

  function Hyperplane ( cff : Multprec_Complex_Vectors.Vector )
                      return Multprec_Complex_Laurentials.Poly is

    use Multprec_Complex_Laurentials;
    res : Poly := Null_Poly;
    t : Term;
    zero : Multprec_Complex_Numbers.Complex_Number := Create(integer(0));

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..cff'last => 0);
    if not Equal(cff(0),zero) then
      Copy(cff(0),t.cf);
      Add(res,t);
    end if;
    for i in 1..cff'last loop
      if not Equal(cff(i),zero) then
        t.dg(i) := 1;
        Copy(cff(i),t.cf);
        Add(res,t);
        t.dg(i) := 0;
      end if;
    end loop;
    Clear(t);
    Clear(zero);
    return res;
  end Hyperplane;

  function Hyperplane ( cff : Multprec_Complex_Vectors.Vector;
                        tol : double_float )
                      return Multprec_Complex_Polynomials.Poly is

    use Multprec_Complex_Polynomials;
    res : Poly := Null_Poly;
    t : Term;
    abscff : Floating_Number;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..cff'last => 0);
    abscff := AbsVal(cff(0));
    if abscff > tol then
      t.cf := Filter(cff(0),tol);
      Add(res,t);
    end if;
    Clear(abscff);
    for i in 1..cff'last loop
      abscff := AbsVal(cff(i));
      if abscff > tol then
        t.dg(i) := 1;
        t.cf := Filter(cff(i),tol);
        Add(res,t);
        t.dg(i) := 0;
      end if;
      Clear(abscff);
    end loop;
    Clear(t);
    return res;
  end Hyperplane;

  function Polynomial ( p : Standard_Complex_Polynomials.Poly )
                      return Standard_Complex_Vectors.Vector is

    use Standard_Complex_Polynomials;
    n : constant integer32 := integer32(Number_of_Unknowns(p));
    res : Standard_Complex_Vectors.Vector(0..n) := (0..n => Create(0.0));

    procedure Scan_Term ( t : in Term; continue : out boolean ) is

      found : boolean := false;

    begin
      for i in t.dg'range loop
        if t.dg(i) = 1 then
          res(i) := t.cf;
          found := true;
        end if;
        exit when found;
      end loop;
      if not found
       then res(0) := t.cf;
      end if;     
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Polynomial;

  function Polynomial ( p : Standard_Complex_Laurentials.Poly )
                      return Standard_Complex_Vectors.Vector is

    use Standard_Complex_Laurentials;
    n : constant integer32 := integer32(Number_of_Unknowns(p));
    res : Standard_Complex_Vectors.Vector(0..n) := (0..n => Create(0.0));

    procedure Scan_Term ( t : in Term; continue : out boolean ) is

      found : boolean := false;

    begin
      for i in t.dg'range loop
        if t.dg(i) = 1 then
          res(i) := t.cf;
          found := true;
        end if;
        exit when found;
      end loop;
      if not found
       then res(0) := t.cf;
      end if;     
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Polynomial;

  function Polynomial ( p : DoblDobl_Complex_Polynomials.Poly )
                      return DoblDobl_Complex_Vectors.Vector is

    use DoblDobl_Complex_Polynomials;
    n : constant integer32 := integer32(Number_of_Unknowns(p));
    ddzero : constant double_double := Double_Double_Numbers.Create(0.0);
    res : DoblDobl_Complex_Vectors.Vector(0..n)
        := (0..n => DoblDobl_Complex_Numbers.Create(ddzero));

    procedure Scan_Term ( t : in Term; continue : out boolean ) is

      found : boolean := false;

    begin
      for i in t.dg'range loop
        if t.dg(i) = 1 then
          res(i) := t.cf;
          found := true;
        end if;
        exit when found;
      end loop;
      if not found
       then res(0) := t.cf;
      end if;     
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Polynomial;

  function Polynomial ( p : DoblDobl_Complex_Laurentials.Poly )
                      return DoblDobl_Complex_Vectors.Vector is

    use DoblDobl_Complex_Laurentials;
    n : constant integer32 := integer32(Number_of_Unknowns(p));
    ddzero : constant double_double := Double_Double_Numbers.Create(0.0);
    res : DoblDobl_Complex_Vectors.Vector(0..n)
        := (0..n => DoblDobl_Complex_Numbers.Create(ddzero));

    procedure Scan_Term ( t : in Term; continue : out boolean ) is

      found : boolean := false;

    begin
      for i in t.dg'range loop
        if t.dg(i) = 1 then
          res(i) := t.cf;
          found := true;
        end if;
        exit when found;
      end loop;
      if not found
       then res(0) := t.cf;
      end if;     
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Polynomial;

  function Polynomial ( p : QuadDobl_Complex_Polynomials.Poly )
                      return QuadDobl_Complex_Vectors.Vector is

    use QuadDobl_Complex_Polynomials;
    n : constant integer32 := integer32(Number_of_Unknowns(p));
    qdzero : constant quad_double := Quad_Double_Numbers.Create(0.0);
    res : QuadDobl_Complex_Vectors.Vector(0..n)
        := (0..n => QuadDobl_Complex_Numbers.Create(qdzero));

    procedure Scan_Term ( t : in Term; continue : out boolean ) is

      found : boolean := false;

    begin
      for i in t.dg'range loop
        if t.dg(i) = 1 then
          res(i) := t.cf;
          found := true;
        end if;
        exit when found;
      end loop;
      if not found
       then res(0) := t.cf;
      end if;     
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Polynomial;

  function Polynomial ( p : QuadDobl_Complex_Laurentials.Poly )
                      return QuadDobl_Complex_Vectors.Vector is

    use QuadDobl_Complex_Laurentials;
    n : constant integer32 := integer32(Number_of_Unknowns(p));
    qdzero : constant quad_double := Quad_Double_Numbers.Create(0.0);
    res : QuadDobl_Complex_Vectors.Vector(0..n)
        := (0..n => QuadDobl_Complex_Numbers.Create(qdzero));

    procedure Scan_Term ( t : in Term; continue : out boolean ) is

      found : boolean := false;

    begin
      for i in t.dg'range loop
        if t.dg(i) = 1 then
          res(i) := t.cf;
          found := true;
        end if;
        exit when found;
      end loop;
      if not found
       then res(0) := t.cf;
      end if;     
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Polynomial;

  function Polynomial ( p : Multprec_Complex_Polynomials.Poly )
                      return Multprec_Complex_Vectors.Vector is

    use Multprec_Complex_Polynomials;
    n : constant integer32 := integer32(Number_of_Unknowns(p));
    res : Multprec_Complex_Vectors.Vector(0..n) := (0..n => Create(integer(0)));

    procedure Scan_Term ( t : in Term; continue : out boolean ) is

      found : boolean := false;

    begin
      for i in t.dg'range loop
        if t.dg(i) = 1 then
          Copy(t.cf,res(i));
          found := true;
        end if;
        exit when found;
      end loop;
      if not found
       then Copy(t.cf,res(0));
      end if;     
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Polynomial;

  function Polynomial ( p : Multprec_Complex_Laurentials.Poly )
                      return Multprec_Complex_Vectors.Vector is

    use Multprec_Complex_Laurentials;
    n : constant integer32 := integer32(Number_of_Unknowns(p));
    res : Multprec_Complex_Vectors.Vector(0..n) := (0..n => Create(integer(0)));

    procedure Scan_Term ( t : in Term; continue : out boolean ) is

      found : boolean := false;

    begin
      for i in t.dg'range loop
        if t.dg(i) = 1 then
          Copy(t.cf,res(i));
          found := true;
        end if;
        exit when found;
      end loop;
      if not found
       then Copy(t.cf,res(0));
      end if;     
      continue := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Polynomial;

-- ELIMINATORS FOR POLYNOMIALS :

  function Remove_Variable
             ( p : Standard_Complex_Polynomials.Poly; k : integer32 )
             return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;
    res : Standard_Complex_Polynomials.Poly := Null_Poly;

    procedure Remove_Variable_in_Term ( t : in Term; cont : out boolean ) is

      rt : Term;

    begin
      rt.cf := t.cf;
      rt.dg := new Standard_Natural_Vectors.Vector(t.dg'first..t.dg'last-1);
      for i in rt.dg'range loop
        if i < k
         then rt.dg(i) := t.dg(i);
         else rt.dg(i) := t.dg(i+1);
        end if;
      end loop;
      Add(res,rt);
      Clear(rt);
      cont := true;
    end Remove_Variable_in_Term;
    procedure Remove_Variable_in_Terms is
      new Visiting_Iterator(Remove_Variable_in_Term);

  begin
    Remove_Variable_in_Terms(p);
    return res;
  end Remove_Variable;

  function Remove_Variable
             ( p : Multprec_Complex_Polynomials.Poly; k : integer32 )
             return Multprec_Complex_Polynomials.Poly is

    use Multprec_Complex_Polynomials;
    res : Multprec_Complex_Polynomials.Poly := Null_Poly;

    procedure Remove_Variable_in_Term ( t : in Term; cont : out boolean ) is

      rt : Term;

    begin
      Copy(t.cf,rt.cf);
      rt.dg := new Standard_Natural_Vectors.Vector(t.dg'first..t.dg'last-1);
      for i in rt.dg'range loop
        if i < k
         then rt.dg(i) := t.dg(i);
         else rt.dg(i) := t.dg(i+1);
        end if;
      end loop;
      Add(res,rt);
      Clear(rt);
      cont := true;
    end Remove_Variable_in_Term;
    procedure Remove_Variable_in_Terms is
      new Visiting_Iterator(Remove_Variable_in_Term);

  begin
    Remove_Variable_in_Terms(p);
    return res;
  end Remove_Variable;

  function Substituting_Polynomial
             ( nvar : integer32; pivots : Standard_Integer_Vectors.Vector;
               hypcff : Standard_Complex_Vectors.Vector; tol : double_float )
             return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;
    res : Standard_Complex_Polynomials.Poly := Null_Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..nvar => 0);
    if AbsVal(hypcff(0)) > tol then
      t.cf := -Filter(hypcff(0),tol);
      Add(res,t);
    end if;
    for i in pivots'range loop
      if AbsVal(hypcff(pivots(i))) > tol then
        t.cf := -Filter(hypcff(pivots(i)),tol);
        t.dg(i) := 1;
        Add(res,t);
        t.dg(i) := 0;
      end if;
    end loop;
    Clear(t);
    return res;
  end Substituting_Polynomial;

  function Substituting_Polynomial
             ( nvar : integer32; pivots : Standard_Integer_Vectors.Vector;
               hypcff : Multprec_Complex_Vectors.Vector; tol : double_float )
             return Multprec_Complex_Polynomials.Poly is

    use Multprec_Complex_Polynomials;
    res : Multprec_Complex_Polynomials.Poly := Null_Poly;
    abscff : Floating_Number;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..nvar => 0);
    abscff := AbsVal(hypcff(0));
    if abscff > tol then
      t.cf := Filter(hypcff(0),tol);
      --Copy(hypcff(0),t.cf);
      Min(t.cf);
      Add(res,t);
    end if;
    Clear(abscff);
    for i in pivots'range loop
      abscff := AbsVal(hypcff(pivots(i)));
      if abscff > tol then
        t.cf := Filter(hypcff(pivots(i)),tol);
        --Copy(hypcff(pivots(i)),t.cf);
        Min(t.cf);
        t.dg(i) := 1;
        Add(res,t);
        t.dg(i) := 0;
      end if;
      Clear(abscff);
    end loop;
    Clear(t);
    return res;
  end Substituting_Polynomial;

  function Substituting_Polynomials
             ( nequ,nvar : integer32; pivots : Standard_Integer_Vectors.Vector;
               hypcff : Standard_Complex_VecVecs.VecVec; tol : double_float )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(1..nequ);

  begin
    for i in res'range loop
      res(i) := Substituting_Polynomial(nvar,pivots,hypcff(i).all,tol);
    end loop;
    return res;
  end Substituting_Polynomials;

  function Substituting_Polynomials
             ( nequ,nvar : integer32; pivots : Standard_Integer_Vectors.Vector;
               hypcff : Multprec_Complex_VecVecs.VecVec; tol : double_float )
             return Multprec_Complex_Poly_Systems.Poly_Sys is

    res : Multprec_Complex_Poly_Systems.Poly_Sys(1..nequ);

  begin
    for i in res'range loop
      res(i) := Substituting_Polynomial(nvar,pivots,hypcff(i).all,tol);
    end loop;
    return res;
  end Substituting_Polynomials;

  function Substitute
             ( p : Standard_Complex_Polynomials.Poly; nvar,level : integer32;
               pivots : Standard_Integer_Vectors.Vector;
               subpols : Standard_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Polynomials.Poly is

  -- DESCRIPTION : 
  --   Given a polynomial for each nonpivot variable, the polynomial
  --   on return will have these nonpivot variables eliminated.

    use Standard_Complex_Polynomials;
    res : Standard_Complex_Polynomials.Poly := Null_Poly;

    procedure Substitute_Term ( t : in Term; cont : out boolean ) is

      st : Term;
      accfac,fac : Poly;
      ind : integer32 := 0;
      empty : boolean := false;

    begin
      for i in t.dg'first..t.dg'last-level loop       -- substitute nonpivots
        if not Is_In(pivots,i) then
          ind := ind+1;
          if t.dg(i) /= 0 then
            if subpols(ind) = Null_Poly then
              empty := true;
              Clear(accfac);
            else
              Copy(subpols(ind),fac);
              for j in 2..t.dg(i) loop
                Mul(fac,subpols(ind));
              end loop;
              if accfac = Null_Poly
               then Copy(fac,accfac);
               else Mul(accfac,fac);
              end if;
              Clear(fac);
            end if;
          end if;
        end if;
        exit when empty;
      end loop;
      if not empty then
        st.cf := t.cf;
        st.dg := new Standard_Natural_Vectors.Vector'(1..nvar => 0);
        for i in pivots'range loop                         -- copy pivots
          st.dg(i) := t.dg(pivots(i));
        end loop;
        for i in 1..level loop                     -- copy last variables
          st.dg(pivots'last+i) := t.dg(t.dg'last-level+i);
        end loop;
        if accfac /= Null_Poly then
          Mul(accfac,st);
          Add(res,accfac);
          Clear(accfac);
        else
          Add(res,st);
        end if;
        Clear(st);
      end if;
      cont := true;
    end Substitute_Term;
    procedure Substitute_Terms is new Visiting_Iterator(Substitute_Term);

  begin
    Substitute_Terms(p);
    return res;
  end Substitute;

  function Substitute
             ( p : Multprec_Complex_Polynomials.Poly; nvar,level : integer32;
               pivots : Standard_Integer_Vectors.Vector;
               subpols : Multprec_Complex_Poly_Systems.Poly_Sys )
             return Multprec_Complex_Polynomials.Poly is

  -- DESCRIPTION : 
  --   Given a polynomial for each nonpivot variable, the polynomial
  --   on return will have these nonpivot variables eliminated.

    use Multprec_Complex_Polynomials;
    res : Multprec_Complex_Polynomials.Poly := Null_Poly;

    procedure Substitute_Term ( t : in Term; cont : out boolean ) is

      st : Term;
      accfac,fac : Poly;
      ind : integer32 := 0;
      empty : boolean := false;

    begin
      for i in t.dg'first..t.dg'last-level loop       -- substitute nonpivots
        if not Is_In(pivots,i) then
          ind := ind+1;
          if t.dg(i) /= 0 then
            if subpols(ind) = Null_Poly then
              empty := true;
              Clear(accfac);
            else
              Copy(subpols(ind),fac);
              for j in 2..t.dg(i) loop
                Mul(fac,subpols(ind));
              end loop;
              if accfac = Null_Poly
               then Copy(fac,accfac);
               else Mul(accfac,fac);
              end if;
              Clear(fac);
            end if;
          end if;
        end if;
      end loop;
      if not empty then
        Copy(t.cf,st.cf);
        st.dg := new Standard_Natural_Vectors.Vector'(1..nvar => 0);
        for i in pivots'range loop                         -- copy pivots
          st.dg(i) := t.dg(pivots(i));
        end loop;
        for i in 1..level loop                     -- copy last variables
          st.dg(pivots'last+i) := t.dg(t.dg'last-level+i);
        end loop;
        if accfac /= Null_Poly then
          Mul(accfac,st);
          Add(res,accfac);
          Clear(accfac);
        else
          Add(res,st);
        end if;
        Clear(st);
      end if;
      cont := true;
    end Substitute_Term;
    procedure Substitute_Terms is new Visiting_Iterator(Substitute_Term);

  begin
    Substitute_Terms(p);
    return res;
  end Substitute;

  function Restrict_to_Linear_Space
             ( p : Standard_Complex_Polynomials.Poly; level : integer32;
               pivots : Standard_Integer_Vectors.Vector;
               subpols : Standard_Complex_Poly_Systems.Poly_Sys;
               tol : double_float ) return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly;
    use Standard_Complex_Polynomials;
    dim : constant integer32 := pivots'length + level;

  begin
    res := Substitute(p,dim,level,pivots,subpols);
    Filter(res,tol);
    return res;
  end Restrict_to_Linear_Space;

  function Restrict_to_Linear_Space
             ( p : Multprec_Complex_Polynomials.Poly; level : integer32;
               pivots : Standard_Integer_Vectors.Vector;
               subpols : Multprec_Complex_Poly_Systems.Poly_Sys;
               tol : double_float ) return Multprec_Complex_Polynomials.Poly is

    res : Multprec_Complex_Polynomials.Poly;
    use Multprec_Complex_Polynomials;
    dim : constant integer32 := pivots'length + level;

  begin
    res := Substitute(p,dim,level,pivots,subpols);
    Filter(res,tol);
    return res;
  end Restrict_to_Linear_Space;

  function Restrict_to_Linear_Space
             ( p : Standard_Complex_Poly_Systems.Poly_Sys; level : integer32;
               pivots : Standard_Integer_Vectors.Vector;
               hypcff : Standard_Complex_VecVecs.VecVec; tol : double_float )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    dim : constant integer32 := pivots'length + level;
    neq : constant integer32 := hypcff'length;
    subpols : Standard_Complex_Poly_Systems.Poly_Sys(1..neq);

  begin
    subpols := Substituting_Polynomials(neq,dim,pivots,hypcff,tol);
    for i in p'range loop
      res(i) := Substitute(p(i),dim,level,pivots,subpols);
      Filter(res(i),tol);
    end loop;
    return res;
  end Restrict_to_Linear_Space;

  function Restrict_to_Linear_Space
             ( p : Multprec_Complex_Poly_Systems.Poly_Sys; level : integer32;
               pivots : Standard_Integer_Vectors.Vector;
               hypcff : Multprec_Complex_VecVecs.VecVec; tol : double_float )
             return Multprec_Complex_Poly_Systems.Poly_Sys is

    res : Multprec_Complex_Poly_Systems.Poly_Sys(p'range);
    dim : constant integer32 := pivots'length + level;
    neq : constant integer32 := hypcff'length;
    subpols : Multprec_Complex_Poly_Systems.Poly_Sys(1..neq);

  begin
    subpols := Substituting_Polynomials(neq,dim,pivots,hypcff,tol);
    for i in p'range loop
      res(i) := Substitute(p(i),dim,level,pivots,subpols);
      Filter(res(i),tol);
    end loop;
    return res;
  end Restrict_to_Linear_Space;

-- ELIMINATORS FOR SOLUTIONS :

  function Remove_Variables
             ( v : Standard_Complex_Vectors.Vector; level,newdim : integer32;
               pivots : Standard_Integer_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is
    
    res : Standard_Complex_Vectors.Vector(1..newdim);

  begin
    for i in pivots'range loop
      res(i) := v(pivots(i));
    end loop;
    for i in 1..level loop
      res(pivots'last+i) := v(v'last+i-level);
    end loop;
    return res;
  end Remove_Variables;

  function Remove_Variables
             ( v : Multprec_Complex_Vectors.Vector; level,newdim : integer32;
               pivots : Standard_Integer_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector is
    
    res : Multprec_Complex_Vectors.Vector(1..newdim);

  begin
    for i in pivots'range loop
      Copy(v(pivots(i)),res(i));
    end loop;
    for i in 1..level loop
      Copy(v(v'last+i-level),res(pivots'last+i));
    end loop;
    return res;
  end Remove_Variables;

  function Remove_Variables
             ( s : Standard_Complex_Solutions.Solution;
               level : integer32; pivots : Standard_Integer_Vectors.Vector )
             return Standard_Complex_Solutions.Solution is

    nn : constant integer32 := pivots'length+level;
    res : Standard_Complex_Solutions.Solution(nn);

  begin
    res.t := s.t;
    res.m := s.m;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    res.v := Remove_Variables(s.v,level,nn,pivots);
    return res;
  end Remove_Variables;

  function Remove_Variables
             ( s : Multprec_Complex_Solutions.Solution;
               level : integer32; pivots : Standard_Integer_Vectors.Vector )
             return Multprec_Complex_Solutions.Solution is

    nn : constant integer32 := pivots'length+level;
    res : Multprec_Complex_Solutions.Solution(nn);

  begin
    res.t := s.t;
    res.m := s.m;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    res.v := Remove_Variables(s.v,level,nn,pivots);
    return res;
  end Remove_Variables;

  function Remove_Variables
             ( sols : Standard_Complex_Solutions.Solution_List;
               level : integer32; pivots : Standard_Integer_Vectors.Vector )
             return Standard_Complex_Solutions.Solution_List is

    res,res_last : Standard_Complex_Solutions.Solution_List;
    use Standard_Complex_Solutions;
    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Remove_Variables(Head_Of(tmp).all,level,pivots));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Remove_Variables;

  function Remove_Variables
             ( sols : Multprec_Complex_Solutions.Solution_List;
               level : integer32; pivots : Standard_Integer_Vectors.Vector )
             return Multprec_Complex_Solutions.Solution_List is

    res,res_last : Multprec_Complex_Solutions.Solution_List;
    use Multprec_Complex_Solutions;
    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Remove_Variables(Head_Of(tmp).all,level,pivots));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Remove_Variables;

  function Restrict_Solution
             ( sols : Standard_Complex_Solutions.Solution_List;
               ind,level : integer32;
               pivots : Standard_Integer_Vectors.Vector ) 
             return Standard_Complex_Solutions.Solution_List is

    use Standard_Complex_Solutions;
    res,res_last,tmp : Solution_List;
    ls : Link_to_Solution;

  begin
    tmp := sols;
    for i in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      if integer32(i) = ind
       then Append(res,res_last,Remove_Variables(ls.all,level,pivots));
       else Append(res,res_last,ls);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Restrict_Solution;

end Planes_and_Polynomials;
