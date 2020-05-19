with text_io;                             use text_io;
with Communications_with_User;            use Communications_with_User;
with Timing_Package;                      use Timing_Package;
with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;        use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;         use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;         use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;         use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;         use Standard_Complex_VecVecs_io;
with Standard_Complex_Matrices;
with Standard_Random_Vectors;
with Standard_Vector_Splitters;           use Standard_Vector_Splitters;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;    use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with Exponent_Indices;
with Standard_Complex_Circuits;
with Standard_Coefficient_Circuits;
with Evaluation_Differentiation_Errors;

procedure ts_perfade is

-- DESCRIPTION :
--   Tests better performing algorithmic differentiation and evaluation.

  procedure Test_Forward ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random vector of dimension dim
  --   and tests the computation of the forward products.

    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    cf : constant Standard_Complex_Vectors.Vector(1..dim-1)
       := Standard_Complex_Vectors.Vector'(1..dim-1 => zero);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    f : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cf);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    fr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(f);
    fi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(f);
    v : Standard_Complex_Vectors.Link_to_Vector;
    err : double_float;

  begin
    Standard_Complex_Circuits.Forward(x,f);
    put_line("the result : "); put_line(f);
    Standard_Coefficient_Circuits.Forward(xr,xi,fr,fi);
    v := Make_Complex(fr,fi);
    put_line("recomputed : "); put_line(v);
    err := Evaluation_Differentiation_Errors.Difference(f,v);
    put("The error :"); put(err,3); new_line;
  end Test_Forward;

  procedure Test_Forward_Backward ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random vector of dimension dim
  --   and tests the computation of the forward/backward products.

    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    cf : constant Standard_Complex_Vectors.Vector(1..dim-1)
       := Standard_Complex_Vectors.Vector'(1..dim-1 => zero);
    cf2 : constant Standard_Complex_Vectors.Vector(1..dim-1)
        := Standard_Complex_Vectors.Vector'(1..dim-1 => zero);
    cb : constant Standard_Complex_Vectors.Vector(1..dim-2)
        := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    cb2 : constant Standard_Complex_Vectors.Vector(1..dim-2)
        := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    f : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cf);
    f2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cf2);
    b : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cb);
    b2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cb2);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    fr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(f);
    fi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(f);
    br : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(b);
    bi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(b);
    fr2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(f2);
    fi2 : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(f2);
    br2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(b2);
    bi2 : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(b2);
    v,w,v2,w2 : Standard_Complex_Vectors.Link_to_Vector;
    err,sumerr : double_float;

  begin
    Standard_Complex_Circuits.Forward_Backward(x,f,b);
    Standard_Complex_Circuits.Fused_Forward_Backward(x,f2,b2);
    Standard_Coefficient_Circuits.Forward_Backward(xr,xi,fr,fi,br,bi);
    Standard_Coefficient_Circuits.Fused_Forward_Backward(xr,xi,fr2,fi2,br2,bi2);
    v := Make_Complex(fr,fi); v2 := Make_Complex(fr2,fi2);
    w := Make_Complex(br,bi); w2 := Make_Complex(br2,bi2);
    put_line("the forward products : "); put_line(f);
    put_line("the forward products with loop fusion : "); put_line(f2);
    err := Evaluation_Differentiation_Errors.Difference(f,f2);
    put("The error :"); put(err,3); new_line;
    sumerr := err;
    put_line("forward products recomputed : "); put_line(v);
    err := Evaluation_Differentiation_Errors.Difference(f,v);
    put("The error :"); put(err,3); new_line;
    sumerr := sumerr + err;
    put_line("forward products recomputed with loop fusion : "); put_line(v2);
    err := Evaluation_Differentiation_Errors.Difference(f,v2);
    put("The error :"); put(err,3); new_line;
    sumerr := sumerr + err;
    put_line("the backward products : "); put_line(b);
    put_line("the backward products with loop fusion : "); put_line(b2);
    err := Evaluation_Differentiation_Errors.Difference(b,b2);
    put("The error :"); put(err,3); new_line;
    sumerr := sumerr + err;
    put_line("backward products recomputed : "); put_line(w);
    err := Evaluation_Differentiation_Errors.Difference(b,w);
    put("The error :"); put(err,3); new_line;
    sumerr := sumerr + err;
    put_line("backward products recomputed with loop fusion : "); put_line(w2);
    err := Evaluation_Differentiation_Errors.Difference(b,w2);
    put("The error :"); put(err,3); new_line;
    sumerr := sumerr + err;
    put("The sum of all errors :"); put(sumerr,3); new_line;
  end Test_Forward_Backward;

  procedure Test_Forward_Backward_Cross ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random vector of dimension dim
  --   and tests the computation of the forward/backward/cross products.

    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    cf : constant Standard_Complex_Vectors.Vector(1..dim-1)
       := Standard_Complex_Vectors.Vector'(1..dim-1 => zero);
    cf2 : constant Standard_Complex_Vectors.Vector(1..dim-1)
        := Standard_Complex_Vectors.Vector'(1..dim-1 => zero);
    cb : constant Standard_Complex_Vectors.Vector(1..dim-2)
        := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    cb2 : constant Standard_Complex_Vectors.Vector(1..dim-2)
        := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    cc : constant Standard_Complex_Vectors.Vector(1..dim-2)
        := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    cc2 : constant Standard_Complex_Vectors.Vector(1..dim-2)
        := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    f : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cf);
    f2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cf2);
    b : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cb);
    b2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cb2);
    c : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cc);
    c2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cc2);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    fr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(f);
    fi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(f);
    br : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(b);
    bi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(b);
    cr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(c);
    ci : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(c);
    fr2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(f2);
    fi2 : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(f2);
    br2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(b2);
    bi2 : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(b2);
    cr2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(c2);
    ci2 : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(c2);
    u,v,w,u2,v2,w2 : Standard_Complex_Vectors.Link_to_Vector;
    err,sumerr : double_float;

  begin
    Standard_Complex_Circuits.Forward_Backward_Cross(x,f,b,c);
    Standard_Complex_Circuits.Fused_Forward_Backward_Cross(x,f2,b2,c2);
    Standard_Coefficient_Circuits.Forward_Backward_Cross
      (xr,xi,fr,fi,br,bi,cr,ci);
    Standard_Coefficient_Circuits.Fused_Forward_Backward_Cross
      (xr,xi,fr2,fi2,br2,bi2,cr2,ci2);
    u := Make_Complex(cr,ci); u2 := Make_Complex(cr2,ci2);
    v := Make_Complex(fr,fi); v2 := Make_Complex(fr2,fi2);
    w := Make_Complex(br,bi); w2 := Make_Complex(br2,bi2);
    put_line("the forward products : "); put_line(f);
    put_line("the forward products with loop fusion : "); put_line(f2);
    err := Evaluation_Differentiation_Errors.Difference(f,f2);
    put("The error :"); put(err,3); new_line;
    sumerr := err;
    put_line("forward products recomputed : "); put_line(v);
    err := Evaluation_Differentiation_Errors.Difference(f,v);
    put("The error :"); put(err,3); new_line;
    sumerr := sumerr + err;
    put_line("forward products recomputed with loop fusion : "); put_line(v2);
    err := Evaluation_Differentiation_Errors.Difference(f,v2);
    put("The error :"); put(err,3); new_line;
    sumerr := sumerr + err;
    put_line("the backward products : "); put_line(b);
    put_line("the backward products with loop fusion : "); put_line(b2);
    err := Evaluation_Differentiation_Errors.Difference(b,b2);
    put("The error :"); put(err,3); new_line;
    sumerr := sumerr + err;
    put_line("backward products recomputed : "); put_line(w);
    err := Evaluation_Differentiation_Errors.Difference(b,w);
    put("The error :"); put(err,3); new_line;
    sumerr := sumerr + err;
    put_line("backward products recomputed with loop fusion : "); put_line(w2);
    err := Evaluation_Differentiation_Errors.Difference(b,w2);
    put("The error :"); put(err,3); new_line;
    sumerr := sumerr + err;
    put_line("the cross products : "); put_line(c);
    put_line("the cross products with loop fusion : "); put_line(c2);
    err := Evaluation_Differentiation_Errors.Difference(c,c2);
    put("The error :"); put(err,3); new_line;
    sumerr := sumerr + err;
    put_line("cross products recomputed : "); put_line(u);
    err := Evaluation_Differentiation_Errors.Difference(c,u);
    put("The error :"); put(err,3); new_line;
    sumerr := sumerr + err;
    put_line("cross products recomputed wth loop fusion: "); put_line(u2);
    err := Evaluation_Differentiation_Errors.Difference(c,u2);
    put("The error :"); put(err,3); new_line;
    sumerr := sumerr + err;
    put("The sum of all errors :"); put(sumerr,3); new_line;
  end Test_Forward_Backward_Cross;

  function Random_Indices
             ( dim : integer32 )
             return Standard_Integer_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of indices between 1 and dim
  --   of variables participating in a product.

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

  function Random_Complex_Circuit
             ( nbr,dim : integer32 )
             return Standard_Complex_Circuits.Circuit is

  -- DESCRIPTION :
  --   Returns a random complex circuit with products of dimension dim
  --   and as many nonconstant coefficients as the number nbr.

    res : Standard_Complex_Circuits.Circuit(nbr)
        := Standard_Complex_Circuits.Allocate(nbr,dim);

  begin
    for k in 1..nbr loop
      res.xps(k) := new Standard_Integer_Vectors.Vector'(Random_Indices(dim));
    end loop;
    res.cff := Standard_Random_Vectors.Random_Vector(1,nbr);
    res.cst := Standard_Random_Numbers.Random1;
    return res;
  end Random_Complex_Circuit;

  function Random_Complex_Circuit
             ( nbr,dim,pwr : integer32 )
             return Standard_Complex_Circuits.Circuit is

  -- DESCRIPTION :
  --   Returns a random complex circuit with products of dimension dim,
  --   as many nonconstant coefficients as the number nbr,
  --   and pwr as the value for the higest power.

    res : Standard_Complex_Circuits.Circuit(nbr)
        := Standard_Complex_Circuits.Allocate(nbr,dim);
    xpk : Standard_Integer_Vectors.Vector(1..dim);

  begin
    for k in 1..nbr loop
      xpk := Standard_Random_Vectors.Random_Vector(1,dim,0,pwr);
      res.xps(k) := new Standard_Integer_Vectors.Vector'(xpk);
      res.idx(k) := Exponent_Indices.Exponent_Index(res.xps(k));
      res.fac(k) := Exponent_Indices.Factor_Index(res.xps(k));
    end loop;
    res.cff := Standard_Random_Vectors.Random_Vector(1,nbr);
    res.cst := Standard_Random_Numbers.Random1;
    return res;
  end Random_Complex_Circuit;

  function Split ( c : Standard_Complex_Circuits.Circuit )
                 return Standard_Coefficient_Circuits.Circuit is

  -- DESCRIPTION :
  --   Returns the circuit c with complex coefficients split into
  --   real and imaginary parts.

    res : Standard_Coefficient_Circuits.Circuit(c.nbr)
        := Standard_Coefficient_Circuits.Allocate(c.nbr,c.dim);

  begin
    for k in 1..c.nbr loop
      res.xps(k) := new Standard_Integer_Vectors.Vector'(c.xps(k).all);
      res.idx(k) := Exponent_Indices.Exponent_Index(res.xps(k));
      res.fac(k) := Exponent_Indices.Factor_Index(res.xps(k));
    end loop;
    Split_Complex(c.cff,res.rcf,res.icf);
    res.rcst := Standard_Complex_Numbers.REAL_PART(c.cst);
    res.icst := Standard_Complex_Numbers.IMAG_PART(c.cst);
    return res;
  end Split;

  function Split ( c : Standard_Complex_Circuits.Circuits )
                 return Standard_Coefficient_Circuits.Circuits is

    res : Standard_Coefficient_Circuits.Circuits(c'range);

  begin
    for k in c'range loop
      declare
        ck : constant Standard_Complex_Circuits.Circuit := c(k).all;
      begin
        res(k) := new Standard_Coefficient_Circuits.Circuit'(Split(ck));
      end;
    end loop;
    return res;
  end Split;

  function Split ( s : Standard_Complex_Circuits.System )
                 return Standard_Coefficient_Circuits.System is

    crc : constant Standard_Coefficient_Circuits.Circuits(1..s.neq)
        := Split(s.crc);
    res : constant Standard_Coefficient_Circuits.System(s.neq,s.dim)
        := Standard_Coefficient_Circuits.Create(crc,s.dim);

  begin
    return res;
  end Split;

  function Split ( s : Standard_Complex_Circuits.Link_to_System )
                 return Standard_Coefficient_Circuits.Link_to_System is

    res : constant Standard_Coefficient_Circuits.Link_to_System
        := new Standard_Coefficient_Circuits.System'(Split(s.all));

  begin
    return res;
  end Split;

  procedure Test_Indexed_Forward_Backward_Cross ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random vector of dimension dim
  --   and tests the computation of the forward/backward/cross products.

    idx : constant Standard_Integer_Vectors.Vector := Random_Indices(dim);
    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    size : constant integer32 := idx'last;
    cf : constant Standard_Complex_Vectors.Vector(1..size-1)
       := Standard_Complex_Vectors.Vector'(1..size-1 => zero);
    cf2 : constant Standard_Complex_Vectors.Vector(1..size-1)
        := Standard_Complex_Vectors.Vector'(1..size-1 => zero);
    cb : constant Standard_Complex_Vectors.Vector(1..size-2)
        := Standard_Complex_Vectors.Vector'(1..size-2 => zero);
    cb2 : constant Standard_Complex_Vectors.Vector(1..size-2)
        := Standard_Complex_Vectors.Vector'(1..size-2 => zero);
    cc : constant Standard_Complex_Vectors.Vector(1..size-2)
        := Standard_Complex_Vectors.Vector'(1..size-2 => zero);
    cc2 : constant Standard_Complex_Vectors.Vector(1..size-2)
        := Standard_Complex_Vectors.Vector'(1..size-2 => zero);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    f : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cf);
    f2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cf2);
    b : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cb);
    b2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cb2);
    c : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cc);
    c2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cc2);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    fr2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(f2);
    fi2 : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(f2);
    br2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(b2);
    bi2 : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(b2);
    cr2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(c2);
    ci2 : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(c2);
    u2,v2,w2 : Standard_Complex_Vectors.Link_to_Vector;
    err,sumerr : double_float;

  begin
    put("The indices in the product : "); put(idx); new_line;
    put("idx'last : "); put(idx'last,1); new_line;
    Standard_Complex_Circuits.Forward_Backward_Cross(idx,x,f,b,c);
    Standard_Coefficient_Circuits.Fused_Forward_Backward_Cross
      (idx,xr,xi,fr2,fi2,br2,bi2,cr2,ci2);
    u2 := Make_Complex(cr2,ci2);
    v2 := Make_Complex(fr2,fi2);
    w2 := Make_Complex(br2,bi2);
    put_line("the forward products : "); put_line(f);
    put_line("forward products recomputed with loop fusion : "); put_line(v2);
    err := Evaluation_Differentiation_Errors.Difference(f,v2);
    put("The error :"); put(err,3); new_line;
    sumerr := err;
    put_line("the backward products : "); put_line(b);
    put_line("backward products recomputed with loop fusion : "); put_line(w2);
    err := Evaluation_Differentiation_Errors.Difference(b,w2);
    put("The error :"); put(err,3); new_line;
    sumerr := sumerr + err;
    put_line("the cross products : "); put_line(c);
    put_line("cross products recomputed wth loop fusion: "); put_line(u2);
    err := Evaluation_Differentiation_Errors.Difference(c,u2);
    put("The error :"); put(err,3); new_line;
    sumerr := sumerr + err;
    put("The sum of all errors :"); put(sumerr,3); new_line;
  end Test_Indexed_Forward_Backward_Cross;

  procedure Test_Power_Table ( dim,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Given the dimension dim and the highest power pwr,
  --   tests the computation of the power table for random values.

    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim) := (1..dim => pwr);
    pwt : constant Standard_Complex_VecVecs.VecVec(x'range)
        := Standard_Complex_Circuits.Allocate(mxe);
    rpwt : constant Standard_Floating_VecVecs.VecVec(x'range)
         := Standard_Coefficient_Circuits.Allocate(mxe);
    ipwt : constant Standard_Floating_VecVecs.VecVec(x'range)
         := Standard_Coefficient_Circuits.Allocate(mxe);
    v : Standard_Complex_VecVecs.VecVec(x'range);
    err : double_float;

  begin
    Standard_Complex_Circuits.Power_Table(mxe,x,pwt);
    put_line("The power table : "); put_line(pwt);
    Standard_Coefficient_Circuits.Power_Table(mxe,xr,xi,rpwt,ipwt);
    v := Standard_Vector_Splitters.Make_Complex(rpwt,ipwt);
    put_line("The recomputed power table : "); put_line(v);
    err := Evaluation_Differentiation_Errors.Difference(pwt,v);
    put("The error :"); put(err,3); new_line;
  end Test_Power_Table;

  function Make_Polynomial
             ( c : Standard_Complex_Circuits.Circuit;
               index : boolean := false )
             return Standard_Complex_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns the polynomial equivalent to the circuit c.
  --   If index, then c.xps should be considered as an index.

    use Standard_Complex_Polynomials;

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

  function Gradient ( p : Standard_Complex_Polynomials.Poly;
                      x : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Straighforward computation of the gradient of p at x,
  --   for testing purposes.

    res : Standard_Complex_Vectors.Vector(x'range);
    dp : Standard_Complex_Polynomials.Poly;

  begin
    for k in x'range loop
      dp := Standard_Complex_Polynomials.Diff(p,k);
      res(k) := Standard_Complex_Poly_Functions.Eval(dp,x);
      Standard_Complex_Polynomials.Clear(dp);
    end loop;
    return res;
  end Gradient;

  procedure Test_Circuit ( nbr,dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random circuit with nbr terms and dimension dim
  --   and tests the differentiation and evaluation.

    c1 : constant Standard_Complex_Circuits.Circuit
       := Random_Complex_Circuit(nbr,dim);
    p : constant Standard_Complex_Polynomials.Poly := Make_Polynomial(c1,true);
    c2 : constant Standard_Coefficient_Circuits.Circuit := Split(c1);
    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    y : constant Standard_Complex_Vectors.Vector(0..dim)
      := Standard_Complex_Vectors.Vector'(0..dim => zero);
    yd : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(y);
    ryd : constant Standard_Floating_Vectors.Link_to_Vector
        := Real_Part(yd);
    iyd : constant Standard_Floating_Vectors.Link_to_Vector
        := Imag_Part(yd);
    yd2 : Standard_Complex_Vectors.Link_to_Vector;
    z : Standard_Complex_Numbers.Complex_Number;
    zd : Standard_Complex_Vectors.Vector(1..dim);
    err : double_float;

  begin
    Standard_Complex_Circuits.Speel(c1,x,yd);
    put_line("The value at a random point :"); put(yd(0)); new_line;
    z := Standard_Complex_Poly_Functions.Eval(p,cx);
    put_line("The value recomputed for testing :"); put(z); new_line;
    put_line("The gradient :"); put_line(yd(1..yd'last));
    zd := Gradient(p,cx);
    put_line("The gradient recomputed for testing :");
    put_line(zd);
   -- err := Evaluation_Differentiation_Errors.Difference(yd(1..yd'last),zd);
   -- put("The error :"); put(err,3); new_line;
    Standard_Coefficient_Circuits.Speel(c2,xr,xi,ryd,iyd);
    yd2 := Make_Complex(ryd,iyd);
    put_line("The recomputed value :"); put(yd2(0)); new_line;
    put_line("The recomputed gradient :"); put_line(yd2(1..yd2'last));
    err := Evaluation_Differentiation_Errors.Difference(yd,yd2);
    put("The error :"); put(err,3); new_line;
  end Test_Circuit;

  procedure Test_Multiply_Factor ( dim,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Given the dimension dim and the highest power pwr,
  --   tests the multiplication of the common factor,
  --   using the power table for random values.

    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim) := (1..dim => pwr);
    pwt : constant Standard_Complex_VecVecs.VecVec(x'range)
        := Standard_Complex_Circuits.Allocate(mxe);
    rpwt : constant Standard_Floating_VecVecs.VecVec(x'range)
         := Standard_Coefficient_Circuits.Allocate(mxe);
    ipwt : constant Standard_Floating_VecVecs.VecVec(x'range)
         := Standard_Coefficient_Circuits.Allocate(mxe);
    xps : constant Standard_Integer_Vectors.Vector(1..dim)
        := Standard_Random_Vectors.Random_Vector(1,dim,0,pwr);
    xp : constant Standard_Integer_Vectors.Link_to_Vector
       := new Standard_Integer_Vectors.Vector'(xps);
    fac : constant Standard_Integer_Vectors.Link_to_Vector
        := Exponent_Indices.Factor_Index(xp);
    cst : constant Standard_Complex_Numbers.Complex_Number
        := Standard_Random_Numbers.Random1;
    rcf : constant double_float := Standard_Complex_Numbers.REAL_PART(cst);
    icf : constant double_float := Standard_Complex_Numbers.IMAG_PART(cst);
    res,res2 : Standard_Complex_Numbers.Complex_Number;
    rpf,ipf : double_float;

  begin
    put("The random exponents : "); put(xps); new_line;
    put("The factor indices : "); put(fac); new_line;
    Standard_Complex_Circuits.Power_Table(mxe,x,pwt);
    Standard_Coefficient_Circuits.Power_Table(mxe,xr,xi,rpwt,ipwt);
    Standard_Complex_Circuits.Multiply_Factor(xp,fac,x,cst,pwt,res);
    put_line("The multiplied factor :"); put(res); new_line;
    Standard_Coefficient_Circuits.Multiply_Factor
      (xp,fac,xr,xi,rcf,icf,rpwt,ipwt,rpf,ipf);
    res2 := Standard_Complex_Numbers.Create(rpf,ipf);
    put_line("The recomputed multiplied factor :"); put(res2); new_line;
  end Test_Multiply_Factor;

  procedure Test_Power_Circuit ( nbr,dim,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random circuit with nbr terms, dimension dim,
  --   highest power pwr, and tests the differentiation and evaluation.

    c1 : constant Standard_Complex_Circuits.Circuit
       := Random_Complex_Circuit(nbr,dim,pwr);
    p : constant Standard_Complex_Polynomials.Poly := Make_Polynomial(c1);
    c2 : constant Standard_Coefficient_Circuits.Circuit := Split(c1);
    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    y : constant Standard_Complex_Vectors.Vector(0..dim)
      := Standard_Complex_Vectors.Vector'(0..dim => zero);
    yd : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(y);
    ryd : constant Standard_Floating_Vectors.Link_to_Vector
        := Real_Part(yd);
    iyd : constant Standard_Floating_Vectors.Link_to_Vector
        := Imag_Part(yd);
    yd2 : Standard_Complex_Vectors.Link_to_Vector;
    mxe : constant Standard_Integer_Vectors.Vector(1..dim)
        := Exponent_Indices.Maxima(c1.xps);
    pwt : constant Standard_Complex_VecVecs.VecVec(x'range)
        := Standard_Complex_Circuits.Allocate(mxe);
    rpwt : constant Standard_Floating_VecVecs.VecVec(x'range)
         := Standard_Coefficient_Circuits.Allocate(mxe);
    ipwt : constant Standard_Floating_VecVecs.VecVec(x'range)
         := Standard_Coefficient_Circuits.Allocate(mxe);
    err : double_float;
    z : Standard_Complex_Numbers.Complex_Number;
    zd : Standard_Complex_Vectors.Vector(cx'range);

  begin
    Standard_Complex_Circuits.Power_Table(mxe,x,pwt);
    Standard_Coefficient_Circuits.Power_Table(mxe,xr,xi,rpwt,ipwt);
    Standard_Complex_Circuits.Speel(c1,x,yd,pwt);
    put_line("The value at a random point :"); put(yd(0)); new_line;
    z := Standard_Complex_Poly_Functions.Eval(p,cx);
    put_line("The value recomputed for testing :"); put(z); new_line;
    put_line("The gradient :"); put_line(yd(1..yd'last));
    zd := Gradient(p,cx);
    put_line("The gradient recomputed for testing :"); put_line(zd);
    Standard_Coefficient_Circuits.Speel(c2,xr,xi,ryd,iyd,rpwt,ipwt);
    yd2 := Make_Complex(ryd,iyd);
    put_line("The recomputed value :"); put(yd2(0)); new_line;
    put_line("The recomputed gradient :"); put_line(yd2(1..yd2'last));
    err := Evaluation_Differentiation_Errors.Difference(yd,yd2);
    put("The error :"); put(err,3); new_line;
  end Test_Power_Circuit;

  procedure Timing_Forward ( dim,frq : in integer32 ) is

  -- DESCRIPTION :
  --   Does as many forward product computations as freq
  --   on random vectors of dimension dim.

    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    cf : constant Standard_Complex_Vectors.Vector(1..dim-1)
       := Standard_Complex_Vectors.Vector'(1..dim-1 => zero);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    f : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cf);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    fr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(f);
    fi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(f);
    timer : Timing_Widget;

  begin
    tstart(timer);
    for k in 1..frq loop
      Standard_Complex_Circuits.Forward(x,f);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex forward products");
    tstart(timer);
    for k in 1..frq loop
      Standard_Coefficient_Circuits.Forward(xr,xi,fr,fi);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"real forward products");
  end Timing_Forward;

  procedure Timing_Forward_Backward ( dim,frq : in integer32 ) is

  -- DESCRIPTION :
  --   Does as many forward/backward product computations as freq
  --   on random vectors of dimension dim.

    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    cf : constant Standard_Complex_Vectors.Vector(1..dim-1)
       := Standard_Complex_Vectors.Vector'(1..dim-1 => zero);
    cf2 : constant Standard_Complex_Vectors.Vector(1..dim-1)
        := Standard_Complex_Vectors.Vector'(1..dim-1 => zero);
    cb : constant Standard_Complex_Vectors.Vector(1..dim-2)
       := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    cb2 : constant Standard_Complex_Vectors.Vector(1..dim-2)
        := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    f : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cf);
    f2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cf2);
    b : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cb);
    b2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cb2);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    fr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(f);
    fi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(f);
    br : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(b);
    bi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(b);
    fr2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(f2);
    fi2 : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(f2);
    br2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(b2);
    bi2 : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(b2);
    timer : Timing_Widget;

  begin
    tstart(timer);
    for k in 1..frq loop
      Standard_Complex_Circuits.Forward_Backward(x,f,b);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex forward & backward products");
    tstart(timer);
    for k in 1..frq loop
      Standard_Complex_Circuits.Fused_Forward_Backward(x,f2,b2);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex with loop fusion");
    tstart(timer);
    for k in 1..frq loop
      Standard_Coefficient_Circuits.Forward_Backward(xr,xi,fr,fi,br,bi);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"real forward & backward products");
    tstart(timer);
    for k in 1..frq loop
      Standard_Coefficient_Circuits.Fused_Forward_Backward
        (xr,xi,fr2,fi2,br2,bi2);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"real with loop fusion");
  end Timing_Forward_Backward;

  procedure Timing_Forward_Backward_Cross ( dim,frq : in integer32 ) is

  -- DESCRIPTION :
  --   Does as many forward/backward/cross product computations as freq
  --   on random vectors of dimension dim.

    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    cf : constant Standard_Complex_Vectors.Vector(1..dim-1)
       := Standard_Complex_Vectors.Vector'(1..dim-1 => zero);
    cf2 : constant Standard_Complex_Vectors.Vector(1..dim-1)
        := Standard_Complex_Vectors.Vector'(1..dim-1 => zero);
    cb : constant Standard_Complex_Vectors.Vector(1..dim-2)
       := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    cb2 : constant Standard_Complex_Vectors.Vector(1..dim-2)
        := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    cc : constant Standard_Complex_Vectors.Vector(1..dim-2)
       := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    cc2 : constant Standard_Complex_Vectors.Vector(1..dim-2)
        := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    f : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cf);
    f2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cf2);
    b : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cb);
    b2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cb2);
    c : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cc);
    c2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cc2);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    fr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(f);
    fi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(f);
    br : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(b);
    bi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(b);
    cr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(c);
    ci : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(c);
    fr2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(f2);
    fi2 : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(f2);
    br2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(b2);
    bi2 : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(b2);
    cr2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(c2);
    ci2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(c2);
    timer : Timing_Widget;

  begin
    tstart(timer);
    for k in 1..frq loop
      Standard_Complex_Circuits.Forward_Backward_Cross(x,f,b,c);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex forward, backward, cross ");
    tstart(timer);
    for k in 1..frq loop
      Standard_Complex_Circuits.Fused_Forward_Backward_Cross(x,f2,b2,c2);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"fused complex forward, backward, cross");
    tstart(timer);
    for k in 1..frq loop
      Standard_Coefficient_Circuits.Forward_Backward_Cross
        (xr,xi,fr,fi,br,bi,cr,ci);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"real forward, backward, cross");
    tstart(timer);
    for k in 1..frq loop
      Standard_Coefficient_Circuits.Fused_Forward_Backward_Cross
        (xr,xi,fr2,fi2,br2,bi2,cr2,ci2);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"fused real forward, backward, cross");
  end Timing_Forward_Backward_Cross;

  procedure Timing_Indexed_Forward_Backward_Cross
              ( dim,frq : in integer32 ) is

  -- DESCRIPTION :
  --   Does as many forward/backward/cross product computations as freq
  --   on random vectors of dimension dim.

    idx : constant Standard_Integer_Vectors.Vector := Random_Indices(dim);
    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    cf : constant Standard_Complex_Vectors.Vector(1..dim-1)
       := Standard_Complex_Vectors.Vector'(1..dim-1 => zero);
    cf2 : constant Standard_Complex_Vectors.Vector(1..dim-1)
        := Standard_Complex_Vectors.Vector'(1..dim-1 => zero);
    cb : constant Standard_Complex_Vectors.Vector(1..dim-2)
       := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    cb2 : constant Standard_Complex_Vectors.Vector(1..dim-2)
        := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    cc : constant Standard_Complex_Vectors.Vector(1..dim-2)
       := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    cc2 : constant Standard_Complex_Vectors.Vector(1..dim-2)
        := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    f : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cf);
    f2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cf2);
    b : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cb);
    b2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cb2);
    c : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cc);
    c2 : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(cc2);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    fr2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(f2);
    fi2 : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(f2);
    br2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(b2);
    bi2 : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(b2);
    cr2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(c2);
    ci2 : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(c2);
    timer : Timing_Widget;

  begin
    tstart(timer);
    for k in 1..frq loop
      Standard_Complex_Circuits.Forward_Backward_Cross(idx,x,f,b,c);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex forward, backward, cross ");
    tstart(timer);
    for k in 1..frq loop
      Standard_Coefficient_Circuits.Fused_Forward_Backward_Cross
        (idx,xr,xi,fr2,fi2,br2,bi2,cr2,ci2);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"fused real forward, backward, cross");
  end Timing_Indexed_Forward_Backward_Cross;

  procedure Timing_Power_Table ( dim,pwr,frq : in integer32 ) is

  -- DESCRIPTION :
  --   Given the dimension dim, the highest power pwr,
  --   and the frequency frq, times the computation of the power table.

    timer : Timing_Widget;
    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim) := (1..dim => pwr);
    pwt : constant Standard_Complex_VecVecs.VecVec(x'range)
        := Standard_Complex_Circuits.Allocate(mxe);
    rpwt : constant Standard_Floating_VecVecs.VecVec(x'range)
         := Standard_Coefficient_Circuits.Allocate(mxe);
    ipwt : constant Standard_Floating_VecVecs.VecVec(x'range)
         := Standard_Coefficient_Circuits.Allocate(mxe);

  begin
    tstart(timer);
    for k in 1..frq loop
      Standard_Complex_Circuits.Power_Table(mxe,x,pwt);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex power table");
    tstart(timer);
    for k in 1..frq loop
      Standard_Coefficient_Circuits.Power_Table(mxe,xr,xi,rpwt,ipwt);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"real power table");
  end Timing_Power_Table;

  procedure Timing_Circuit ( nbr,dim,frq: in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random circuit with nbr terms and dimension dim
  --   and runs the differentiation and evaluation procedures
  --   as many times as the frequency frq.

    timer : Timing_Widget;
    c1 : constant Standard_Complex_Circuits.Circuit
       := Random_Complex_Circuit(nbr,dim);
    c2 : constant Standard_Coefficient_Circuits.Circuit := Split(c1);
    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    y : constant Standard_Complex_Vectors.Vector(0..dim)
      := Standard_Complex_Vectors.Vector'(0..dim => zero);
    yd : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(y);
    ryd : constant Standard_Floating_Vectors.Link_to_Vector
        := Real_Part(yd);
    iyd : constant Standard_Floating_Vectors.Link_to_Vector
        := Imag_Part(yd);

  begin
    tstart(timer);
    for k in 1..frq loop
      Standard_Complex_Circuits.Speel(c1,x,yd);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex Speel");
    tstart(timer);
    for k in 1..frq loop
      Standard_Coefficient_Circuits.Speel(c2,xr,xi,ryd,iyd);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"real Speel");
  end Timing_Circuit;

  procedure Timing_Multiply_Factor ( dim,pwr,frq : in integer32 ) is

  -- DESCRIPTION :
  --   Given the dimension dim and the highest power pwr,
  --   times the multiplication of the common factor,
  --   using the power table for random values.

    timer : Timing_Widget;
    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim) := (1..dim => pwr);
    pwt : constant Standard_Complex_VecVecs.VecVec(x'range)
        := Standard_Complex_Circuits.Allocate(mxe);
    rpwt : constant Standard_Floating_VecVecs.VecVec(x'range)
         := Standard_Coefficient_Circuits.Allocate(mxe);
    ipwt : constant Standard_Floating_VecVecs.VecVec(x'range)
         := Standard_Coefficient_Circuits.Allocate(mxe);
    xps : constant Standard_Integer_Vectors.Vector(1..dim)
        := Standard_Random_Vectors.Random_Vector(1,dim,0,pwr);
    xp : constant Standard_Integer_Vectors.Link_to_Vector
       := new Standard_Integer_Vectors.Vector'(xps);
    fac : constant Standard_Integer_Vectors.Link_to_Vector
        := Exponent_Indices.Factor_Index(xp);
    cst : Standard_Complex_Numbers.Complex_Number;
    res,res2 : Standard_Complex_Numbers.Complex_Number;
    rcf,icf,rpf,ipf : double_float;

  begin
    put("The random exponents : "); put(xps); new_line;
    put("The factor indices : "); put(fac); new_line;
    Standard_Complex_Circuits.Power_Table(mxe,x,pwt);
    Standard_Coefficient_Circuits.Power_Table(mxe,xr,xi,rpwt,ipwt);
    tstart(timer);
    for k in 1..frq loop
      cst := Standard_Random_Numbers.Random1;
      Standard_Complex_Circuits.Multiply_Factor(xp,fac,x,cst,pwt,res);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex multiply factor");
    tstart(timer);
    for k in 1..frq loop
      cst := Standard_Random_Numbers.Random1;
      rpf := Standard_Complex_Numbers.REAL_PART(cst);
      ipf := Standard_Complex_Numbers.IMAG_PART(cst);
      Standard_Coefficient_Circuits.Multiply_Factor
        (xp,fac,xr,xi,rcf,icf,rpwt,ipwt,rpf,ipf);
      res2 := Standard_Complex_Numbers.Create(rpf,ipf);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"real multiply factor");
  end Timing_Multiply_Factor;

  procedure Timing_Power_Circuit ( nbr,dim,pwr,frq: in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random circuit with nbr terms, dimension dim, highest
  --   power pwr, and runs the differentiation and evaluation procedures
  --   as many times as the frequency frq.

    timer : Timing_Widget;
    c1 : constant Standard_Complex_Circuits.Circuit
       := Random_Complex_Circuit(nbr,dim,pwr);
    c2 : constant Standard_Coefficient_Circuits.Circuit := Split(c1);
    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    y : constant Standard_Complex_Vectors.Vector(0..dim)
      := Standard_Complex_Vectors.Vector'(0..dim => zero);
    yd : constant Standard_Complex_Vectors.Link_to_Vector
       := new Standard_Complex_Vectors.Vector'(y);
    ryd : constant Standard_Floating_Vectors.Link_to_Vector
        := Real_Part(yd);
    iyd : constant Standard_Floating_Vectors.Link_to_Vector
        := Imag_Part(yd);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim)
        := Exponent_Indices.Maxima(c1.xps);
    pwt : constant Standard_Complex_VecVecs.VecVec(x'range)
        := Standard_Complex_Circuits.Allocate(mxe);
    rpwt : constant Standard_Floating_VecVecs.VecVec(x'range)
         := Standard_Coefficient_Circuits.Allocate(mxe);
    ipwt : constant Standard_Floating_VecVecs.VecVec(x'range)
         := Standard_Coefficient_Circuits.Allocate(mxe);

  begin
    Standard_Complex_Circuits.Power_Table(mxe,x,pwt);
    Standard_Coefficient_Circuits.Power_Table(mxe,xr,xi,rpwt,ipwt);
    tstart(timer);
    for k in 1..frq loop
      Standard_Complex_Circuits.Speel(c1,x,yd,pwt);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex Speel");
    tstart(timer);
    for k in 1..frq loop
      Standard_Coefficient_Circuits.Speel(c2,xr,xi,ryd,iyd,rpwt,ipwt);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"real Speel");
  end Timing_Power_Circuit;

  function Constant_Coefficient
             ( p : Standard_Complex_Polynomials.Poly )
             return Standard_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the constant coefficient of the polynomial p.

    use Standard_Complex_Polynomials;

    dim : constant integer32 := integer32(Number_of_Unknowns(p));
    deg : Degrees := new Standard_Natural_Vectors.Vector'(1..dim => 0);
    res : constant Standard_Complex_Numbers.Complex_Number := Coeff(p,deg);

  begin
    Clear(deg);
    return res;
  end Constant_Coefficient;

  function Is_NonZero ( c : Standard_Complex_Numbers.Complex_Number )
                   return integer32 is

  -- DESCRIPTION :
  --   Returns 1 if the coefficient is nonzero,
  --   returns 0 otherwise.

  begin
    if Standard_Complex_Numbers.REAL_PART(c) /= 0.0 then
      return 1;
    elsif Standard_Complex_Numbers.IMAG_PART(c) /= 0.0 then
      return 1;
    else
      return 0;
    end if;
  end Is_NonZero;

  function Make_Complex_Circuit
             ( p : Standard_Complex_Polynomials.Poly )
             return Standard_Complex_Circuits.Circuit is

  -- DESCRIPTION :
  --   Returns the circuit representation of the polynomial p.
    
    use Standard_Complex_Polynomials;

    nbr : constant integer32 := integer32(Number_of_Terms(p));
    dim : constant integer32 := integer32(Number_of_Unknowns(p));
    cst : constant Standard_Complex_Numbers.Complex_Number
        := Constant_Coefficient(p);
    isz : constant integer32 := Is_NonZero(cst);
    res : Standard_Complex_Circuits.Circuit(nbr-isz)
        := Standard_Complex_Circuits.Allocate(nbr-isz,dim);
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
    res.cst := cst;
    Visit_Terms(p);
    return res;
  end Make_Complex_Circuit;

  function Make_Complex_System
             ( p : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
               verbose : in boolean := true )
             return Standard_Complex_Circuits.Link_to_System is

  -- DESCRIPTION :
  --   Returns the system of circuits defined by p.
  --   If verbose, then the tableau format of the system is written.

    use Standard_Complex_Circuits;

    res : Link_to_System;
    c : Circuits(p'range);
    d : integer32;

  begin
    for k in c'range loop
      c(k) := new Circuit'(Make_Complex_Circuit(p(k)));
      if verbose then
        for i in 1..c(k).nbr loop
          put(c(k).cff(i)); put(c(k).xps(i)); new_line;
        end loop;
        put(c(k).cst); new_line;
      end if;
    end loop;
    d := c(c'first).dim;
    res := new System'(Create(c,d));
    return res;
  end Make_Complex_System;

  procedure Write_Matrix ( A : in Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Writes the matrix A component-wise, with explicit indexing.

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("A["); put(i,1); put(","); put(j,1); put("] : ");
        put(A(i,j)); new_line;
      end loop;
    end loop;
  end Write_Matrix;

  function Sum_of_Errors
             ( x,y : in Standard_Complex_Vectors.Vector )
             return double_float is

  -- DESCRIPTION :
  --   Returns the sum of the component-wise differences between
  --   the vectors x and y.

  -- REQUIRED : x'range = y'range.

    use Standard_Complex_numbers;

    res : double_float := 0.0;
    val : Complex_Number;

  begin
    for i in x'range loop
      val := x(i) - y(i);
      res := res + AbsVal(val);
    end loop;
    return res;
  end Sum_of_Errors;

  function Sum_of_Errors
             ( A,B : in Standard_Complex_Matrices.Matrix )
             return double_float is

  -- DESCRIPTION :
  --   Returns the sum of the component-wise differences between
  --   the matrices A and B.

  -- REQUIRED : A'range(1) = B'range(1) and A'range(2) = B'range(2).

    use Standard_Complex_numbers;

    res : double_float := 0.0;
    val : Complex_Number;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        val := A(i,j) - B(i,j);
        res := res + AbsVal(val);
      end loop;
    end loop;
    return res;
  end Sum_of_Errors;

  procedure Test_Evaluation_and_Differentiation
              ( p : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                cs : in Standard_Complex_Circuits.Link_to_System;
                cffsys : in Standard_Coefficient_Circuits.Link_to_System ) is

  -- DESCRIPTION :
  --   Tests the evaluation and differentiation at a random point,
  --   for the polynomial system p and systems of circuits cs and cffsys.

    cx : constant Standard_Complex_Vectors.Vector(1..cs.dim)
       := Standard_Random_Vectors.Random_Vector(1,cs.dim);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    y : constant Standard_Complex_Vectors.Vector
      := Standard_Complex_Poly_SysFun.Eval(p.all,cx);
    jm : Standard_Complex_Jaco_Matrices.Jaco_Mat(p'range,cx'range)
       := Standard_Complex_Jaco_Matrices.Create(p.all);
    jmx : constant Standard_Complex_Matrices.Matrix(p'range,cx'range)
        := Standard_Complex_Jaco_Matrices.Eval(jm,cx);
    err : double_float;

  begin
    Standard_Complex_Circuits.EvalDiff(cs,x);
    Standard_Coefficient_Circuits.EvalDiff(cffsys,xr,xi);
    put_line("The value at a random point :"); put_line(cs.fx);
    put_line("The value recomputed for testing :"); put_line(y);
    err := Sum_of_Errors(cs.fx,y);
    put("Sum of errors :"); put(err,3); new_line;
    put_line("The recomputed value :"); put_line(cffsys.fx);
    err := Sum_of_Errors(cffsys.fx,y);
    put("Sum of errors :"); put(err,3); new_line;
    put_line("The evaluated Jacobian matrix :"); Write_Matrix(cs.jm);
    put_line("The matrix recomputed for testing :"); Write_Matrix(jmx);
    err := Sum_of_Errors(cs.jm,jmx);
    put("Sum of errors :"); put(err,3); new_line;
    put_line("The recomputed matrix :"); Write_Matrix(cffsys.jm);
    err := Sum_of_Errors(cffsys.jm,jmx);
    put("Sum of errors :"); put(err,3); new_line;
    Standard_Complex_Jaco_Matrices.Clear(jm);
  end Test_Evaluation_and_Differentiation;

  procedure Timing_Evaluation_and_Differentiation
              ( p : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                cs : in Standard_Complex_Circuits.Link_to_System;
                cffsys : in Standard_Coefficient_Circuits.Link_to_System;
                frq : in integer32 ) is

  -- DESCRIPTION :
  --   Times the evaluation and differentiation at a random point,
  --   for the polynomial system p and systems of circuits cs and cffsys,
  --   with frequency frq.

    timer : Timing_Widget;
    cx : constant Standard_Complex_Vectors.Vector(1..cs.dim)
       := Standard_Random_Vectors.Random_Vector(1,cs.dim);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    y : Standard_Complex_Vectors.Vector(p'range);
    ep : Standard_Complex_Poly_SysFun.Eval_Poly_Sys(p'range)
       := Standard_Complex_Poly_SysFun.Create(p.all);
    jm : Standard_Complex_Jaco_Matrices.Jaco_Mat(p'range,cx'range)
       := Standard_Complex_Jaco_Matrices.Create(p.all);
    ejm : Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat(p'range,cx'range)
        := Standard_Complex_Jaco_Matrices.Create(jm);
    jmx : Standard_Complex_Matrices.Matrix(p'range,cx'range);

  begin
    tstart(timer);
    for k in 1..frq loop
      y := Standard_Complex_Poly_SysFun.Eval(ep,cx);
      jmx := Standard_Complex_Jaco_Matrices.Eval(ejm,cx);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"nested Horner");
    tstart(timer);
    for k in 1..frq loop
      Standard_Complex_Circuits.EvalDiff(cs,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex circuit ade");
    tstart(timer);
    for k in 1..frq loop
      Standard_Coefficient_Circuits.EvalDiff(cffsys,xr,xi);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"real circuit ade");
    Standard_Complex_Poly_SysFun.Clear(ep);
    Standard_Complex_Jaco_Matrices.Clear(ejm);
    Standard_Complex_Jaco_Matrices.Clear(jm);
  end Timing_Evaluation_and_Differentiation;

  procedure Test_System ( frq : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system
  --   and then evaluates the system at a random point.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    cmpsys : Standard_Complex_Circuits.Link_to_System;
    cffsys : Standard_Coefficient_Circuits.Link_to_System;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(lp);
    cmpsys := Make_Complex_System(lp);
    cffsys := Split(cmpsys);
    if frq = 0
     then Test_Evaluation_and_Differentiation(lp,cmpsys,cffsys);
     else Timing_Evaluation_and_Differentiation(lp,cmpsys,cffsys,frq);
    end if;
  end Test_System;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension,
  --   the type of test, and then launches the test.

    dim,frq,pwr,nbr : integer32 := 0;
    ans,tst : character;

  begin
    new_line;
    put("Test system ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      tst := '9';
    else
      new_line;
      put("Give the dimension of the vectors : "); get(dim);
      new_line;
      put_line("MENU for testing ADE :");
      put_line("  1. forward products");
      put_line("  2. forward and backward products");
      put_line("  3. forward, backward, and cross products");
      put_line("  4. indexed forward, backward, and cross products");
      put_line("  5. power table");
      put_line("  6. circuit differentiation and evaluation");
      put_line("  7. multiplication with common factor");
      put_line("  8. evaluation and differentiaton of random circuit");
      put_line("  9. for given system, test evaluation and differentiation");
      put("Type 1, 2, 3, 4, 5, 6, 7, 8, or 9 to select the test : ");
      Ask_Alternative(tst,"123456789");
    end if;
    if tst = '5' or tst = '7' or tst = '8' then
      new_line;
      put("Give the highest power : "); get(pwr);
    end if;
    if tst = '6' or tst = '8' then
      new_line;
      put("Give the number of terms in the circuit : "); get(nbr);
    end if;
    new_line;
    put("Interactive tests ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      case tst is
        when '1' => Test_Forward(dim);
        when '2' => Test_Forward_Backward(dim);
        when '3' => Test_Forward_Backward_Cross(dim);
        when '4' => Test_Indexed_Forward_Backward_Cross(dim);
        when '5' => Test_Power_Table(dim,pwr);
        when '6' => Test_Circuit(nbr,dim);
        when '7' => Test_Multiply_Factor(dim,pwr);
        when '8' => Test_Power_Circuit(nbr,dim,pwr);
        when '9' => Test_System;
        when others => null;
      end case;
    else
      new_line;
      put("Give the frequency of the tests : "); get(frq);
      case tst is
        when '1' => Timing_Forward(dim,frq);
        when '2' => Timing_Forward_Backward(dim,frq);
        when '3' => Timing_Forward_Backward_Cross(dim,frq);
        when '4' => Timing_Indexed_Forward_Backward_Cross(dim,frq);
        when '5' => Timing_Power_Table(dim,pwr,frq);
        when '6' => Timing_Circuit(nbr,dim,frq);
        when '7' => Timing_Multiply_Factor(dim,pwr,frq);
        when '8' => Timing_Power_Circuit(nbr,dim,pwr,frq);
        when '9' => Test_System(frq);
        when others => null;
      end case;
    end if;
  end Main;

begin
  Main;
end ts_perfade;
