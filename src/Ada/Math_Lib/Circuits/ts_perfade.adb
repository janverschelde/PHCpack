with text_io;                             use text_io;
with Communications_with_User;            use Communications_with_User;
with Timing_Package;                      use Timing_Package;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;        use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;         use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;
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
with Standard_Circuit_Makers;

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

  procedure Test_Indexed_Forward_Backward_Cross ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random vector of dimension dim
  --   and tests the computation of the forward/backward/cross products.

    idx : constant Standard_Integer_Vectors.Vector
        := Standard_Circuit_Makers.Random_Indices(dim);
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

  procedure Test_Circuit ( nbr,dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random circuit with nbr terms and dimension dim
  --   and tests the differentiation and evaluation.

    c1 : constant Standard_Complex_Circuits.Circuit
       := Standard_Circuit_Makers.Random_Complex_Circuit(nbr,dim);
    p : constant Standard_Complex_Polynomials.Poly
      := Standard_Circuit_Makers.Make_Polynomial(c1,true);
    c2 : constant Standard_Coefficient_Circuits.Circuit
       := Standard_Circuit_Makers.Split(c1);
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
    Standard_Complex_Circuits.Indexed_Speel(c1,x,yd);
    put_line("The value at a random point :"); put(yd(0)); new_line;
    z := Standard_Complex_Poly_Functions.Eval(p,cx);
    put_line("The value recomputed for testing :"); put(z); new_line;
    put_line("The gradient :"); put_line(yd(1..yd'last));
    zd := Standard_Circuit_Makers.Gradient(p,cx);
    put_line("The gradient recomputed for testing :");
    put_line(zd);
   -- err := Evaluation_Differentiation_Errors.Difference(yd(1..yd'last),zd);
   -- put("The error :"); put(err,3); new_line;
    Standard_Coefficient_Circuits.Indexed_Speel(c2,xr,xi,ryd,iyd);
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
       := Standard_Circuit_Makers.Random_Complex_Circuit(nbr,dim,pwr);
    p : constant Standard_Complex_Polynomials.Poly
      := Standard_Circuit_Makers.Make_Polynomial(c1);
    c2 : constant Standard_Coefficient_Circuits.Circuit
       := Standard_Circuit_Makers.Split(c1);
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
    zd := Standard_Circuit_Makers.Gradient(p,cx);
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

    idx : constant Standard_Integer_Vectors.Vector
        := Standard_Circuit_Makers.Random_Indices(dim);
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
       := Standard_Circuit_Makers.Random_Complex_Circuit(nbr,dim);
    c2 : constant Standard_Coefficient_Circuits.Circuit
       := Standard_Circuit_Makers.Split(c1);
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
      Standard_Complex_Circuits.Indexed_Speel(c1,x,yd);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex Speel");
    tstart(timer);
    for k in 1..frq loop
      Standard_Coefficient_Circuits.Indexed_Speel(c2,xr,xi,ryd,iyd);
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
      rcf := Standard_Complex_Numbers.REAL_PART(cst);
      icf := Standard_Complex_Numbers.IMAG_PART(cst);
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
       := Standard_Circuit_Makers.Random_Complex_Circuit(nbr,dim,pwr);
    c2 : constant Standard_Coefficient_Circuits.Circuit
       := Standard_Circuit_Makers.Split(c1);
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

  function Merge_Jacobian_Matrix
             ( s : Standard_Coefficient_Circuits.Link_to_System )
             return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the columns of s.jrc and s.jic,
  --   merged into a complex matrix, for testing purposes.

    res : Standard_Complex_Matrices.Matrix(1..s.neq,1..s.dim);
    rlnk,ilnk : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for k in 1..s.dim loop -- make the k-th column of the matrix
      rlnk := s.jrc(k);
      ilnk := s.jic(k);
      for i in 1..s.neq loop -- run over all the row elements
        res(i,k) := Standard_Complex_Numbers.Create(rlnk(i),ilnk(i));
      end loop;
    end loop;
    return res;
  end Merge_Jacobian_Matrix;

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
    mrg : Standard_Complex_Matrices.Matrix(p'range,cx'range);
    err : double_float;

  begin
    Standard_Complex_Circuits.EvalDiff(cs,x);
    Standard_Coefficient_Circuits.EvalDiff(cffsys,xr,xi);
    put_line("The value at a random point :"); put_line(cs.fx);
    put_line("The value recomputed for testing :"); put_line(y);
    err := Evaluation_Differentiation_Errors.Sum_of_Errors(cs.fx,y);
    put("Sum of errors :"); put(err,3); new_line;
    put_line("The recomputed value :"); put_line(cffsys.fx);
    err := Evaluation_Differentiation_Errors.Sum_of_Errors(cffsys.fx,y);
    put("Sum of errors :"); put(err,3); new_line;
    Standard_Coefficient_Circuits.Eval(cffsys,xr,xi);
    put_line("The value computed without derivatives :");
    put_line(cffsys.fx);
    err := Evaluation_Differentiation_Errors.Sum_of_Errors(cffsys.fx,y);
    put("Sum of errors :"); put(err,3); new_line;
    put_line("The evaluated Jacobian matrix :");
    Standard_Circuit_Makers.Write_Matrix(cs.jm);
    put_line("The matrix recomputed for testing :"); 
    Standard_Circuit_Makers.Write_Matrix(jmx);
    err := Evaluation_Differentiation_Errors.Sum_of_Errors(cs.jm,jmx);
    put("Sum of errors :"); put(err,3); new_line;
    put_line("The recomputed matrix :"); 
    Standard_Circuit_Makers.Write_Matrix(cffsys.jm);
    err := Evaluation_Differentiation_Errors.Sum_of_Errors(cffsys.jm,jmx);
    put("Sum of errors :"); put(err,3); new_line;
    mrg := Merge_Jacobian_Matrix(cffsys);
    put_line("The merged recomputed matrix :");
    Standard_Circuit_Makers.Write_Matrix(mrg);
    err := Evaluation_Differentiation_Errors.Sum_of_Errors(mrg,jmx);
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
    cmpsys := Standard_Circuit_Makers.Make_Complex_System(lp);
    cffsys := Standard_Circuit_Makers.Split(cmpsys);
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
