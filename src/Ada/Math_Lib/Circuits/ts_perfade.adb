with text_io;                             use text_io;
with Communications_with_User;            use Communications_with_User;
with Timing_Package;                      use Timing_Package;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;        use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;         use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;         use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;         use Standard_Complex_VecVecs_io;
with Standard_Random_Vectors;
with Standard_Vector_Splitters;           use Standard_Vector_Splitters;
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
    print_times(standard_output,timer,"complex power table");
  end Timing_Power_Table;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension,
  --   the type of test, and then launches the test.

    dim,frq,pwr : integer32 := 0;
    ans,tst : character;

  begin
    new_line;
    put("Give the dimension of the vectors : "); get(dim);
    new_line;
    put_line("MENU for testing ADE :");
    put_line("  1. forward products");
    put_line("  2. forward and backward products");
    put_line("  3. forward, backward, and cross products");
    put_line("  4. indexed forward, backward, and cross products");
    put_line("  5. power table");
    put("Type 1, 2, 3, 4, or 5 to select the test : ");
    Ask_Alternative(tst,"12345");
    if tst = '5' then
      new_line;
      put("Give the highest power : "); get(pwr);
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
        when others => null;
      end case;
    end if;
  end Main;

begin
  Main;
end ts_perfade;
