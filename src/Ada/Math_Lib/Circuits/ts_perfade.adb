with text_io;                             use text_io;
with Communications_with_User;            use Communications_with_User;
with Timing_Package;                      use Timing_Package;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;         use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;
with Standard_Vector_Splitters;           use Standard_Vector_Splitters;

procedure ts_perfade is

-- DESCRIPTION :
--   Tests better performing algorithmic differentiation and evaluation.

  procedure Forward ( x : in Standard_Complex_Vectors.Link_to_Vector;
                      f : in Standard_Complex_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Computes all forward products of the values in x
  --   and stores the products in f.

  -- REQUIRED : f'first = x'first = 1 and f'last >= x'last-1.

    use Standard_Complex_Numbers;

  begin
    f(f'first) := x(x'first)*x(x'first+1);
    for k in 2..x'last-1 loop
      f(k) := f(k-1)*x(k+1);
    end loop;
  end Forward;

  procedure Forward ( xr : in Standard_Floating_Vectors.Link_to_Vector;
                      xi : in Standard_Floating_Vectors.Link_to_Vector;
                      fr : in Standard_Floating_Vectors.Link_to_Vector;
                      fi : in Standard_Floating_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Computes all forward products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in fr and the imaginary parts in fi.

  -- REQUIRED :
  --   xr'range = xi'range, fr'first = xr'first = 1,
  --   and fi'last >= xi'last-1.

    zr,zi,pr,pi,qr,qi : double_float;
    idx : integer32;
    dim : constant integer32 := xr'last-1;

  begin
   -- f(f'first) := x(x'first)*x(x'first+1);
    pr := xr(1); pi := xi(1);
    idx := xr'first+1;
    qr := xr(idx);  qi := xi(idx);
    zr := pr*qr - pi*qi;
    zi := pr*qi + pi*qr;
    fr(1) := zr; fi(1) := zi;
    for k in 2..dim loop 
     -- f(k) := f(k-1)*x(k+1);
      pr := zr; pi := zi;
      idx := k+1;
      qr := xr(idx); qi := xi(idx);
      zr := pr*qr - pi*qi;
      zi := pr*qi + pi*qr;
      fr(k) := zr; fi(k) := zi;
    end loop;
  end Forward;

  procedure Forward_Backward
              ( x,f,b : in Standard_Complex_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Computes all forward products of the values in x
  --   and stores the products in f.
  --   Computes all backward products of the values in x
  --   and stores the products in b.

  -- REQUIRED :
  --    f'first = x'first = 1 and f'last >= x'last-1.
  --    b'first = b'first = 1 and b'last >= x'last-2.

    use Standard_Complex_Numbers;

  begin
    f(f'first) := x(x'first)*x(x'first+1);
    for k in 2..x'last-1 loop
      f(k) := f(k-1)*x(k+1);
    end loop;
    b(b'first) := x(x'last)*x(x'last-1);
    for k in 2..x'last-2 loop
      b(k) := b(k-1)*x(x'last-k);
    end loop;
  end Forward_Backward;

  procedure Forward_Backward
              ( xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                fr : in Standard_Floating_Vectors.Link_to_Vector;
                fi : in Standard_Floating_Vectors.Link_to_Vector;
                br : in Standard_Floating_Vectors.Link_to_Vector;
                bi : in Standard_Floating_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Computes all forward products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in fr and the imaginary parts in fi.
  --   Computes all backward products of the values with real parts
  --   in xr and imaginary parts in xi and stores the real parts of
  --   the products in br and the imaginary parts in bi.

  -- REQUIRED :
  --   xr'range = xi'range,
  --   fr'first = xr'first = 1, fi'last >= xi'last-1,
  --   br'first = br'first = 1, bi'last >= bi'last-2.

    zr,zi,pr,pi,qr,qi : double_float;
    idx : integer32;
    dim : constant integer32 := xr'last-1;

  begin
   -- f(f'first) := x(x'first)*x(x'first+1);
    pr := xr(1); pi := xi(1);
    idx := xr'first+1;
    qr := xr(idx);  qi := xi(idx);
    zr := pr*qr - pi*qi;
    zi := pr*qi + pi*qr;
    fr(1) := zr; fi(1) := zi;
    for k in 2..dim loop 
     -- f(k) := f(k-1)*x(k+1);
      pr := zr; pi := zi;
      idx := k+1;
      qr := xr(idx); qi := xi(idx);
      zr := pr*qr - pi*qi;
      zi := pr*qi + pi*qr;
      fr(k) := zr; fi(k) := zi;
    end loop;
   -- b(b'first) := x(x'last)*x(x'last-1);
    pr := xr(xr'last); pi := xi(xr'last);
    idx := xi'last-1;
    qr := xr(idx); qi := xi(idx);
    zr := pr*qr - pi*qi;
    zi := pr*qi + pi*qr;
    br(1) := zr; bi(1) := zi;
    for k in 2..xr'last-2 loop
     -- b(k) := b(k-1)*x(x'last-k);
      pr := zr; pi := zi;
      idx := xr'last-k;
      qr := xr(idx); qi := xi(idx);
      zr := pr*qr - pi*qi;
      zi := pr*qi + pi*qr;
      br(k) := zr; bi(k) := zi;
    end loop;
  end Forward_Backward;

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

  begin
    Forward(x,f);
    put_line("the result : "); put_line(f);
    Forward(xr,xi,fr,fi);
    v := Make_Complex(fr,fi);
    put_line("recomputed : "); put_line(v);
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
    cb : constant Standard_Complex_Vectors.Vector(1..dim-2)
       := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    f : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cf);
    b : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cb);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    fr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(f);
    fi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(f);
    br : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(b);
    bi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(b);
    v,w : Standard_Complex_Vectors.Link_to_Vector;

  begin
    Forward_Backward(x,f,b);
    Forward_Backward(xr,xi,fr,fi,br,bi);
    v := Make_Complex(fr,fi);
    w := Make_Complex(br,bi);
    put_line("the forward products : "); put_line(f);
    put_line("forward products recomputed : "); put_line(v);
    put_line("the backward products : "); put_line(b);
    put_line("backward products recomputed : "); put_line(w);
  end Test_Forward_Backward;

  procedure Timing_Test_Forward ( dim,frq : in integer32 ) is

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
      Forward(x,f);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex forward products");
    tstart(timer);
    for k in 1..frq loop
      Forward(xr,xi,fr,fi);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"real forward products");
  end Timing_Test_Forward;

  procedure Timing_Test_Forward_Backward ( dim,frq : in integer32 ) is

  -- DESCRIPTION :
  --   Does as many forward/backward product computations as freq
  --   on random vectors of dimension dim.

    cx : constant Standard_Complex_Vectors.Vector(1..dim)
       := Standard_Random_Vectors.Random_Vector(1,dim);
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    cf : constant Standard_Complex_Vectors.Vector(1..dim-1)
       := Standard_Complex_Vectors.Vector'(1..dim-1 => zero);
    cb : constant Standard_Complex_Vectors.Vector(1..dim-2)
       := Standard_Complex_Vectors.Vector'(1..dim-2 => zero);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cx);
    f : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cf);
    b : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(cb);
    xr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(x);
    xi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(x);
    fr : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(f);
    fi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(f);
    br : constant Standard_Floating_Vectors.Link_to_Vector := Real_Part(b);
    bi : constant Standard_Floating_Vectors.Link_to_Vector := Imag_Part(b);
    timer : Timing_Widget;

  begin
    tstart(timer);
    for k in 1..frq loop
      Forward_Backward(x,f,b);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex forward & backward products");
    tstart(timer);
    for k in 1..frq loop
      Forward_Backward(xr,xi,fr,fi,br,bi);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"real forward & backward products");
  end Timing_Test_Forward_Backward;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension,
  --   the type of test, and then launches the test.

    dim,frq : integer32 := 0;
    ans,tst : character;

  begin
    new_line;
    put("Give the dimension of the vectors : "); get(dim);
    new_line;
    put_line("MENU for testing ADE :");
    put_line("  1. forward products");
    put_line("  2. forward and backward products");
    put("Type 1 or 2 to select the test : "); Ask_Alternative(tst,"12");
    new_line;
    put("Interactive tests ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      case tst is
        when '1' => Test_Forward(dim);
        when '2' => Test_Forward_Backward(dim);
        when others => null;
      end case;
    else
      new_line;
      put("Give the frequency of the tests : "); get(frq);
      case tst is
        when '1' => Timing_Test_Forward(dim,frq);
        when '2' => Timing_Test_Forward_Backward(dim,frq);
        when others => null;
      end case;
    end if;
  end Main;

begin
  Main;
end ts_perfade;
