with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Functions;
with Standard_Random_Polynomials;        use Standard_Random_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_Root_Refiners;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with DoblDobl_Random_Numbers;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Polynomials_io;    use DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Functions;
with DoblDobl_Random_Polynomials;        use DoblDobl_Random_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Solutions;
with DoblDobl_Root_Refiners;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with QuadDobl_Random_Numbers;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials_io;    use QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Poly_Functions;
with QuadDobl_Random_Polynomials;        use QuadDobl_Random_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Solutions;
with QuadDobl_Root_Refiners;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;
with Multprec_Random_Numbers;
with Multprec_Complex_Vectors;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;    use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Functions;
with Multprec_Random_Polynomials;        use Multprec_Random_Polynomials;
with Standard_to_Multprec_Convertors;
with Multprec_Complex_Solutions;
with Multprec_Root_Refiners;
with Standard_Extrapolators;
with DoblDobl_Extrapolators;
with QuadDobl_Extrapolators;
with Multprec_Extrapolators;
with Sample_Plane_Curves;

procedure ts_extrapts is

-- DESCRIPTION :
--   Development of numerical extrapolation for use in predictors.
--   While the sampling of the points in ts_extrapol happened vertically,
--   in this test we cut the algebraic curve with a skew line.

  function Line ( a,b,c : Standard_Complex_Numbers.Complex_Number )
                return Standard_Complex_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns the line a*x + b*y + c as a polynomial in two variables. 

    use Standard_Complex_Polynomials;
    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..2 => 0);
    t.cf := c;
    res := Create(t);
    t.dg(1) := 1; t.cf := a;
    Add(res,t);
    t.dg(1) := 0; t.dg(2) := 1; t.cf := b;
    Add(res,t);
    Clear(t);
    return res;
  end Line;

  function Line ( a,b,c : DoblDobl_Complex_Numbers.Complex_Number )
                return DoblDobl_Complex_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns the line a*x + b*y + c as a polynomial in two variables. 

    use DoblDobl_Complex_Polynomials;
    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..2 => 0);
    t.cf := c;
    res := Create(t);
    t.dg(1) := 1; t.cf := a;
    Add(res,t);
    t.dg(1) := 0; t.dg(2) := 1; t.cf := b;
    Add(res,t);
    Clear(t);
    return res;
  end Line;

  function Line ( a,b,c : QuadDobl_Complex_Numbers.Complex_Number )
                return QuadDobl_Complex_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns the line a*x + b*y + c as a polynomial in two variables. 

    use QuadDobl_Complex_Polynomials;
    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..2 => 0);
    t.cf := c;
    res := Create(t);
    t.dg(1) := 1; t.cf := a;
    Add(res,t);
    t.dg(1) := 0; t.dg(2) := 1; t.cf := b;
    Add(res,t);
    Clear(t);
    return res;
  end Line;

  function Line ( a,b,c : Multprec_Complex_Numbers.Complex_Number )
                return Multprec_Complex_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns the line a*x + b*y + c as a polynomial in two variables. 

    use Multprec_Complex_Polynomials;
    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..2 => 0);
    Multprec_Complex_Numbers.Copy(c,t.cf);
    res := Create(t);
    t.dg(1) := 1;
    Multprec_Complex_Numbers.Copy(a,t.cf);
    Add(res,t);
    t.dg(1) := 0; t.dg(2) := 1;
    Multprec_Complex_Numbers.Copy(b,t.cf);
    Add(res,t);
    Clear(t);
    return res;
  end Line;

  procedure Standard_Newton
              ( f : in Standard_Complex_Poly_Systems.Poly_Sys;
                s : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies Newton's method with standard double arithmetic to the
  --   solutions in s of the system in f.

    use Standard_Root_Refiners;

    epsxa : constant double_float := 1.0E-12;
    epsfa : constant double_float := 1.0E-12;
    tolsing : constant double_float := 1.0E-8;
    nb : natural32 := 0;
    deflate : boolean := false;

  begin
    Silent_Root_Refiner(f,s,epsxa,epsfa,tolsing,nb,5,deflate);
  end Standard_Newton;

  procedure DoblDobl_Newton
              ( f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies Newton's method with standard double arithmetic to the
  --   solutions in s of the system in f.

    use DoblDobl_Root_Refiners;

    epsxa : constant double_float := 1.0E-24;
    epsfa : constant double_float := 1.0E-24;
    tolsing : constant double_float := 1.0E-8;
    nb : natural32 := 0;

  begin
    Silent_Root_Refiner(f,s,epsxa,epsfa,tolsing,nb,5);
  end DoblDobl_Newton;

  procedure QuadDobl_Newton
              ( f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies Newton's method with standard double arithmetic to the
  --   solutions in s of the system in f.

    use QuadDobl_Root_Refiners;

    epsxa : constant double_float := 1.0E-48;
    epsfa : constant double_float := 1.0E-48;
    tolsing : constant double_float := 1.0E-8;
    nb : natural32 := 0;

  begin
    Silent_Root_Refiner(f,s,epsxa,epsfa,tolsing,nb,9);
  end QuadDobl_Newton;

  procedure Multprec_Newton
              ( f : in Multprec_Complex_Poly_Systems.Poly_Sys;
                deci : in natural32;
                s : in out Multprec_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies Newton's method with multiprecsion arithmetic in the
  --   working precision defined by the number of decimal places in deci
  --   to the solutions in s of the system in f.

    use Multprec_Root_Refiners;

    d_eps : constant double_float := 10.0**natural(-deci + 8);
    d_tol : constant double_float := 10.0**natural(-deci/2);
    epsxa : Floating_Number := create(d_eps);
    epsfa : Floating_Number := create(d_eps);
    tolsing : Floating_Number := create(d_tol);
    nb : natural32 := 0;

  begin
    Silent_Root_Refiner(f,s,epsxa,epsfa,tolsing,nb,5+deci/5,false);
    Clear(epsxa); Clear(epsfa); Clear(tolsing);
  end Multprec_Newton;

  procedure Standard_Extrapolate
              ( p : in Standard_Complex_Polynomials.Poly;
                b : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies various extrapolators to the points in the list b.
    
    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;
    use Standard_Extrapolators;

    d : constant integer32 := integer32(Length_Of(b));
    s : Solution_Array(1..d) := Create(b);
    t : constant Complex_Number := s(d).t + Create(0.1);
    v : Standard_Complex_Vectors.Vector(1..2);
    y : Complex_Number;
    abs_y : double_float;

  begin
    v := Extrapolate(t,s(d-1).t,s(d).t,s(d-1).v,s(d).v);
    y := Standard_Complex_Poly_Functions.Eval(p,v);
    put("order 1 : "); put(y);
    abs_y := AbsVal(y); put(" : "); put(abs_y,3); new_line;
    if d > 2 then
      v := Extrapolate(t,s(d-2).t,s(d-1).t,s(d).t,s(d-2).v,s(d-1).v,s(d).v);
      y := Standard_Complex_Poly_Functions.Eval(p,v);
      put("order 2 : "); put(y);
      abs_y := AbsVal(y); put(" : "); put(abs_y,3); new_line;
    end if;
    if d > 3 then
      v := Extrapolate(t,s(d-3).t,s(d-2).t,s(d-1).t,s(d).t,
                         s(d-3).v,s(d-2).v,s(d-1).v,s(d).v);
      y := Standard_Complex_Poly_Functions.Eval(p,v);
      put("order 3 : "); put(y);
      abs_y := AbsVal(y); put(" : "); put(abs_y,3); new_line;
    end if;
    if d > 4 then
      v := Extrapolate(t,s(d-4).t,s(d-3).t,s(d-2).t,s(d-1).t,s(d).t,
                         s(d-4).v,s(d-3).v,s(d-2).v,s(d-1).v,s(d).v);
      y := Standard_Complex_Poly_Functions.Eval(p,v);
      put("order 4 : "); put(y);
      abs_y := AbsVal(y); put(" : "); put(abs_y,3); new_line;
    end if;
    if d > 5 then
      v := Extrapolate(t,s(d-5).t,s(d-4).t,s(d-3).t,s(d-2).t,s(d-1).t,s(d).t,
                         s(d-5).v,s(d-4).v,s(d-3).v,s(d-2).v,s(d-1).v,s(d).v);
      y := Standard_Complex_Poly_Functions.Eval(p,v);
      put("order 5 : "); put(y);
      abs_y := AbsVal(y); put(" : "); put(abs_y,3); new_line;
    end if;
    Clear(s);
  end Standard_Extrapolate;

  procedure DoblDobl_Extrapolate
              ( p : in DoblDobl_Complex_Polynomials.Poly;
                b : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies various extrapolators to the points in the list b.
    
    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Extrapolators;

    d : constant integer32 := integer32(Length_Of(b));
    s : Solution_Array(1..d) := Create(b);
    ddt : constant double_double := create(0.1);
    t : constant Complex_Number := s(d).t + Create(ddt);
    v : DoblDobl_Complex_Vectors.Vector(1..2);
    y : Complex_Number;
    abs_y : double_double;

  begin
    v := Extrapolate(t,s(d-1).t,s(d).t,s(d-1).v,s(d).v);
    y := DoblDobl_Complex_Poly_Functions.Eval(p,v);
    put("order 1 : "); put(y);
    abs_y := AbsVal(y); put(" : "); put(abs_y,3); new_line;
    if d > 2 then
      v := Extrapolate(t,s(d-2).t,s(d-1).t,s(d).t,s(d-2).v,s(d-1).v,s(d).v);
      y := DoblDobl_Complex_Poly_Functions.Eval(p,v);
      put("order 2 : "); put(y);
      abs_y := AbsVal(y); put(" : "); put(abs_y,3); new_line;
    end if;
    if d > 3 then
      v := Extrapolate(t,s(d-3).t,s(d-2).t,s(d-1).t,s(d).t,
                         s(d-3).v,s(d-2).v,s(d-1).v,s(d).v);
      y := DoblDobl_Complex_Poly_Functions.Eval(p,v);
      put("order 3 : "); put(y);
      abs_y := AbsVal(y); put(" : "); put(abs_y,3); new_line;
    end if;
    if d > 4 then
      v := Extrapolate(t,s(d-4).t,s(d-3).t,s(d-2).t,s(d-1).t,s(d).t,
                         s(d-4).v,s(d-3).v,s(d-2).v,s(d-1).v,s(d).v);
      y := DoblDobl_Complex_Poly_Functions.Eval(p,v);
      put("order 4 : "); put(y);
      abs_y := AbsVal(y); put(" : "); put(abs_y,3); new_line;
    end if;
    if d > 5 then
      v := Extrapolate(t,s(d-5).t,s(d-4).t,s(d-3).t,s(d-2).t,s(d-1).t,s(d).t,
                         s(d-5).v,s(d-4).v,s(d-3).v,s(d-2).v,s(d-1).v,s(d).v);
      y := DoblDobl_Complex_Poly_Functions.Eval(p,v);
      put("order 5 : "); put(y);
      abs_y := AbsVal(y); put(" : "); put(abs_y,3); new_line;
    end if;
    Clear(s);
  end DoblDobl_Extrapolate;

  procedure QuadDobl_Extrapolate
              ( p : in QuadDobl_Complex_Polynomials.Poly;
                b : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies various extrapolators to the points in the list b.
    
    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Extrapolators;

    d : constant integer32 := integer32(Length_Of(b));
    s : Solution_Array(1..d) := Create(b);
    ddt : constant quad_double := create(0.1);
    t : constant Complex_Number := s(d).t + Create(ddt);
    v : QuadDobl_Complex_Vectors.Vector(1..2);
    y : Complex_Number;
    abs_y : quad_double;

  begin
    v := Extrapolate(t,s(d-1).t,s(d).t,s(d-1).v,s(d).v);
    y := QuadDobl_Complex_Poly_Functions.Eval(p,v);
    put("order 1 : "); put(y);
    abs_y := AbsVal(y); put(" : "); put(abs_y,3); new_line;
    if d > 2 then
      v := Extrapolate(t,s(d-2).t,s(d-1).t,s(d).t,s(d-2).v,s(d-1).v,s(d).v);
      y := QuadDobl_Complex_Poly_Functions.Eval(p,v);
      put("order 2 : "); put(y);
      abs_y := AbsVal(y); put(" : "); put(abs_y,3); new_line;
    end if;
    if d > 3 then
      v := Extrapolate(t,s(d-3).t,s(d-2).t,s(d-1).t,s(d).t,
                         s(d-3).v,s(d-2).v,s(d-1).v,s(d).v);
      y := QuadDobl_Complex_Poly_Functions.Eval(p,v);
      put("order 3 : "); put(y);
      abs_y := AbsVal(y); put(" : "); put(abs_y,3); new_line;
    end if;
    if d > 4 then
      v := Extrapolate(t,s(d-4).t,s(d-3).t,s(d-2).t,s(d-1).t,s(d).t,
                         s(d-4).v,s(d-3).v,s(d-2).v,s(d-1).v,s(d).v);
      y := QuadDobl_Complex_Poly_Functions.Eval(p,v);
      put("order 4 : "); put(y);
      abs_y := AbsVal(y); put(" : "); put(abs_y,3); new_line;
    end if;
    if d > 5 then
      v := Extrapolate(t,s(d-5).t,s(d-4).t,s(d-3).t,s(d-2).t,s(d-1).t,s(d).t,
                         s(d-5).v,s(d-4).v,s(d-3).v,s(d-2).v,s(d-1).v,s(d).v);
      y := QuadDobl_Complex_Poly_Functions.Eval(p,v);
      put("order 5 : "); put(y);
      abs_y := AbsVal(y); put(" : "); put(abs_y,3); new_line;
    end if;
    Clear(s);
  end QuadDobl_Extrapolate;

  procedure Multprec_Extrapolate
              ( p : in Multprec_Complex_Polynomials.Poly;
                b : in Multprec_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies various extrapolators to the points in the list b.
    
    use Multprec_Complex_Numbers;
    use Multprec_Complex_Solutions;
    use Multprec_Extrapolators;

    d : constant integer32 := integer32(Length_Of(b));
    s : Solution_Array(1..d) := Create(b);
    ddt : Floating_Number := create(0.1);
    t : Complex_Number := s(d).t + Create(ddt);
    v : Multprec_Complex_Vectors.Vector(1..2);
    y : Complex_Number;
    abs_y : Floating_Number;

  begin
    v := Extrapolate(t,s(d-1).t,s(d).t,s(d-1).v,s(d).v);
    y := Multprec_Complex_Poly_Functions.Eval(p,v);
    put("order 1 : "); put(y);
    abs_y := AbsVal(y); put(" : "); put(abs_y,3); new_line;
    if d > 2 then
      v := Extrapolate(t,s(d-2).t,s(d-1).t,s(d).t,s(d-2).v,s(d-1).v,s(d).v);
      y := Multprec_Complex_Poly_Functions.Eval(p,v);
      put("order 2 : "); put(y);
      abs_y := AbsVal(y); put(" : "); put(abs_y,3); new_line;
    end if;
    if d > 3 then
      v := Extrapolate(t,s(d-3).t,s(d-2).t,s(d-1).t,s(d).t,
                         s(d-3).v,s(d-2).v,s(d-1).v,s(d).v);
      y := Multprec_Complex_Poly_Functions.Eval(p,v);
      put("order 3 : "); put(y);
      abs_y := AbsVal(y); put(" : "); put(abs_y,3); new_line;
    end if;
    if d > 4 then
      v := Extrapolate(t,s(d-4).t,s(d-3).t,s(d-2).t,s(d-1).t,s(d).t,
                         s(d-4).v,s(d-3).v,s(d-2).v,s(d-1).v,s(d).v);
      y := Multprec_Complex_Poly_Functions.Eval(p,v);
      put("order 4 : "); put(y);
      abs_y := AbsVal(y); put(" : "); put(abs_y,3); new_line;
    end if;
    if d > 5 then
      v := Extrapolate(t,s(d-5).t,s(d-4).t,s(d-3).t,s(d-2).t,s(d-1).t,s(d).t,
                         s(d-5).v,s(d-4).v,s(d-3).v,s(d-2).v,s(d-1).v,s(d).v);
      y := Multprec_Complex_Poly_Functions.Eval(p,v);
      put("order 5 : "); put(y);
      abs_y := AbsVal(y); put(" : "); put(abs_y,3); new_line;
    end if;
    Clear(s);
  end Multprec_Extrapolate;

  procedure Standard_Slice_and_Extrapolate
              ( p : in Standard_Complex_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Builds a grid of points starting at a vertical slice of
  --   the plane algebraic curve defined by the polynomial p.

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    d : constant integer32 := Degree(p);
    x : constant Complex_Number := Standard_Random_Numbers.Random1;
    s : Array_of_Solution_Lists(0..d);
    one : constant Complex_Number := Create(integer(1));
    t : Complex_Number := Create(0.01);
    h : Poly;
    f : Poly_Sys(1..2);
    b : Solution_List;

  begin
    put("A random polynomial of degree "); put(d,1); put(" :");
    put_line(p);
    s(0) := Sample_Plane_Curves.Sample(p,x);
    Sample_Plane_Curves.Standard_Evaluate(p,s(0));
    for k in 1..d loop
      h := Line(one,t,-x);
      put("The line : "); put_line(h);
      f(1) := p; f(2) := h;
      Copy(s(k-1),s(k));
      Set_Continuation_Parameter(s(k),t);
      Standard_Newton(f,s(k));
      Sample_Plane_Curves.Standard_Evaluate(p,s(k));
      t := t + Create(0.01);
    end loop;
    for k in 1..d loop
      b := Sample_Plane_Curves.Branch(s,natural32(k));
      new_line;
      put("extrapolation at branch "); put(k,1); put_line(" :");
      Standard_Extrapolate(p,b);
      Clear(b);
    end loop;
    Clear(s); 
  end Standard_Slice_and_Extrapolate;

  procedure DoblDobl_Slice_and_Extrapolate
              ( p : in DoblDobl_Complex_Polynomials.Poly ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    d : constant integer32 := Degree(p);
    x : constant Complex_Number := DoblDobl_Random_Numbers.Random1;
    s : Array_of_Solution_Lists(0..d);
    one : constant Complex_Number := Create(integer(1));
    ddt : constant double_double := create(0.01);
    t : Complex_Number := create(ddt);
    h : Poly;
    f : Poly_Sys(1..2);
    b : Solution_List;

  begin
    put("A random polynomial of degree "); put(d,1); put(" :");
    put_line(p);
    s(0) := Sample_Plane_Curves.Sample(p,x);
    Sample_Plane_Curves.DoblDobl_Evaluate(p,s(0));
    for k in 1..d loop
      h := Line(one,t,-x);
      put("The line : "); put_line(h);
      f(1) := p; f(2) := h;
      Copy(s(k-1),s(k));
      Set_Continuation_Parameter(s(k),t);
      DoblDobl_Newton(f,s(k));
      Sample_Plane_Curves.DoblDobl_Evaluate(p,s(k));
      t := t + Create(ddt);
    end loop;
    for k in 1..d loop
      b := Sample_Plane_Curves.Branch(s,natural32(k));
      new_line;
      put("extrapolation at branch "); put(k,1); put_line(" :");
      DoblDobl_Extrapolate(p,b);
      Clear(b);
    end loop;
    Clear(s);
  end DoblDobl_Slice_and_Extrapolate;

  procedure QuadDobl_Slice_and_Extrapolate
              ( p : in QuadDobl_Complex_Polynomials.Poly ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    d : constant integer32 := Degree(p);
    x : constant Complex_Number := QuadDobl_Random_Numbers.Random1;
    s : Array_of_Solution_Lists(0..d);
    one : constant Complex_Number := Create(integer(1));
    ddt : constant quad_double := create(0.01);
    t : Complex_Number := create(ddt);
    h : Poly;
    f : Poly_Sys(1..2);
    b : Solution_List;

  begin
    put("A random polynomial of degree "); put(d,1); put(" :");
    put_line(p);
    s(0) := Sample_Plane_Curves.Sample(p,x);
    Sample_Plane_Curves.QuadDobl_Evaluate(p,s(0));
    for k in 1..d loop
      h := Line(one,t,-x);
      put("The line : "); put_line(h);
      f(1) := p; f(2) := h;
      Copy(s(k-1),s(k));
      Set_Continuation_Parameter(s(k),t);
      QuadDobl_Newton(f,s(k));
      Sample_Plane_Curves.QuadDobl_Evaluate(p,s(k));
      t := t + Create(ddt);
    end loop;
    for k in 1..d loop
      b := Sample_Plane_Curves.Branch(s,natural32(k));
      new_line;
      put("extrapolation at branch "); put(k,1); put_line(" :");
      QuadDobl_Extrapolate(p,b);
      Clear(b);
    end loop;
    Clear(s);
  end QuadDobl_Slice_and_Extrapolate;

  procedure Multprec_Slice_and_Extrapolate
              ( p : in Multprec_Complex_Polynomials.Poly;
                deci,size : in natural32 ) is

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Polynomials;
    use Multprec_Complex_Poly_Systems;
    use Multprec_Complex_Solutions;

    d : constant integer32 := Degree(p);
    x : Complex_Number := Multprec_Random_Numbers.Random(size);
    s : Array_of_Solution_Lists(0..d);
    one : Complex_Number := Create(integer(1));
    ddt : Floating_Number := create(0.01);
    t : Complex_Number := create(ddt);
    h : Poly;
    f : Poly_Sys(1..2);
    b : Solution_List;

  begin
    put("A random polynomial of degree "); put(d,1); put(" :");
    put_line(p);
    s(0) := Sample_Plane_Curves.Sample(p,x,size);
    Sample_Plane_Curves.Multprec_Evaluate(p,s(0));
    for k in 1..d loop
      h := Line(one,t,-x);
      put("The line : "); put_line(h);
      f(1) := p; f(2) := h;
      Copy(s(k-1),s(k));
      Set_Continuation_Parameter(s(k),t);
      Multprec_Newton(f,deci,s(k));
      Sample_Plane_Curves.Multprec_Evaluate(p,s(k));
      t := t + Create(ddt);
    end loop;
    for k in 1..d loop
      b := Sample_Plane_Curves.Branch(s,natural32(k));
      new_line;
      put("extrapolation at branch "); put(k,1); put_line(" :");
      Multprec_Extrapolate(p,b);
      Clear(b);
    end loop;
    Clear(s);
  end Multprec_Slice_and_Extrapolate;

  procedure Main is

    d,deci,size : natural32 := 0;
    ans : character;
    st_p : Standard_Complex_Polynomials.Poly;
    dd_p : DoblDobl_Complex_Polynomials.Poly;
    qd_p : QuadDobl_Complex_Polynomials.Poly;
    mp_p : Multprec_Complex_Polynomials.Poly;

  begin
    new_line;
    put_line("Extrapolation on a plane algebraic curve ...");
    put("Give the degree : "); get(d);
    new_line;
    put_line("MENU to select the precision :");
    put_line("  1. execute in standard floating point arithmetic;");
    put_line("  2. execute in double double floating point arithmetic;");
    put_line("  3. execute in quad double floating point arithmetic;");
    put_line("  4. execute in multiprecision floating point arithmetic.");
    put("Type 1, 2, 3, or 4 to select the precision : ");
    Ask_Alternative(ans,"1234");
    case ans is
      when '1' =>
        st_p := Random_Dense_Poly(2,d,0);
        Standard_Slice_and_Extrapolate(st_p);
      when '2' =>
        dd_p := Random_Dense_Poly(2,d,0);
        DoblDobl_Slice_and_Extrapolate(dd_p);
      when '3' =>
        qd_p := Random_Dense_Poly(2,d,0);
        QuadDobl_Slice_and_Extrapolate(qd_p);
      when '4' =>
        new_line;
        put("Give the number of decimal places in the working precision : ");
        get(deci);
        size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
        mp_p := Random_Dense_Poly(2,d,0);
        Standard_to_Multprec_Convertors.Set_Size(mp_p,size);
        Multprec_Slice_and_Extrapolate(mp_p,deci,size);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_extrapts;
