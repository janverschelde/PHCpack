with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Functions;
with Standard_Random_Polynomials;        use Standard_Random_Polynomials;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Polynomials_io;    use DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Functions;
with DoblDobl_Random_Polynomials;        use DoblDobl_Random_Polynomials;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials_io;    use QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Poly_Functions;
with QuadDobl_Random_Polynomials;        use QuadDobl_Random_Polynomials;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Multprec_Complex_Solutions;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;        use Multprec_Complex_Vectors_io;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;    use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Functions;
with Multprec_Random_Polynomials;        use Multprec_Random_Polynomials;
with Standard_to_Multprec_Convertors;
with Standard_Extrapolators;
with DoblDobl_Extrapolators;
with QuadDobl_Extrapolators;
with Multprec_Extrapolators;
with Sample_Plane_Curves;                use Sample_Plane_Curves;

procedure ts_extrapol is

-- DESCRIPTION :
--   Development of numerical extrapolation for use in predictors.
--   We generate a plane algebraic curve defined by a random polynomial
--   in two variables, take samples at the curve and then extrapolate.
--   The higher the order of the extrapolator, the smaller the value
--   of the extrapolation when evaluated at the polynomial, which implies
--   that with higher-order extrapolators we stay closer to the curve.
--   The samples are taken "vertically" via new values for x.

  procedure Linear_Extrapolate
              ( p : in Standard_Complex_Polynomials.Poly;
                s : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies linear extrapolation to the first two points in the list s.

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    x,y,x0,x1,y0,y1,e : Standard_Complex_Numbers.Complex_Number;
    ls1,ls2 : Link_to_Solution;
    v : Standard_Complex_Vectors.Vector(1..2);
    abs_e : double_float;

  begin
    ls1 := Head_Of(s);
    ls2 := Head_Of(Tail_Of(s));
    x0 := ls1.v(1); y0 := ls1.v(2);
    x1 := ls2.v(1); y1 := ls2.v(2);
    x := Perturb(x1,0.1);
    y := Standard_Extrapolators.Extrapolate(x,x0,x1,y0,y1);
    v(1) := x; v(2) := y;
    e := Standard_Complex_Poly_Functions.Eval(p,v);
    abs_e := AbsVal(e);
    put_line("The extrapolated value at "); put_line(v);
    put("is "); put(e);
    put(" : "); put(abs_e);  new_line; 
  end Linear_Extrapolate;

  procedure Linear_Extrapolate
              ( p : in DoblDobl_Complex_Polynomials.Poly;
                s : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies linear extrapolation to the first two points in the list s.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    x,y,x0,x1,y0,y1,e : DoblDobl_Complex_Numbers.Complex_Number;
    ls1,ls2 : Link_to_Solution;
    v : DoblDobl_Complex_Vectors.Vector(1..2);
    abs_e : double_double;

  begin
    ls1 := Head_Of(s);
    ls2 := Head_Of(Tail_Of(s));
    x0 := ls1.v(1); y0 := ls1.v(2);
    x1 := ls2.v(1); y1 := ls2.v(2);
    x := Perturb(x1,0.1);
    y := DoblDobl_Extrapolators.Extrapolate(x,x0,x1,y0,y1);
    v(1) := x; v(2) := y;
    e := DoblDobl_Complex_Poly_Functions.Eval(p,v);
    abs_e := AbsVal(e);
    put_line("The extrapolated value at "); put_line(v);
    put("is "); put(e);
    put(" : "); put(abs_e);  new_line; 
  end Linear_Extrapolate;

  procedure Linear_Extrapolate
              ( p : in QuadDobl_Complex_Polynomials.Poly;
                s : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies linear extrapolation to the first two points in the list s.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    x,y,x0,x1,y0,y1,e : QuadDobl_Complex_Numbers.Complex_Number;
    ls1,ls2 : Link_to_Solution;
    v : QuadDobl_Complex_Vectors.Vector(1..2);
    abs_e : quad_double;

  begin
    ls1 := Head_Of(s);
    ls2 := Head_Of(Tail_Of(s));
    x0 := ls1.v(1); y0 := ls1.v(2);
    x1 := ls2.v(1); y1 := ls2.v(2);
    x := Perturb(x1,0.1);
    y := QuadDobl_Extrapolators.Extrapolate(x,x0,x1,y0,y1);
    v(1) := x; v(2) := y;
    e := QuadDobl_Complex_Poly_Functions.Eval(p,v);
    abs_e := AbsVal(e);
    put_line("The extrapolated value at "); put_line(v);
    put("is "); put(e);
    put(" : "); put(abs_e);  new_line; 
  end Linear_Extrapolate;

  procedure Linear_Extrapolate
              ( p : in Multprec_Complex_Polynomials.Poly;
                s : in Multprec_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies linear extrapolation to the first two points in the list s.

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Solutions;

    x,y,x0,x1,y0,y1,e : Multprec_Complex_Numbers.Complex_Number;
    ls1,ls2 : Link_to_Solution;
    v : Multprec_Complex_Vectors.Vector(1..2);
    abs_e : Floating_Number;

  begin
    ls1 := Head_Of(s);
    ls2 := Head_Of(Tail_Of(s));
    x0 := ls1.v(1); y0 := ls1.v(2);
    x1 := ls2.v(1); y1 := ls2.v(2);
    x := Perturb(x1,0.1);
    y := Multprec_Extrapolators.Extrapolate(x,x0,x1,y0,y1);
    v(1) := x; v(2) := y;
    e := Multprec_Complex_Poly_Functions.Eval(p,v);
    abs_e := AbsVal(e);
    put_line("The extrapolated value at "); put_line(v);
    put("is "); put(e);
    put(" : "); put(abs_e);  new_line; 
    Clear(x); Clear(y); Clear(e); Clear(abs_e);
  end Linear_Extrapolate;

  procedure Quadratic_Extrapolate
              ( p : in Standard_Complex_Polynomials.Poly;
                s : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies quadratic extrapolation to the first 3 points in the list s.

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    x,y,x0,x1,x2,y0,y1,y2,e : Standard_Complex_Numbers.Complex_Number;
    ls1,ls2,ls3 : Link_to_Solution;
    v : Standard_Complex_Vectors.Vector(1..2);
    abs_e : double_float;

  begin
    ls1 := Head_Of(s);
    ls2 := Head_Of(Tail_Of(s));
    ls3 := Head_Of(Tail_Of(Tail_Of(s)));
    x0 := ls1.v(1); y0 := ls1.v(2);
    x1 := ls2.v(1); y1 := ls2.v(2);
    x2 := ls3.v(1); y2 := ls3.v(2);
    x := Perturb(x2,0.1);
    y := Standard_Extrapolators.Extrapolate(x,x0,x1,x2,y0,y1,y2);
    v(1) := x; v(2) := y;
    e := Standard_Complex_Poly_Functions.Eval(p,v);
    abs_e := AbsVal(e);
    put_line("The extrapolated value at "); put_line(v);
    put("is "); put(e);
    put(" : "); put(abs_e);  new_line; 
  end Quadratic_Extrapolate;

  procedure Quadratic_Extrapolate
              ( p : in DoblDobl_Complex_Polynomials.Poly;
                s : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies quadratic extrapolation to the first 3 points in the list s.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    x,y,x0,x1,x2,y0,y1,y2,e : DoblDobl_Complex_Numbers.Complex_Number;
    ls1,ls2,ls3 : Link_to_Solution;
    v : DoblDobl_Complex_Vectors.Vector(1..2);
    abs_e : double_double;

  begin
    ls1 := Head_Of(s);
    ls2 := Head_Of(Tail_Of(s));
    ls3 := Head_Of(Tail_Of(Tail_Of(s)));
    x0 := ls1.v(1); y0 := ls1.v(2);
    x1 := ls2.v(1); y1 := ls2.v(2);
    x2 := ls3.v(1); y2 := ls3.v(2);
    x := Perturb(x2,0.1);
    y := DoblDobl_Extrapolators.Extrapolate(x,x0,x1,x2,y0,y1,y2);
    v(1) := x; v(2) := y;
    e := DoblDobl_Complex_Poly_Functions.Eval(p,v);
    abs_e := AbsVal(e);
    put_line("The extrapolated value at "); put_line(v);
    put("is "); put(e);
    put(" : "); put(abs_e);  new_line; 
  end Quadratic_Extrapolate;

  procedure Quadratic_Extrapolate
              ( p : in QuadDobl_Complex_Polynomials.Poly;
                s : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies quadratic extrapolation to the first 3 points in the list s.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    x,y,x0,x1,x2,y0,y1,y2,e : QuadDobl_Complex_Numbers.Complex_Number;
    ls1,ls2,ls3 : Link_to_Solution;
    v : QuadDobl_Complex_Vectors.Vector(1..2);
    abs_e : quad_double;

  begin
    ls1 := Head_Of(s);
    ls2 := Head_Of(Tail_Of(s));
    ls3 := Head_Of(Tail_Of(Tail_Of(s)));
    x0 := ls1.v(1); y0 := ls1.v(2);
    x1 := ls2.v(1); y1 := ls2.v(2);
    x2 := ls3.v(1); y2 := ls3.v(2);
    x := Perturb(x2,0.1);
    y := QuadDobl_Extrapolators.Extrapolate(x,x0,x1,x2,y0,y1,y2);
    v(1) := x; v(2) := y;
    e := QuadDobl_Complex_Poly_Functions.Eval(p,v);
    abs_e := AbsVal(e);
    put_line("The extrapolated value at "); put_line(v);
    put("is "); put(e);
    put(" : "); put(abs_e);  new_line; 
  end Quadratic_Extrapolate;

  procedure Quadratic_Extrapolate
              ( p : in Multprec_Complex_Polynomials.Poly;
                s : in Multprec_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies quadratic extrapolation to the first 3 points in the list s.

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Solutions;

    x,y,x0,x1,x2,y0,y1,y2,e : Multprec_Complex_Numbers.Complex_Number;
    ls1,ls2,ls3 : Link_to_Solution;
    v : Multprec_Complex_Vectors.Vector(1..2);
    abs_e : Floating_Number;

  begin
    ls1 := Head_Of(s);
    ls2 := Head_Of(Tail_Of(s));
    ls3 := Head_Of(Tail_Of(Tail_Of(s)));
    x0 := ls1.v(1); y0 := ls1.v(2);
    x1 := ls2.v(1); y1 := ls2.v(2);
    x2 := ls3.v(1); y2 := ls3.v(2);
    x := Perturb(x2,0.1);
    y := Multprec_Extrapolators.Extrapolate(x,x0,x1,x2,y0,y1,y2);
    v(1) := x; v(2) := y;
    e := Multprec_Complex_Poly_Functions.Eval(p,v);
    abs_e := AbsVal(e);
    put_line("The extrapolated value at "); put_line(v);
    put("is "); put(e);
    put(" : "); put(abs_e);  new_line; 
    Clear(x); Clear(y); Clear(e); Clear(abs_e);
  end Quadratic_Extrapolate;

  procedure Cubic_Extrapolate
              ( p : in Standard_Complex_Polynomials.Poly;
                s : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies cubic extrapolation to the first 4 points in the list s.

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    x,y,x0,x1,x2,x3,y0,y1,y2,y3,e : Standard_Complex_Numbers.Complex_Number;
    ls1,ls2,ls3,ls4 : Link_to_Solution;
    v : Standard_Complex_Vectors.Vector(1..2);
    tmp : Solution_List := s;
    abs_e : double_float;

  begin
    ls1 := Head_Of(s); tmp := Tail_Of(tmp);
    ls2 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls3 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls4 := Head_Of(tmp);
    x0 := ls1.v(1); y0 := ls1.v(2);
    x1 := ls2.v(1); y1 := ls2.v(2);
    x2 := ls3.v(1); y2 := ls3.v(2);
    x3 := ls4.v(1); y3 := ls4.v(2);
    x := Perturb(x3,0.1);
    y := Standard_Extrapolators.Extrapolate(x,x0,x1,x2,x3,y0,y1,y2,y3);
    v(1) := x; v(2) := y;
    e := Standard_Complex_Poly_Functions.Eval(p,v);
    abs_e := AbsVal(e);
    put_line("The extrapolated value at "); put_line(v);
    put("is "); put(e); 
    put(" : "); put(abs_e); new_line; 
  end Cubic_Extrapolate;

  procedure Cubic_Extrapolate
              ( p : in DoblDobl_Complex_Polynomials.Poly;
                s : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies cubic extrapolation to the first 4 points in the list s.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    x,y,x0,x1,x2,x3,y0,y1,y2,y3,e : DoblDobl_Complex_Numbers.Complex_Number;
    ls1,ls2,ls3,ls4 : Link_to_Solution;
    v : DoblDobl_Complex_Vectors.Vector(1..2);
    tmp : Solution_List := s;
    abs_e : double_double;

  begin
    ls1 := Head_Of(s); tmp := Tail_Of(tmp);
    ls2 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls3 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls4 := Head_Of(tmp);
    x0 := ls1.v(1); y0 := ls1.v(2);
    x1 := ls2.v(1); y1 := ls2.v(2);
    x2 := ls3.v(1); y2 := ls3.v(2);
    x3 := ls4.v(1); y3 := ls4.v(2);
    x := Perturb(x3,0.1);
    y := DoblDobl_Extrapolators.Extrapolate(x,x0,x1,x2,x3,y0,y1,y2,y3);
    v(1) := x; v(2) := y;
    e := DoblDobl_Complex_Poly_Functions.Eval(p,v);
    abs_e := AbsVal(e);
    put_line("The extrapolated value at "); put_line(v);
    put("is "); put(e); 
    put(" : "); put(abs_e); new_line; 
  end Cubic_Extrapolate;

  procedure Cubic_Extrapolate
              ( p : in QuadDobl_Complex_Polynomials.Poly;
                s : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies cubic extrapolation to the first 4 points in the list s.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    x,y,x0,x1,x2,x3,y0,y1,y2,y3,e : QuadDobl_Complex_Numbers.Complex_Number;
    ls1,ls2,ls3,ls4 : Link_to_Solution;
    v : QuadDobl_Complex_Vectors.Vector(1..2);
    tmp : Solution_List := s;
    abs_e : quad_double;

  begin
    ls1 := Head_Of(s); tmp := Tail_Of(tmp);
    ls2 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls3 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls4 := Head_Of(tmp);
    x0 := ls1.v(1); y0 := ls1.v(2);
    x1 := ls2.v(1); y1 := ls2.v(2);
    x2 := ls3.v(1); y2 := ls3.v(2);
    x3 := ls4.v(1); y3 := ls4.v(2);
    x := Perturb(x3,0.1);
    y := QuadDobl_Extrapolators.Extrapolate(x,x0,x1,x2,x3,y0,y1,y2,y3);
    v(1) := x; v(2) := y;
    e := QuadDobl_Complex_Poly_Functions.Eval(p,v);
    abs_e := AbsVal(e);
    put_line("The extrapolated value at "); put_line(v);
    put("is "); put(e); 
    put(" : "); put(abs_e); new_line; 
  end Cubic_Extrapolate;

  procedure Cubic_Extrapolate
              ( p : in Multprec_Complex_Polynomials.Poly;
                s : in Multprec_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies cubic extrapolation to the first 4 points in the list s.

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Solutions;

    x,y,x0,x1,x2,x3,y0,y1,y2,y3,e : Multprec_Complex_Numbers.Complex_Number;
    ls1,ls2,ls3,ls4 : Link_to_Solution;
    v : Multprec_Complex_Vectors.Vector(1..2);
    tmp : Solution_List := s;
    abs_e : Floating_Number;

  begin
    ls1 := Head_Of(s); tmp := Tail_Of(tmp);
    ls2 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls3 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls4 := Head_Of(tmp);
    x0 := ls1.v(1); y0 := ls1.v(2);
    x1 := ls2.v(1); y1 := ls2.v(2);
    x2 := ls3.v(1); y2 := ls3.v(2);
    x3 := ls4.v(1); y3 := ls4.v(2);
    x := Perturb(x3,0.1);
    y := Multprec_Extrapolators.Extrapolate(x,x0,x1,x2,x3,y0,y1,y2,y3);
    v(1) := x; v(2) := y;
    e := Multprec_Complex_Poly_Functions.Eval(p,v);
    abs_e := AbsVal(e);
    put_line("The extrapolated value at "); put_line(v);
    put("is "); put(e); 
    put(" : "); put(abs_e); new_line; 
  end Cubic_Extrapolate;

  procedure Quartic_Extrapolate
              ( p : in Standard_Complex_Polynomials.Poly;
                s : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies quartic extrapolation to the first 4 points in the list s.

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    x,y,e : Standard_Complex_Numbers.Complex_Number;
    x0,x1,x2,x3,x4,y0,y1,y2,y3,y4 : Standard_Complex_Numbers.Complex_Number;
    ls1,ls2,ls3,ls4,ls5 : Link_to_Solution;
    v : Standard_Complex_Vectors.Vector(1..2);
    tmp : Solution_List := s;
    abs_e : double_float;

  begin
    ls1 := Head_Of(s); tmp := Tail_Of(tmp);
    ls2 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls3 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls4 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls5 := Head_Of(tmp);
    x0 := ls1.v(1); y0 := ls1.v(2);
    x1 := ls2.v(1); y1 := ls2.v(2);
    x2 := ls3.v(1); y2 := ls3.v(2);
    x3 := ls4.v(1); y3 := ls4.v(2);
    x4 := ls5.v(1); y4 := ls5.v(2);
    x := Perturb(x4,0.1);
    y := Standard_Extrapolators.Extrapolate(x,x0,x1,x2,x3,x4,y0,y1,y2,y3,y4);
    v(1) := x; v(2) := y;
    e := Standard_Complex_Poly_Functions.Eval(p,v);
    abs_e := AbsVal(e);
    put_line("The extrapolated value at "); put_line(v);
    put("is "); put(e); 
    put(" : "); put(abs_e); new_line; 
  end Quartic_Extrapolate;

  procedure Quartic_Extrapolate
              ( p : in DoblDobl_Complex_Polynomials.Poly;
                s : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies quartic extrapolation to the first 4 points in the list s.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    x,y,e : DoblDobl_Complex_Numbers.Complex_Number;
    x0,x1,x2,x3,x4,y0,y1,y2,y3,y4 : DoblDobl_Complex_Numbers.Complex_Number;
    ls1,ls2,ls3,ls4,ls5 : Link_to_Solution;
    v : DoblDobl_Complex_Vectors.Vector(1..2);
    tmp : Solution_List := s;
    abs_e : double_double;

  begin
    ls1 := Head_Of(s); tmp := Tail_Of(tmp);
    ls2 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls3 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls4 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls5 := Head_Of(tmp);
    x0 := ls1.v(1); y0 := ls1.v(2);
    x1 := ls2.v(1); y1 := ls2.v(2);
    x2 := ls3.v(1); y2 := ls3.v(2);
    x3 := ls4.v(1); y3 := ls4.v(2);
    x4 := ls5.v(1); y4 := ls5.v(2);
    x := Perturb(x4,0.1);
    y := DoblDobl_Extrapolators.Extrapolate(x,x0,x1,x2,x3,x4,y0,y1,y2,y3,y4);
    v(1) := x; v(2) := y;
    e := DoblDobl_Complex_Poly_Functions.Eval(p,v);
    abs_e := AbsVal(e);
    put_line("The extrapolated value at "); put_line(v);
    put("is "); put(e); 
    put(" : "); put(abs_e); new_line; 
  end Quartic_Extrapolate;

  procedure Quartic_Extrapolate
              ( p : in QuadDobl_Complex_Polynomials.Poly;
                s : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies quartic extrapolation to the first 4 points in the list s.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    x,y,e : QuadDobl_Complex_Numbers.Complex_Number;
    x0,x1,x2,x3,x4,y0,y1,y2,y3,y4 : QuadDobl_Complex_Numbers.Complex_Number;
    ls1,ls2,ls3,ls4,ls5 : Link_to_Solution;
    v : QuadDobl_Complex_Vectors.Vector(1..2);
    tmp : Solution_List := s;
    abs_e : quad_double;

  begin
    ls1 := Head_Of(s); tmp := Tail_Of(tmp);
    ls2 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls3 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls4 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls5 := Head_Of(tmp);
    x0 := ls1.v(1); y0 := ls1.v(2);
    x1 := ls2.v(1); y1 := ls2.v(2);
    x2 := ls3.v(1); y2 := ls3.v(2);
    x3 := ls4.v(1); y3 := ls4.v(2);
    x4 := ls5.v(1); y4 := ls5.v(2);
    x := Perturb(x4,0.1);
    y := QuadDobl_Extrapolators.Extrapolate(x,x0,x1,x2,x3,x4,y0,y1,y2,y3,y4);
    v(1) := x; v(2) := y;
    e := QuadDobl_Complex_Poly_Functions.Eval(p,v);
    abs_e := AbsVal(e);
    put_line("The extrapolated value at "); put_line(v);
    put("is "); put(e); 
    put(" : "); put(abs_e); new_line; 
  end Quartic_Extrapolate;

  procedure Quartic_Extrapolate
              ( p : in Multprec_Complex_Polynomials.Poly;
                s : in Multprec_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies quartic extrapolation to the first 4 points in the list s.

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Solutions;

    x,y,e : Multprec_Complex_Numbers.Complex_Number;
    x0,x1,x2,x3,x4,y0,y1,y2,y3,y4 : Multprec_Complex_Numbers.Complex_Number;
    ls1,ls2,ls3,ls4,ls5 : Link_to_Solution;
    v : Multprec_Complex_Vectors.Vector(1..2);
    tmp : Solution_List := s;
    abs_e : Floating_Number;

  begin
    ls1 := Head_Of(s); tmp := Tail_Of(tmp);
    ls2 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls3 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls4 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls5 := Head_Of(tmp);
    x0 := ls1.v(1); y0 := ls1.v(2);
    x1 := ls2.v(1); y1 := ls2.v(2);
    x2 := ls3.v(1); y2 := ls3.v(2);
    x3 := ls4.v(1); y3 := ls4.v(2);
    x4 := ls5.v(1); y4 := ls5.v(2);
    x := Perturb(x4,0.1);
    y := Multprec_Extrapolators.Extrapolate(x,x0,x1,x2,x3,x4,y0,y1,y2,y3,y4);
    v(1) := x; v(2) := y;
    e := Multprec_Complex_Poly_Functions.Eval(p,v);
    abs_e := AbsVal(e);
    put_line("The extrapolated value at "); put_line(v);
    put("is "); put(e); 
    put(" : "); put(abs_e); new_line; 
  end Quartic_Extrapolate;

  procedure Quintic_Extrapolate
              ( p : in Standard_Complex_Polynomials.Poly;
                s : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies quintic extrapolation to the first 4 points in the list s.

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;
    use Standard_Extrapolators;

    x,y,x0,x1,x2,x3,x4,x5 : Standard_Complex_Numbers.Complex_Number;
    y0,y1,y2,y3,y4,y5,e : Standard_Complex_Numbers.Complex_Number;
    ls1,ls2,ls3,ls4,ls5,ls6 : Link_to_Solution;
    v : Standard_Complex_Vectors.Vector(1..2);
    tmp : Solution_List := s;
    abs_e : double_float;

  begin
    ls1 := Head_Of(s); tmp := Tail_Of(tmp);
    ls2 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls3 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls4 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls5 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls6 := Head_Of(tmp);
    x0 := ls1.v(1); y0 := ls1.v(2);
    x1 := ls2.v(1); y1 := ls2.v(2);
    x2 := ls3.v(1); y2 := ls3.v(2);
    x3 := ls4.v(1); y3 := ls4.v(2);
    x4 := ls5.v(1); y4 := ls5.v(2);
    x5 := ls6.v(1); y5 := ls6.v(2);
    x := Perturb(x5,0.1);
    y := Extrapolate(x,x0,x1,x2,x3,x4,x5,y0,y1,y2,y3,y4,y5);
    v(1) := x; v(2) := y;
    e := Standard_Complex_Poly_Functions.Eval(p,v);
    abs_e := AbsVal(e);
    put_line("The extrapolated value at "); put_line(v);
    put("is "); put(e); 
    put(" : "); put(abs_e); new_line; 
  end Quintic_Extrapolate;

  procedure Quintic_Extrapolate
              ( p : in DoblDobl_Complex_Polynomials.Poly;
                s : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies quintic extrapolation to the first 4 points in the list s.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Extrapolators;

    x,y,x0,x1,x2,x3,x4,x5 : DoblDobl_Complex_Numbers.Complex_Number;
    y0,y1,y2,y3,y4,y5,e : DoblDobl_Complex_Numbers.Complex_Number;
    ls1,ls2,ls3,ls4,ls5,ls6 : Link_to_Solution;
    v : DoblDobl_Complex_Vectors.Vector(1..2);
    tmp : Solution_List := s;
    abs_e : double_double;

  begin
    ls1 := Head_Of(s); tmp := Tail_Of(tmp);
    ls2 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls3 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls4 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls5 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls6 := Head_Of(tmp);
    x0 := ls1.v(1); y0 := ls1.v(2);
    x1 := ls2.v(1); y1 := ls2.v(2);
    x2 := ls3.v(1); y2 := ls3.v(2);
    x3 := ls4.v(1); y3 := ls4.v(2);
    x4 := ls5.v(1); y4 := ls5.v(2);
    x5 := ls6.v(1); y5 := ls6.v(2);
    x := Perturb(x5,0.1);
    y := Extrapolate(x,x0,x1,x2,x3,x4,x5,y0,y1,y2,y3,y4,y5);
    v(1) := x; v(2) := y;
    e := DoblDobl_Complex_Poly_Functions.Eval(p,v);
    abs_e := AbsVal(e);
    put_line("The extrapolated value at "); put_line(v);
    put("is "); put(e); 
    put(" : "); put(abs_e); new_line; 
  end Quintic_Extrapolate;

  procedure Quintic_Extrapolate
              ( p : in QuadDobl_Complex_Polynomials.Poly;
                s : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies quintic extrapolation to the first 4 points in the list s.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Extrapolators;

    x,y,x0,x1,x2,x3,x4,x5 : QuadDobl_Complex_Numbers.Complex_Number;
    y0,y1,y2,y3,y4,y5,e : QuadDobl_Complex_Numbers.Complex_Number;
    ls1,ls2,ls3,ls4,ls5,ls6 : Link_to_Solution;
    v : QuadDobl_Complex_Vectors.Vector(1..2);
    tmp : Solution_List := s;
    abs_e : quad_double;

  begin
    ls1 := Head_Of(s); tmp := Tail_Of(tmp);
    ls2 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls3 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls4 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls5 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls6 := Head_Of(tmp);
    x0 := ls1.v(1); y0 := ls1.v(2);
    x1 := ls2.v(1); y1 := ls2.v(2);
    x2 := ls3.v(1); y2 := ls3.v(2);
    x3 := ls4.v(1); y3 := ls4.v(2);
    x4 := ls5.v(1); y4 := ls5.v(2);
    x5 := ls6.v(1); y5 := ls6.v(2);
    x := Perturb(x5,0.1);
    y := Extrapolate(x,x0,x1,x2,x3,x4,x5,y0,y1,y2,y3,y4,y5);
    v(1) := x; v(2) := y;
    e := QuadDobl_Complex_Poly_Functions.Eval(p,v);
    abs_e := AbsVal(e);
    put_line("The extrapolated value at "); put_line(v);
    put("is "); put(e); 
    put(" : "); put(abs_e); new_line; 
  end Quintic_Extrapolate;

  procedure Quintic_Extrapolate
              ( p : in Multprec_Complex_Polynomials.Poly;
                s : in Multprec_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Applies quintic extrapolation to the first 4 points in the list s.

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Solutions;
    use Multprec_Extrapolators;

    x,y,x0,x1,x2,x3,x4,x5 : Multprec_Complex_Numbers.Complex_Number;
    y0,y1,y2,y3,y4,y5,e : Multprec_Complex_Numbers.Complex_Number;
    ls1,ls2,ls3,ls4,ls5,ls6 : Link_to_Solution;
    v : Multprec_Complex_Vectors.Vector(1..2);
    tmp : Solution_List := s;
    abs_e : Floating_Number;

  begin
    ls1 := Head_Of(s); tmp := Tail_Of(tmp);
    ls2 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls3 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls4 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls5 := Head_Of(tmp); tmp := Tail_Of(tmp);
    ls6 := Head_Of(tmp);
    x0 := ls1.v(1); y0 := ls1.v(2);
    x1 := ls2.v(1); y1 := ls2.v(2);
    x2 := ls3.v(1); y2 := ls3.v(2);
    x3 := ls4.v(1); y3 := ls4.v(2);
    x4 := ls5.v(1); y4 := ls5.v(2);
    x5 := ls6.v(1); y5 := ls6.v(2);
    x := Perturb(x5,0.1);
    y := Extrapolate(x,x0,x1,x2,x3,x4,x5,y0,y1,y2,y3,y4,y5);
    v(1) := x; v(2) := y;
    e := Multprec_Complex_Poly_Functions.Eval(p,v);
    abs_e := AbsVal(e);
    put_line("The extrapolated value at "); put_line(v);
    put("is "); put(e); 
    put(" : "); put(abs_e); new_line; 
  end Quintic_Extrapolate;

  procedure Standard_Slice_and_Extrapolate
              ( p : in Standard_Complex_Polynomials.Poly ) is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Solutions;

    d : constant integer32 := Degree(p);
    s : Array_of_Solution_Lists(0..d);
    b : Solution_List;

  begin
    put("A random polynomial of degree "); put(d,1); put(" :");
    put_line(p);
    Standard_Slice(p,natural32(d),s);
    for k in 1..d loop
      new_line;
      put("Extrapolating along branch "); put(k,1); put_line(" ...");
      b := Branch(s,natural32(k));
      put_line("-> linear extrapolation :");
      Linear_Extrapolate(p,b);
      if d > 1 then
        put_line("-> quadratic extrapolation :");
        Quadratic_Extrapolate(p,b);
      end if;
      if d > 2 then
        put_line("-> cubic extrapolation :");
        Cubic_Extrapolate(p,b);
      end if;
      if d > 3 then
        put_line("-> quartic extrapolation :");
        Quartic_Extrapolate(p,b);
      end if;
      if d > 4 then
        put_line("-> quintic extrapolation :");
        Quintic_Extrapolate(p,b);
      end if;
      Clear(b);
    end loop;
  end Standard_Slice_and_Extrapolate;

  procedure DoblDobl_Slice_and_Extrapolate
              ( p : in DoblDobl_Complex_Polynomials.Poly ) is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Solutions;

    d : constant integer32 := Degree(p);
    s : Array_of_Solution_Lists(0..d);
    b : Solution_List;

  begin
    put("A random polynomial of degree "); put(d,1); put(" :");
    put_line(p);
    DoblDobl_Slice(p,natural32(d),s);
    for k in 1..d loop
      new_line;
      put("Extrapolating along branch "); put(k,1); put_line(" ...");
      b := Branch(s,natural32(k));
      put_line("-> linear extrapolation :");
      Linear_Extrapolate(p,b);
      if d > 1 then
        put_line("-> quadratic extrapolation :");
        Quadratic_Extrapolate(p,b);
      end if;
      if d > 2 then
        put_line("-> cubic extrapolation :");
        Cubic_Extrapolate(p,b);
      end if;
      if d > 3 then
        put_line("-> quartic extrapolation :");
        Quartic_Extrapolate(p,b);
      end if;
      if d > 4 then
        put_line("-> quintic extrapolation :");
        Quintic_Extrapolate(p,b);
      end if;
    end loop;
  end DoblDobl_Slice_and_Extrapolate;

  procedure QuadDobl_Slice_and_Extrapolate
              ( p : in QuadDobl_Complex_Polynomials.Poly ) is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Solutions;

    d : constant integer32 := Degree(p);
    s : Array_of_Solution_Lists(0..d);
    b : Solution_List;

  begin
    put("A random polynomial of degree "); put(d,1); put(" :");
    put_line(p);
    QuadDobl_Slice(p,natural32(d),s);
    for k in 1..d loop
      new_line;
      put("Extrapolating along branch "); put(k,1); put_line(" ...");
      b := Branch(s,natural32(k));
      put_line("-> linear extrapolation :");
      Linear_Extrapolate(p,b);
      if d > 1 then
        put_line("-> quadratic extrapolation :");
        Quadratic_Extrapolate(p,b);
      end if;
      if d > 2 then
        put_line("-> cubic extrapolation :");
        Cubic_Extrapolate(p,b);
      end if;
      if d > 3 then
        put_line("-> quartic extrapolation :");
        Quartic_Extrapolate(p,b);
      end if;
      if d > 4 then
        put_line("-> quintic extrapolation :");
        Quintic_Extrapolate(p,b);
      end if;
    end loop;
  end QuadDobl_Slice_and_Extrapolate;

  procedure Multprec_Slice_and_Extrapolate
              ( p : in out Multprec_Complex_Polynomials.Poly ) is

    use Multprec_Complex_Polynomials;
    use Multprec_Complex_Solutions;

    d : constant integer32 := Degree(p);
    s : Array_of_Solution_Lists(0..d);
    deci : natural32 := 0;
    size : natural32;

  begin
    new_line;
    put("Give the number of decimal places : "); get(deci);
    size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
    Standard_to_Multprec_Convertors.Set_Size(p,size);
    put("A random polynomial of degree "); put(d,1); put_line(" :");
    put_line(p);
    Multprec_Slice(p,natural32(d),size,s);
   -- for k in 1..d loop
   --   declare
   --     b : Solution_List := Branch(s,natural32(k));
   --   begin
   --     put("The solutions in branch "); put(k,1); put_line(" :");
   --     Multprec_Evaluate(p,b);
   --   end;
   -- end loop;
    for k in 1..d loop
      new_line;
      put("Extrapolating along branch "); put(k,1); put_line(" ...");
      declare
        b : Solution_List := Branch(s,natural32(k));
      begin
        put_line("-> linear extrapolation :");
       -- put_line("The solutions in the branch before Linear Extrapolate :");
       -- Multprec_Evaluate(p,b);
        Linear_Extrapolate(p,b);
       -- new_line;
       -- put_line("The solutions in the branch after Linear_Extrapolate :");
       -- Multprec_Evaluate(p,b);
       -- new_line;
        if d > 1 then
          put_line("-> quadratic extrapolation :");
          Quadratic_Extrapolate(p,b);
        end if;
        if d > 2 then
          put_line("-> cubic extrapolation :");
          Cubic_Extrapolate(p,b);
        end if;
        if d > 3 then
          put_line("-> quartic extrapolation :");
          Quartic_Extrapolate(p,b);
        end if;
        if d > 4 then
          put_line("-> quintic extrapolation :");
          Quintic_Extrapolate(p,b);
        end if;
      end;
    end loop;
  end Multprec_Slice_and_Extrapolate;

  procedure Main is

    d : natural32 := 0;
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
      when '1' => st_p := Random_Dense_Poly(2,d,0);
                  Standard_Slice_and_Extrapolate(st_p);
      when '2' => dd_p := Random_Dense_Poly(2,d,0);
                  DoblDobl_Slice_and_Extrapolate(dd_p);
      when '3' => qd_p := Random_Dense_Poly(2,d,0);
                  QuadDobl_Slice_and_Extrapolate(qd_p);
      when '4' => mp_p := Random_Dense_Poly(2,d,0);
                  Multprec_Slice_and_Extrapolate(mp_p);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_extrapol;
