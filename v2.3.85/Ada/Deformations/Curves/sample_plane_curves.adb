with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Poly_Functions;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with DoblDobl_Random_Numbers;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Poly_Functions;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with QuadDobl_Random_Numbers;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_Poly_Functions;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;
with Multprec_Random_Numbers;
with Multprec_Complex_Number_Tools;
with Multprec_Complex_Vectors_io;        use Multprec_Complex_Vectors_io;
with Multprec_Complex_Vector_Tools;
with Standard_to_Multprec_Convertors;
with Multprec_Complex_Poly_Functions;
--with Multprec_Complex_Solutions_io;      use Multprec_Complex_Solutions_io; 
with Standard_Durand_Kerner;
with DoblDobl_Durand_Kerner;
with QuadDobl_Durand_Kerner;
with Multprec_Durand_Kerner;
with Black_Box_Univariate_Solvers;

package body Sample_Plane_Curves is

-- PART I : points on curves as solutions with same first coordinate

  function Insert_First_Coordinate
              ( s : Standard_Complex_Solutions.Solution;
                x : Standard_Complex_Numbers.Complex_Number )
              return Standard_Complex_Solutions.Solution is

    res : Standard_Complex_Solutions.Solution(s.n+1);

  begin
    res.t := s.t;
    res.m := s.m;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    res.v(1) := x;
    res.v(2) := s.v(1);
    return res;
  end Insert_First_Coordinate;

  function Insert_First_Coordinate
              ( s : DoblDobl_Complex_Solutions.Solution;
                x : DoblDobl_Complex_Numbers.Complex_Number )
              return DoblDobl_Complex_Solutions.Solution is

    res : DoblDobl_Complex_Solutions.Solution(s.n+1);

  begin
    res.t := s.t;
    res.m := s.m;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    res.v(1) := x;
    res.v(2) := s.v(1);
    return res;
  end Insert_First_Coordinate;

  function Insert_First_Coordinate
              ( s : QuadDobl_Complex_Solutions.Solution;
                x : QuadDobl_Complex_Numbers.Complex_Number )
              return QuadDobl_Complex_Solutions.Solution is

    res : QuadDobl_Complex_Solutions.Solution(s.n+1);

  begin
    res.t := s.t;
    res.m := s.m;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    res.v(1) := x;
    res.v(2) := s.v(1);
    return res;
  end Insert_First_Coordinate;

  function Insert_First_Coordinate
              ( s : Multprec_Complex_Solutions.Solution;
                x : Multprec_Complex_Numbers.Complex_Number )
              return Multprec_Complex_Solutions.Solution is

    res : Multprec_Complex_Solutions.Solution(s.n+1);

  begin
    Multprec_Complex_Numbers.Copy(s.t,res.t);
    res.m := s.m;
    Multprec_Floating_Numbers.Copy(s.err,res.err);
    Multprec_Floating_Numbers.Copy(s.rco,res.rco);
    Multprec_Floating_Numbers.Copy(s.res,res.res);
    Multprec_Complex_Numbers.Copy(x,res.v(1));
    Multprec_Complex_Numbers.Copy(s.v(1),res.v(2));
    return res;
  end Insert_First_Coordinate;

  function Insert_First_Coordinate
              ( s : Standard_Complex_Solutions.Solution_List;
                x : Standard_Complex_Numbers.Complex_Number )
              return Standard_Complex_Solutions.Solution_List is

    use Standard_Complex_Solutions;
    tmp : Standard_Complex_Solutions.Solution_List := s;
    res,res_last : Standard_Complex_Solutions.Solution_List;
    ls : Standard_Complex_Solutions.Link_to_Solution;

  begin
    while not Standard_Complex_Solutions.Is_Null(tmp) loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      Append(res,res_last,Insert_First_Coordinate(ls.all,x));
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Insert_First_Coordinate;

  function Insert_First_Coordinate
              ( s : DoblDobl_Complex_Solutions.Solution_List;
                x : DoblDobl_Complex_Numbers.Complex_Number )
              return DoblDobl_Complex_Solutions.Solution_List is

    use DoblDobl_Complex_Solutions;
    tmp : Solution_List := s;
    res,res_last : Solution_List;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Insert_First_Coordinate(ls.all,x));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Insert_First_Coordinate;

  function Insert_First_Coordinate
              ( s : QuadDobl_Complex_Solutions.Solution_List;
                x : QuadDobl_Complex_Numbers.Complex_Number )
              return QuadDobl_Complex_Solutions.Solution_List is

    use QuadDobl_Complex_Solutions;
    tmp : Solution_List := s;
    res,res_last : Solution_List;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Insert_First_Coordinate(ls.all,x));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Insert_First_Coordinate;

  function Insert_First_Coordinate
              ( s : Multprec_Complex_Solutions.Solution_List;
                x : Multprec_Complex_Numbers.Complex_Number )
              return Multprec_Complex_Solutions.Solution_List is

    use Multprec_Complex_Solutions;
    tmp : Solution_List := s;
    res,res_last : Solution_List;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Insert_First_Coordinate(ls.all,x));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Insert_First_Coordinate;

  function Second_Coordinate
             ( s : Standard_Complex_Solutions.Solution_List )
             return Standard_Complex_Vectors.Vector is

    use Standard_Complex_Solutions;
    res : Standard_Complex_Vectors.Vector(1..integer32(Length_Of(s)));
    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
    for i in res'range loop
      ls := Head_Of(tmp);
      res(i) := ls.v(2);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Second_Coordinate;

  function Second_Coordinate
             ( s : DoblDobl_Complex_Solutions.Solution_List )
             return DoblDobl_Complex_Vectors.Vector is

    use DoblDobl_Complex_Solutions;
    res : DoblDobl_Complex_Vectors.Vector(1..integer32(Length_Of(s)));
    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
    for i in res'range loop
      ls := Head_Of(tmp);
      res(i) := ls.v(2);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Second_Coordinate;

  function Second_Coordinate
             ( s : QuadDobl_Complex_Solutions.Solution_List )
             return QuadDobl_Complex_Vectors.Vector is

    use QuadDobl_Complex_Solutions;
    res : QuadDobl_Complex_Vectors.Vector(1..integer32(Length_Of(s)));
    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
    for i in res'range loop
      ls := Head_Of(tmp);
      res(i) := ls.v(2);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Second_Coordinate;

  function Second_Coordinate
             ( s : Multprec_Complex_Solutions.Solution_List )
             return Multprec_Complex_Vectors.Vector is

    use Multprec_Complex_Solutions;
    res : Multprec_Complex_Vectors.Vector(1..integer32(Length_Of(s)));
    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
    for i in res'range loop
      ls := Head_Of(tmp);
      Multprec_Complex_Numbers.Copy(ls.v(2),res(i));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Second_Coordinate;

  procedure Update_Samples
              ( s : in out Standard_Complex_Solutions.Solution_List;
                x : in Standard_Complex_Numbers.Complex_Number;
                y : in Standard_Complex_Vectors.Vector ) is

    use Standard_Complex_Solutions;
    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
    for i in y'range loop
      ls := Head_Of(tmp);
      ls.v(1) := x;
      ls.v(2) := y(i);
      tmp := Tail_Of(tmp);
    end loop;
  end Update_Samples;

  procedure Update_Samples
              ( s : in out DoblDobl_Complex_Solutions.Solution_List;
                x : in DoblDobl_Complex_Numbers.Complex_Number;
                y : in DoblDobl_Complex_Vectors.Vector ) is

    use DoblDobl_Complex_Solutions;
    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
    for i in y'range loop
      ls := Head_Of(tmp);
      ls.v(1) := x;
      ls.v(2) := y(i);
      tmp := Tail_Of(tmp);
    end loop;
  end Update_Samples;

  procedure Update_Samples
              ( s : in out QuadDobl_Complex_Solutions.Solution_List;
                x : in QuadDobl_Complex_Numbers.Complex_Number;
                y : in QuadDobl_Complex_Vectors.Vector ) is

    use QuadDobl_Complex_Solutions;
    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
    for i in y'range loop
      ls := Head_Of(tmp);
      ls.v(1) := x;
      ls.v(2) := y(i);
      tmp := Tail_Of(tmp);
    end loop;
  end Update_Samples;

  procedure Update_Samples
              ( s : in out Multprec_Complex_Solutions.Solution_List;
                x : in Multprec_Complex_Numbers.Complex_Number;
                y : in Multprec_Complex_Vectors.Vector ) is

    use Multprec_Complex_Solutions;
    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
   -- new_line;
   -- put_line("The solution list before the update : ");
   -- put(standard_output,Length_Of(s),natural32(Head_Of(s).n),s);
    for i in y'range loop
      ls := Head_Of(tmp);
      Multprec_Complex_Numbers.Copy(x,ls.v(1));
      Multprec_Complex_Numbers.Copy(y(i),ls.v(2));
      tmp := Tail_Of(tmp);
    end loop;
   -- new_line;
   -- put_line("The solution list after the update : ");
   -- put(standard_output,Length_Of(s),natural32(Head_Of(s).n),s);
  end Update_Samples;

-- PART II : calling univariate polynomial solvers

  function Sample ( p : Standard_Complex_Polynomials.Poly;
                    x : Standard_Complex_Numbers.Complex_Number )
                  return Standard_Complex_Solutions.Solution_List is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Solutions;
    use Black_Box_Univariate_Solvers;

    res,qsols : Solution_List;
    q : Poly := Standard_Complex_Poly_Functions.Eval(p,x,1);

  begin
    Black_Box_Durand_Kerner(q,qsols);
    res := Insert_First_Coordinate(qsols,x);
    Clear(q); Clear(qsols);
    return res;
  end Sample;

  function Sample ( p : DoblDobl_Complex_Polynomials.Poly;
                    x : DoblDobl_Complex_Numbers.Complex_Number )
                  return DoblDobl_Complex_Solutions.Solution_List is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Solutions;
    use Black_Box_Univariate_Solvers;

    res,qsols : Solution_List;
    q : Poly := DoblDobl_Complex_Poly_Functions.Eval(p,x,1);

  begin
    Black_Box_Durand_Kerner(q,qsols);
    res := Insert_First_Coordinate(qsols,x);
    Clear(q); Clear(qsols);
    return res;
  end Sample;

  function Sample ( p : QuadDobl_Complex_Polynomials.Poly;
                    x : QuadDobl_Complex_Numbers.Complex_Number )
                  return QuadDobl_Complex_Solutions.Solution_List is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Solutions;
    use Black_Box_Univariate_Solvers;

    res,qsols : Solution_List;
    q : Poly := QuadDobl_Complex_Poly_Functions.Eval(p,x,1);

  begin
    Black_Box_Durand_Kerner(q,qsols);
    res := Insert_First_Coordinate(qsols,x);
    Clear(q); Clear(qsols);
    return res;
  end Sample;

  function Sample ( p : Multprec_Complex_Polynomials.Poly;
                    x : Multprec_Complex_Numbers.Complex_Number;
                    size : natural32 )
                  return Multprec_Complex_Solutions.Solution_List is

    use Multprec_Complex_Polynomials;
    use Multprec_Complex_Solutions;
    use Black_Box_Univariate_Solvers;

    res,qsols : Solution_List;
    q : Poly := Multprec_Complex_Poly_Functions.Eval(p,x,1);

  begin
    Black_Box_Durand_Kerner(q,size,qsols);
    res := Insert_First_Coordinate(qsols,x);
    Clear(q); Clear(qsols);
    return res;
  end Sample;

  procedure Sample ( p : in Standard_Complex_Polynomials.Poly;
                     x : in Standard_Complex_Numbers.Complex_Number; 
                     sols : in out Standard_Complex_Solutions.Solution_List;
                     fail : out boolean ) is

    use Standard_Complex_Vectors;
    use Standard_Complex_Polynomials;
    use Standard_Complex_Solutions;
    use Standard_Durand_Kerner;
    use Black_Box_Univariate_Solvers;

    z : Vector(1..integer32(Length_Of(sols))) := Second_Coordinate(sols);
    r : Vector(z'range);
    q : Poly := Standard_Complex_Poly_Functions.Eval(p,x,1);
    d : constant integer32 := Degree(q);
    c : constant Vector := Coefficient_Vector(natural32(d),q);
    max : constant natural32 := 4;
    eps : constant double_float := 1.0E-12;
    nb : natural32;

  begin
    Silent_Durand_Kerner(c,z,r,max,eps,nb,fail);
    Update_Samples(sols,x,z);
    Clear(q);
  end Sample;

  procedure Sample ( p : in DoblDobl_Complex_Polynomials.Poly;
                     x : in DoblDobl_Complex_Numbers.Complex_Number; 
                     sols : in out DoblDobl_Complex_Solutions.Solution_List;
                     fail : out boolean ) is

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Durand_Kerner;
    use Black_Box_Univariate_Solvers;

    z : Vector(1..integer32(Length_Of(sols))) := Second_Coordinate(sols);
    r : Vector(z'range);
    q : Poly := DoblDobl_Complex_Poly_Functions.Eval(p,x,1);
    d : constant integer32 := Degree(q);
    c : constant Vector := Coefficient_Vector(natural32(d),q);
    max : constant natural32 := 8;
    eps : constant double_float := 1.0E-24;
    nb : natural32;

  begin
    DoblDobl_Durand_Kerner.Silent_Durand_Kerner(c,z,r,max,eps,nb,fail);
    Update_Samples(sols,x,z);
    Clear(q);
  end Sample;

  procedure Sample ( p : in QuadDobl_Complex_Polynomials.Poly;
                     x : in QuadDobl_Complex_Numbers.Complex_Number; 
                     sols : in out QuadDobl_Complex_Solutions.Solution_List;
                     fail : out boolean ) is

    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Durand_Kerner;
    use Black_Box_Univariate_Solvers;

    z : Vector(1..integer32(Length_Of(sols))) := Second_Coordinate(sols);
    r : Vector(z'range);
    q : Poly := QuadDobl_Complex_Poly_Functions.Eval(p,x,1);
    d : constant integer32 := Degree(q);
    c : constant Vector := Coefficient_Vector(natural32(d),q);
    max : constant natural32 := 16;
    eps : constant double_float := 1.0E-48;
    nb : natural32;

  begin
    QuadDobl_Durand_Kerner.Silent_Durand_Kerner(c,z,r,max,eps,nb,fail);
    Update_Samples(sols,x,z);
    Clear(q);
  end Sample;

  procedure Sample ( p : in Multprec_Complex_Polynomials.Poly;
                     x : in Multprec_Complex_Numbers.Complex_Number; 
                     size : in natural32;
                     sols : in out Multprec_Complex_Solutions.Solution_List;
                     fail : out boolean ) is

    use Multprec_Complex_Vectors;
    use Multprec_Complex_Polynomials;
    use Multprec_Complex_Solutions;
    use Multprec_Durand_Kerner;
    use Black_Box_Univariate_Solvers;

    z : Vector(1..integer32(Length_Of(sols))) := Second_Coordinate(sols);
    r : Vector(z'range);
    q : Poly := Multprec_Complex_Poly_Functions.Eval(p,x,1);
    d : constant integer32 := Degree(q);
    c : Vector(0..d) := Coefficient_Vector(natural32(d),q);
    max : constant natural32 := 10*size*natural32(d);
    deci : constant natural32
         := Multprec_Floating_Numbers.Size_to_Decimal(size);
    d_eps : constant double_float := 10.0**natural(-deci + 8);
    eps : constant Floating_Number := create(d_eps);
    nb : natural32;

   -- procedure Write ( step : in natural32; 
   --                   z,res : in Multprec_Complex_Vectors.Vector ) is
   -- begin
   --   put("step : "); put(step,1); put_line(" :");
   --   put("z : "); put_line(z);
   --   put("res : "); put_line(res);
   -- end Write;

   -- procedure RDK is
   --   new Multprec_Durand_Kerner.Reporting_Durand_Kerner(Write);

  begin
    Multprec_Complex_Vector_Tools.Set_Size(c,size);
    Multprec_Durand_Kerner.Silent_Durand_Kerner(c,z,r,max,eps,nb,fail);
   -- RDK(c,z,r,max,eps,nb,fail);
    Update_Samples(sols,x,z);
    Clear(q); Clear(z); Clear(r); Clear(c);
  end Sample;

  procedure Standard_Evaluate
              ( p : in Standard_Complex_Polynomials.Poly;
                s : in Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    tmp : Solution_List := s;
    ls : Link_to_Solution;
    sum : double_float := 0.0;
    y : Complex_Number;

  begin
    for i in 1..Length_Of(s) loop
      ls := Head_Of(tmp);
      put("Solution "); put(i,1); put_line(" :"); put_line(ls.v);
      y := Standard_Complex_Poly_Functions.Eval(p,ls.v);
      put("Value at solution "); put(i,1); put(" : "); put(y); new_line;
      sum := sum + Absval(y);
      tmp := Tail_Of(tmp);
    end loop;
    put("check sum : "); put(sum); new_line;
  end Standard_Evaluate;

  procedure DoblDobl_Evaluate
              ( p : in DoblDobl_Complex_Polynomials.Poly;
                s : in DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    tmp : Solution_List := s;
    ls : Link_to_Solution;
    sum : double_double := Double_Double_Numbers.create(0.0);
    y : Complex_Number;

  begin
    for i in 1..Length_Of(s) loop
      ls := Head_Of(tmp);
      put("Solution "); put(i,1); put_line(" :"); put_line(ls.v);
      y := DoblDobl_Complex_Poly_Functions.Eval(p,ls.v);
      put("Value at solution "); put(i,1); put(" : "); put(y); new_line;
      sum := sum + Absval(y);
      tmp := Tail_Of(tmp);
    end loop;
    put("check sum : "); put(sum); new_line;
  end DoblDobl_Evaluate;

  procedure QuadDobl_Evaluate
              ( p : in QuadDobl_Complex_Polynomials.Poly;
                s : in QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    tmp : Solution_List := s;
    ls : Link_to_Solution;
    sum : quad_double := Quad_Double_Numbers.create(0.0);
    y : Complex_Number;

  begin
    for i in 1..Length_Of(s) loop
      ls := Head_Of(tmp);
      put("Solution "); put(i,1); put_line(" :"); put_line(ls.v);
      y := QuadDobl_Complex_Poly_Functions.Eval(p,ls.v);
      put("Value at solution "); put(i,1); put(" : "); put(y); new_line;
      sum := sum + Absval(y);
      tmp := Tail_Of(tmp);
    end loop;
    put("check sum : "); put(sum); new_line;
  end QuadDobl_Evaluate;

  procedure Multprec_Evaluate
              ( p : in Multprec_Complex_Polynomials.Poly;
                s : in Multprec_Complex_Solutions.Solution_List ) is

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Solutions;

    tmp : Solution_List := s;
    ls : Link_to_Solution;
    sum : Floating_Number := Multprec_Floating_Numbers.create(0.0);
    y : Complex_Number;
    abs_y : Floating_Number;

  begin
    for i in 1..Length_Of(s) loop
      ls := Head_Of(tmp);
      put("Solution "); put(i,1); put_line(" :"); put_line(ls.v);
      y := Multprec_Complex_Poly_Functions.Eval(p,ls.v);
      put("Value at solution "); put(i,1); put(" : "); put(y); new_line;
      abs_y := AbsVal(y);
      Add(sum,abs_y);
      Clear(y); Clear(abs_y);
      tmp := Tail_Of(tmp);
    end loop;
    put("check sum : "); put(sum); new_line;
  end Multprec_Evaluate;

-- PART III : building a grid of points

  function Perturb ( x : Standard_Complex_Numbers.Complex_Number;
                     m : double_float )
                   return Standard_Complex_Numbers.Complex_Number is

    use Standard_Complex_Numbers;

    z : constant Complex_Number := Standard_Random_Numbers.Random1;
    res : constant Complex_Number := x + m*z;

  begin
    return res;
  end Perturb;

  function Perturb ( x : DoblDobl_Complex_Numbers.Complex_Number;
                     m : double_float )
                   return DoblDobl_Complex_Numbers.Complex_Number is

    use DoblDobl_Complex_Numbers;

    z : constant Complex_Number := DoblDobl_Random_Numbers.Random1;
    ddm : constant double_double := create(m);
    res : constant Complex_Number := x + ddm*z;

  begin
    return res;
  end Perturb;

  function Perturb ( x : QuadDobl_Complex_Numbers.Complex_Number;
                     m : double_float )
                   return QuadDobl_Complex_Numbers.Complex_Number is

    use QuadDobl_Complex_Numbers;

    z : constant Complex_Number := QuadDobl_Random_Numbers.Random1;
    qdm : constant quad_double := create(m);
    res : constant Complex_Number := x + qdm*z;

  begin
    return res;
  end Perturb;

  function Perturb ( x : Multprec_Complex_Numbers.Complex_Number;
                     m : double_float )
                   return Multprec_Complex_Numbers.Complex_Number is

    use Multprec_Complex_Numbers;

    z : Complex_Number := Multprec_Random_Numbers.Random(2);
    qdm : Floating_Number := Multprec_Floating_Numbers.create(m);
    res : Complex_Number;

  begin
    res := qdm*z;
    Add(res,x);
    Clear(z); Clear(qdm);
    return res;
  end Perturb;

  function Branch ( s : Standard_Complex_Solutions.Array_of_Solution_Lists;
                    k : natural32 )
                  return Standard_Complex_Solutions.Solution_List is

    use Standard_Complex_Solutions;
    res,res_last : Solution_List;

  begin
    for i in s'range loop
      Append(res,res_last,Retrieve(s(i),k));
    end loop;
    return res;
  end Branch;

  function Branch ( s : DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                    k : natural32 )
                  return DoblDobl_Complex_Solutions.Solution_List is

    use DoblDobl_Complex_Solutions;
    res,res_last : Solution_List;

  begin
    for i in s'range loop
      Append(res,res_last,Retrieve(s(i),k));
    end loop;
    return res;
  end Branch;

  function Branch ( s : QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                    k : natural32 )
                  return QuadDobl_Complex_Solutions.Solution_List is

    use QuadDobl_Complex_Solutions;
    res,res_last : Solution_List;

  begin
    for i in s'range loop
      Append(res,res_last,Retrieve(s(i),k));
    end loop;
    return res;
  end Branch;

  function Branch ( s : Multprec_Complex_Solutions.Array_of_Solution_Lists;
                    k : natural32 )
                  return Multprec_Complex_Solutions.Solution_List is

    use Multprec_Complex_Solutions;
    res,res_last : Solution_List;

  begin
    for i in s'range loop
      declare
        ls : Link_to_Solution := Retrieve(s(i),k);
        cp : Solution(ls.n);
      begin
        Copy(ls.all,cp);
        Append(res,res_last,cp);
      end;
    end loop;
    return res;
  end Branch;

  procedure Standard_Slice
              ( p : in Standard_Complex_Polynomials.Poly;
                d : in natural32;
                s : out Standard_Complex_Solutions.Array_of_Solution_Lists ) is

    x : Standard_Complex_Numbers.Complex_Number;
    fail : boolean;

  begin
    x := Standard_Random_Numbers.Random1;
    put("A random value for x1 : "); put(x); new_line;
    s(0) := Sample(p,x);
    Standard_Evaluate(p,s(0));
    for i in 1..integer32(d) loop
      x := Perturb(x,0.01);
      Standard_Complex_Solutions.Copy(s(i-1),s(i));
      Sample(p,x,s(i),fail);
      if fail
       then put_line("Failed to reach required accuracy!");
      end if;
      Standard_Evaluate(p,s(i));
    end loop;
  end Standard_Slice;

  procedure DoblDobl_Slice
              ( p : in DoblDobl_Complex_Polynomials.Poly;
                d : in natural32;
                s : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists ) is

    x : DoblDobl_Complex_Numbers.Complex_Number;
    fail : boolean;

  begin
    x := DoblDobl_Random_Numbers.Random1;
    put("A random value for x1 : "); put(x); new_line;
    s(0) := Sample(p,x);
    DoblDobl_Evaluate(p,s(0));
    for i in 1..integer32(d) loop
      x := Perturb(x,0.01);
      DoblDobl_Complex_Solutions.Copy(s(i-1),s(i));
      Sample(p,x,s(i),fail);
      if fail
       then put_line("Failed to reach required accuracy!");
      end if;
      DoblDobl_Evaluate(p,s(i));
    end loop;
  end DoblDobl_Slice;

  procedure QuadDobl_Slice
              ( p : in QuadDobl_Complex_Polynomials.Poly;
                d : in natural32;
                s : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists ) is

    x : QuadDobl_Complex_Numbers.Complex_Number;
    fail : boolean;

  begin
    x := QuadDobl_Random_Numbers.Random1;
    put("A random value for x1 : "); put(x); new_line;
    s(0) := Sample(p,x);
    QuadDobl_Evaluate(p,s(0));
    for i in 1..integer32(d) loop
      x := Perturb(x,0.01);
      QuadDobl_Complex_Solutions.Copy(s(i-1),s(i));
      Sample(p,x,s(i),fail);
      if fail
       then put_line("Failed to reach required accuracy!");
      end if;
      QuadDobl_Evaluate(p,s(i));
    end loop;
  end QuadDobl_Slice;

  function Multprec_Random1 return Multprec_Complex_Numbers.Complex_Number is

  -- DESCRPTION :
  --   Generates a random complex number in standard precision on the
  --   complex unit circle and converts it to multiprecision.

    res : Multprec_Complex_Numbers.Complex_Number;
    stcp : Standard_Complex_Numbers.Complex_Number;
    mpre,mpim : Floating_Number;

  begin
    stcp := Standard_Random_Numbers.Random1;
    mpre := Create(Standard_Complex_Numbers.REAL_PART(stcp));
    mpim := Create(Standard_Complex_Numbers.IMAG_PART(stcp));
    res := Multprec_Complex_Numbers.Create(mpre,mpim);
    Clear(mpre); Clear(mpim);
    return res;
  end Multprec_Random1;

  procedure Multprec_Slice
              ( p : in Multprec_Complex_Polynomials.Poly;
                d,size : in natural32;
                s : out Multprec_Complex_Solutions.Array_of_Solution_Lists ) is

    x : Multprec_Complex_Vectors.Vector(0..integer32(d));
    fail : boolean;

    use Multprec_Complex_Solutions;

  begin
    x(0) := Multprec_Random1;
    Multprec_Complex_Number_Tools.Set_Size(x(0),size);
    new_line;
    put("A random value for x(0) : "); put(x(0)); new_line;
    s(0) := Sample(p,x(0),size);
    Multprec_Evaluate(p,s(0));
    for i in 1..integer32(d) loop
      x(i) := Perturb(x(i-1),0.01);
      Multprec_Complex_Number_Tools.Set_Size(x(i),size);
      put("the value for x("); put(i,1); put(") : "); put(x(i)); new_line;
      Multprec_Complex_Solutions.Copy(s(i-1),s(i));
      Sample(p,x(i),size,s(i),fail);
     -- put("solution list "); put(i,1); put_line(" :");
     -- put(standard_output,Length_Of(s(i)),natural32(Head_Of(s(i)).n),s(i));
      if fail
       then put_line("Failed to reach required accuracy!");
      end if;
      Multprec_Evaluate(p,s(i));
    end loop;
   -- new_line;
   -- for i in s'range loop
   --   put("Solutions in list "); put(i,1); put_line(" :");
   --   put(standard_output,Length_Of(s(i)),natural32(Head_Of(s(i)).n),s(i));
   -- end loop;
   -- new_line;
  end Multprec_Slice;

end Sample_Plane_Curves;
