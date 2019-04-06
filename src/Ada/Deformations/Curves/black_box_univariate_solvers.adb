with Timing_Package;                   use Timing_Package;
with Standard_Natural_Numbers_io;      use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;        use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;     use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;      use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers_io;      use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers_io;      use QuadDobl_Complex_Numbers_io;
with Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;     use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers_io;      use Multprec_Complex_Numbers_io;
with Double_Double_Numbers;
with Double_Double_Numbers_io;         use Double_Double_Numbers_io;
with Quad_Double_Numbers;
with Quad_Double_Numbers_io;           use Quad_Double_Numbers_io;
with Standard_Complex_Vectors_io;      use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors_io;      use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors_io;      use QuadDobl_Complex_Vectors_io;
with Multprec_Complex_Vectors_io;      use Multprec_Complex_Vectors_io;
with Standard_Random_Vectors;          use Standard_Random_Vectors;
with DoblDobl_Random_Vectors;          use DoblDobl_Random_Vectors;
with QuadDobl_Random_Vectors;          use QuadDobl_Random_Vectors;
with Multprec_Complex_Vector_Tools;
with Standard_Complex_Polynomials_io;  use Standard_Complex_Polynomials_io;
with DoblDobl_Complex_Polynomials_io;  use DoblDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Polynomials_io;  use QuadDobl_Complex_Polynomials_io;
with Standard_Complex_Laurentials_io;  use Standard_Complex_Laurentials_io;
with Multprec_Complex_Laurentials_io;  use Multprec_Complex_Laurentials_io;
with Standard_Complex_Solutions_io;    use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;    use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;    use QuadDobl_Complex_Solutions_io;
with Multprec_Complex_Solutions_io;    use Multprec_Complex_Solutions_io;
with Standard_Durand_Kerner;
with DoblDobl_Durand_Kerner;
with QuadDobl_Durand_Kerner;
with Multprec_Durand_Kerner;

package body Black_Box_Univariate_Solvers is

  function Create_Solution_List
             ( z : Standard_Complex_Numbers.Complex_Number )
             return Standard_Complex_Solutions.Solution_List is

    use Standard_Complex_Solutions;

    sols,last : Solution_List;
    sol : Solution(1);

  begin
    sol.t := Standard_Complex_Numbers.Create(0.0);
    sol.m := 1;
    sol.v(1) := z;
    sol.err := 0.0;
    sol.rco := 1.0;
    sol.res := 0.0;
    Append(sols,last,sol);
    return sols;
  end Create_Solution_List;

  function Create_Solution_List
             ( z : DoblDobl_Complex_Numbers.Complex_Number )
             return DoblDobl_Complex_Solutions.Solution_List is

    use DoblDobl_Complex_Solutions;

    sols,last : Solution_List;
    sol : Solution(1);

  begin
    sol.t := DoblDobl_Complex_Numbers.Create(integer(0));
    sol.m := 1;
    sol.v(1) := z;
    sol.err := Double_Double_Numbers.create(0.0);
    sol.rco := Double_Double_Numbers.create(1.0);
    sol.res := Double_Double_Numbers.create(0.0);
    Append(sols,last,sol);
    return sols;
  end Create_Solution_List;

  function Create_Solution_List
             ( z : QuadDobl_Complex_Numbers.Complex_Number )
             return QuadDobl_Complex_Solutions.Solution_List is

    use QuadDobl_Complex_Solutions;

    sols,last : Solution_List;
    sol : Solution(1);

  begin
    sol.t := QuadDobl_Complex_Numbers.Create(integer(0));
    sol.m := 1;
    sol.v(1) := z;
    sol.err := Quad_Double_Numbers.create(0.0);
    sol.rco := Quad_Double_Numbers.create(1.0);
    sol.res := Quad_Double_Numbers.create(0.0);
    Append(sols,last,sol);
    return sols;
  end Create_Solution_List;

  function Create_Solution_List
             ( z : Multprec_Complex_Numbers.Complex_Number )
             return Multprec_Complex_Solutions.Solution_List is

    use Multprec_Complex_Solutions;

    sols,last : Solution_List;
    sol : Solution(1);

  begin
    sol.t := Multprec_Complex_Numbers.Create(integer(0));
    sol.m := 1;
    Multprec_Complex_Numbers.Copy(z,sol.v(1));
    sol.err := Multprec_Floating_Numbers.create(0.0);
    sol.rco := Multprec_Floating_Numbers.create(1.0);
    sol.res := Multprec_Floating_Numbers.create(0.0);
    Append(sols,last,sol);
    return sols;
  end Create_Solution_List;

  function Create_Solution_List
               ( z : Standard_Complex_Vectors.Vector;
                 err,rco,res : Standard_Floating_Vectors.Vector )
               return Standard_Complex_Solutions.Solution_List is

    use Standard_Complex_Solutions;

    sols,last : Solution_List;
    sol : Solution(1);

  begin
    for i in z'range loop
      sol.t := Standard_Complex_Numbers.Create(0.0);
      sol.m := 1;
      sol.v(1) := z(i);
      sol.err := err(i);
      sol.rco := rco(i);
      sol.res := res(i);
      Append(sols,last,sol);
    end loop;
    return sols;
  end Create_Solution_List;

  function Create_Solution_List
               ( z : DoblDobl_Complex_Vectors.Vector;
                 err,rco,res : Standard_Floating_Vectors.Vector )
               return DoblDobl_Complex_Solutions.Solution_List is

    use DoblDobl_Complex_Solutions;

    sols,last : Solution_List;
    sol : Solution(1);

  begin
    for i in z'range loop
      sol.t := DoblDobl_Complex_Numbers.Create(integer(0));
      sol.m := 1;
      sol.v(1) := z(i);
      sol.err := Double_Double_Numbers.create(err(i));
      sol.rco := Double_Double_Numbers.create(rco(i));
      sol.res := Double_Double_Numbers.create(res(i));
      Append(sols,last,sol);
    end loop;
    return sols;
  end Create_Solution_List;

  function Create_Solution_List
               ( z : QuadDobl_Complex_Vectors.Vector;
                 err,rco,res : Standard_Floating_Vectors.Vector )
               return QuadDobl_Complex_Solutions.Solution_List is

    use QuadDobl_Complex_Solutions;

    sols,last : Solution_List;
    sol : Solution(1);

  begin
    for i in z'range loop
      sol.t := QuadDobl_Complex_Numbers.Create(integer(0));
      sol.m := 1;
      sol.v(1) := z(i);
      sol.err := Quad_Double_Numbers.create(err(i));
      sol.rco := Quad_Double_Numbers.create(rco(i));
      sol.res := Quad_Double_Numbers.create(res(i));
      Append(sols,last,sol);
    end loop;
    return sols;
  end Create_Solution_List;

  function Create_Solution_List
               ( z : Multprec_Complex_Vectors.Vector;
                 err,rco,res : Multprec_Floating_Vectors.Vector )
               return Multprec_Complex_Solutions.Solution_List is

    use Multprec_Complex_Solutions;

    sols,last : Solution_List;
    sol : Solution(1);

  begin
    for i in z'range loop
      sol.t := Multprec_Complex_Numbers.Create(integer(0));
      sol.m := 1;
      Multprec_Complex_Numbers.Copy(z(i),sol.v(1));
      Multprec_Floating_Numbers.Copy(err(i),sol.err);
      Multprec_Floating_Numbers.Copy(rco(i),sol.rco);
      Multprec_Floating_Numbers.Copy(res(i),sol.res);
      Append(sols,last,sol);
    end loop;
    return sols;
  end Create_Solution_List;

  procedure Write_Results
               ( file : in file_type; step : in natural32;
                 z,res : in Standard_Complex_Vectors.Vector ) is

    absres : double_float;

  begin
    new_line(file);
    put(file,"Results after "); put(file,step,1);
    put_line(file," iterations :");
    put_line(file,
    "------------------------------------------------------------------------");
    put_line(file,
    "| APPROXIMATED ROOTS (real and imaginary part) |     RESIDUALS         |");
    put_line(file,
    "------------------------------------------------------------------------");
    for i in z'range loop
      absres := Standard_Complex_Numbers.AbsVal(res(i));
      put(file,"| "); put(file,z(i)); put(file," | "); 
      put(file,absres,3); put(file," |");
      new_line(file);
    end loop;
    put_line(file,
    "------------------------------------------------------------------------");
  end Write_Results;

  procedure Write_Results
               ( file : in file_type; step : in natural32;
                 z,res : in DoblDobl_Complex_Vectors.Vector ) is

    use Double_Double_Numbers;
    absres : double_double;

  begin
    new_line(file);
    put(file,"Results after "); put(file,step,1);
    put_line(file," iterations :");
    put_line(file,
    "------------------------------------------------------------------------");
    put_line(file,
    "| APPROXIMATED ROOTS (real and imaginary part) |     RESIDUALS         |");
    put_line(file,
    "------------------------------------------------------------------------");
    for i in z'range loop
      absres := DoblDobl_Complex_Numbers.AbsVal(res(i));
      put(file,"| "); put(file,z(i)); put(file," | "); 
      put(file,absres,3); put(file," |");
      new_line(file);
    end loop;
    put_line(file,
    "------------------------------------------------------------------------");
  end Write_Results;

  procedure Write_Results
               ( file : in file_type; step : in natural32;
                 z,res : in QuadDobl_Complex_Vectors.Vector ) is

    use Quad_Double_Numbers;
    absres : quad_double;

  begin
    new_line(file);
    put(file,"Results after "); put(file,step,1);
    put_line(file," iterations :");
    put_line(file,
    "------------------------------------------------------------------------");
    put_line(file,
    "| APPROXIMATED ROOTS (real and imaginary part) |     RESIDUALS         |");
    put_line(file,
    "------------------------------------------------------------------------");
    for i in z'range loop
      absres := QuadDobl_Complex_Numbers.AbsVal(res(i));
      put(file,"| "); put(file,z(i)); put(file," | "); 
      put(file,absres,3); put(file," |");
      new_line(file);
    end loop;
    put_line(file,
    "------------------------------------------------------------------------");
  end Write_Results;

  procedure Write_Results
               ( file : in file_type; step : in natural32;
                 z,res : in Multprec_Complex_Vectors.Vector ) is

    use Multprec_Floating_Numbers;
    absres : Floating_Number;

  begin
    new_line(file);
    put(file,"Results after "); put(file,step,1);
    put_line(file," iterations :");
    put_line(file,
    "------------------------------------------------------------------------");
    put_line(file,
    "| APPROXIMATED ROOTS (real and imaginary part) |     RESIDUALS         |");
    put_line(file,
    "------------------------------------------------------------------------");
    for i in z'range loop
      absres := Multprec_Complex_Numbers.AbsVal(res(i));
      put(file,"| "); put(file,z(i)); put(file," | "); 
      put(file,absres,3); put(file," |");
      new_line(file);
      Clear(absres);
    end loop;
    put_line(file,
    "------------------------------------------------------------------------");
  end Write_Results;

  procedure Standard_Compute_Roots
              ( p : in Standard_Complex_Vectors.Vector;
                z : out Standard_Complex_Vectors.Vector;
                err,rco,res : out Standard_Floating_Vectors.Vector;
                fail : out boolean ) is

    use Standard_Durand_Kerner;

    deg : constant integer32 := p'last;
    max : constant natural32 := 10*natural32(deg);
    eps : constant double_float := 1.0E-12;
    r : Standard_Complex_Vectors.Vector(1..deg);
    dp : constant Standard_Complex_Vectors.Vector(0..deg-1) := Derivative(p);
    nb : natural32 := 0;

  begin
    z := Random_Vector(1,deg);
    Silent_Durand_Kerner(p,z,r,max,eps,nb,fail);
    Newton(p,dp,z,err,rco,res);
  end Standard_Compute_Roots;

  procedure DoblDobl_Compute_Roots
              ( p : in DoblDobl_Complex_Vectors.Vector;
                z : out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out Standard_Floating_Vectors.Vector;
                fail : out boolean ) is

    use DoblDobl_Durand_Kerner;
  
    deg : constant integer32 := p'last;
    max : constant natural32 := 10*natural32(deg);
    eps : constant double_float := 1.0E-24;
    r : DoblDobl_Complex_Vectors.Vector(1..deg);
    dp : constant DoblDobl_Complex_Vectors.Vector(0..deg-1) := Derivative(p);
    nb : natural32 := 0;

  begin
    z := Random_Vector(1,deg);
    Silent_Durand_Kerner(p,z,r,max,eps,nb,fail);
    Newton(p,dp,z,err,rco,res);
  end DoblDobl_Compute_Roots;

  procedure QuadDobl_Compute_Roots
              ( p : in QuadDobl_Complex_Vectors.Vector;
                z : out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out Standard_Floating_Vectors.Vector;
                fail : out boolean ) is

    use QuadDobl_Durand_Kerner;
  
    deg : constant integer32 := p'last;
    max : constant natural32 := 10*natural32(deg);
    eps : constant double_float := 1.0E-48;
    r : QuadDobl_Complex_Vectors.Vector(1..deg);
    dp : constant QuadDobl_Complex_Vectors.Vector(0..deg-1) := Derivative(p);
    nb : natural32 := 0;

  begin
    z := Random_Vector(1,deg);
    Silent_Durand_Kerner(p,z,r,max,eps,nb,fail);
    Newton(p,dp,z,err,rco,res);
  end QuadDobl_Compute_Roots;

  procedure Standard_Find_Roots
              ( d : in integer32;
                cp : in Standard_Complex_Vectors.Vector;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    z : Standard_Complex_Vectors.Vector(1..d);
    err,rco,res : Standard_Floating_Vectors.Vector(z'range);
    fail : boolean;

  begin
    Standard_Compute_Roots(cp,z,err,rco,res,fail);
    sols := Create_Solution_List(z,err,rco,res);
  end Standard_Find_Roots;

  procedure DoblDobl_Find_Roots
              ( d : in integer32;
                cp : in DoblDobl_Complex_Vectors.Vector;
                sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    z : DoblDobl_Complex_Vectors.Vector(1..d);
    err,rco,res : Standard_Floating_Vectors.Vector(z'range);
    fail : boolean;

  begin
    DoblDobl_Compute_Roots(cp,z,err,rco,res,fail);
    sols := Create_Solution_List(z,err,rco,res);
  end DoblDobl_Find_Roots;

  procedure QuadDobl_Find_Roots
              ( d : in integer32;
                cp : in QuadDobl_Complex_Vectors.Vector;
                sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    z : QuadDobl_Complex_Vectors.Vector(1..d);
    err,rco,res : Standard_Floating_Vectors.Vector(z'range);
    fail : boolean;

  begin
    QuadDobl_Compute_Roots(cp,z,err,rco,res,fail);
    sols := Create_Solution_List(z,err,rco,res);
  end QuadDobl_Find_Roots;

  function Random_Vector1
              ( first,last : integer32 ) 
              return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(first..last);
    stc : constant Standard_Complex_Vectors.Vector(first..last)
        := Random_Vector(first,last);
    mpre,mpim : Multprec_Floating_Numbers.Floating_Number;

  begin
    for i in res'range loop
      mpre := Multprec_Floating_Numbers.Create
                (Standard_Complex_Numbers.REAL_PART(stc(i)));
      mpim := Multprec_Floating_Numbers.Create
                (Standard_Complex_Numbers.IMAG_PART(stc(i)));
      res(i) := Multprec_Complex_Numbers.Create(mpre,mpim);
      Multprec_Floating_Numbers.Clear(mpre);
      Multprec_Floating_Numbers.Clear(mpim);
    end loop;
    return res;
  end Random_Vector1;

  procedure Multprec_Find_Roots
              ( d : in integer32; size : in natural32;
                cp : in Multprec_Complex_Vectors.Vector;
                sols : out Multprec_Complex_Solutions.Solution_List ) is

    use Multprec_Durand_Kerner;
  
    max : constant natural32 := 10*size*natural32(d);
    deci : constant natural32
         := Multprec_Floating_Numbers.Size_to_Decimal(size);
    d_eps : constant double_float := 10.0**natural(-deci + 8);
    eps : Multprec_Floating_Numbers.Floating_Number
        := Multprec_Floating_Numbers.create(d_eps);
    z,r : Multprec_Complex_Vectors.Vector(1..d);
    dcp : Multprec_Complex_Vectors.Vector(0..d-1) := Derivative(cp);
    err,rco,res : Multprec_Floating_Vectors.Vector(z'range);
    nb : natural32 := 0;
    fail : boolean;

  begin
    z := Random_Vector1(1,d);
    Multprec_Complex_Vector_Tools.Set_Size(z,size);
    Multprec_Complex_Vector_Tools.Set_Size(dcp,size);
    Silent_Durand_Kerner(cp,z,r,max,eps,nb,fail);
    Newton(cp,dcp,z,err,rco,res);
    sols := Create_Solution_List(z,err,rco,res);
    Multprec_Floating_Numbers.Clear(eps);
    Multprec_Complex_Vectors.Clear(z);
    Multprec_Complex_Vectors.Clear(r);
  end Multprec_Find_Roots;

  procedure Standard_Find_Roots
              ( file : in file_type; d : in integer32;
                cp : in Standard_Complex_Vectors.Vector;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Durand_Kerner;

    max : constant natural32 := 10*natural32(d);
    eps : constant double_float := 1.0E-13;
    z,r : Standard_Complex_Vectors.Vector(1..d);
    dcp : constant Standard_Complex_Vectors.Vector(0..d-1) := Derivative(cp);
    err,rco,res : Standard_Floating_Vectors.Vector(z'range);
    nb : natural32;
    fail : boolean;

  begin
    z := Random_Vector(1,d);
    Silent_Durand_Kerner(cp,z,r,max,eps,nb,fail);
    if fail
     then put_line(file,"precision insufficient to reach results");
    end if;
    Write_Results(file,nb,z,r);
    Newton(cp,dcp,z,err,rco,res);
    sols := Create_Solution_List(z,err,rco,res);
  end Standard_Find_Roots;

  procedure DoblDobl_Find_Roots
              ( file : in file_type; d : in integer32;
                cp : in DoblDobl_Complex_Vectors.Vector;
                sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Durand_Kerner;

    max : constant natural32 := 10*natural32(d);
    eps : constant double_float := 1.0E-13;
    z,r : DoblDobl_Complex_Vectors.Vector(1..d);
    dcp : constant DoblDobl_Complex_Vectors.Vector(0..d-1) := Derivative(cp);
    err,rco,res : Standard_Floating_Vectors.Vector(z'range);
    nb : natural32;
    fail : boolean;

  begin
    z := Random_Vector(1,d);
    Silent_Durand_Kerner(cp,z,r,max,eps,nb,fail);
    if fail
     then put_line(file,"precision insufficient to reach results");
    end if;
    Write_Results(file,nb,z,r);
    Newton(cp,dcp,z,err,rco,res);
    sols := Create_Solution_List(z,err,rco,res);
  end DoblDobl_Find_Roots;

  procedure QuadDobl_Find_Roots
              ( file : in file_type; d : in integer32;
                cp : in QuadDobl_Complex_Vectors.Vector;
                sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Durand_Kerner;

    max : constant natural32 := 10*natural32(d);
    eps : constant double_float := 1.0E-13;
    z,r : QuadDobl_Complex_Vectors.Vector(1..d);
    dcp : constant QuadDobl_Complex_Vectors.Vector(0..d-1) := Derivative(cp);
    err,rco,res : Standard_Floating_Vectors.Vector(z'range);
    nb : natural32;
    fail : boolean;

  begin
    z := Random_Vector(1,d);
    Silent_Durand_Kerner(cp,z,r,max,eps,nb,fail);
    if fail
     then put_line(file,"precision insufficient to reach results");
    end if;
    Write_Results(file,nb,z,r);
    Newton(cp,dcp,z,err,rco,res);
    sols := Create_Solution_List(z,err,rco,res);
  end QuadDobl_Find_Roots;

  procedure Multprec_Find_Roots
              ( file : in file_type; d : in integer32; size : in natural32;
                cp : in Multprec_Complex_Vectors.Vector;
                sols : out Multprec_Complex_Solutions.Solution_List ) is

    use Multprec_Durand_Kerner;

    max : constant natural32 := 10*size*natural32(d);
    deci : constant natural32
         := Multprec_Floating_Numbers.Size_to_Decimal(size);
    d_eps : constant double_float := 10.0**natural(-deci + 8);
    eps : Multprec_Floating_Numbers.Floating_Number
        := Multprec_Floating_Numbers.create(d_eps);
    z,r : Multprec_Complex_Vectors.Vector(1..d);
    dcp : Multprec_Complex_Vectors.Vector(0..d-1) := Derivative(cp);
    err,rco,res : Multprec_Floating_Vectors.Vector(z'range);
    nb : natural32;
    fail : boolean;

  begin
    z := Random_Vector1(1,d);
    Multprec_Complex_Vector_Tools.Set_Size(z,size);
    Multprec_Complex_Vector_Tools.Set_Size(dcp,size);
    Silent_Durand_Kerner(cp,z,r,max,eps,nb,fail);
    if fail
     then put_line(file,"precision insufficient to reach results");
    end if;
    Write_Results(file,nb,z,r);
    Newton(cp,dcp,z,err,rco,res);
    sols := Create_Solution_List(z,err,rco,res);
    Multprec_Floating_Numbers.Clear(eps);
    Multprec_Complex_Vectors.Clear(z);
    Multprec_Complex_Vectors.Clear(r);
  end Multprec_Find_Roots;

  function Coefficient_Vector
              ( d : natural32; p : Standard_Complex_Polynomials.Poly )
              return Standard_Complex_Vectors.Vector is

  -- DECRIPTION :
  --   Returns the coefficient vector of p, eventually divided 
  --   by its leading coefficient if p is not monic.

  -- REQUIRED : d is the degree of p.

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;

    res : Standard_Complex_Vectors.Vector(0..integer32(d))
        := (0..integer32(d) => Create(0.0));

    procedure Scan_Coefficient ( t : in Term; cont : out boolean ) is
    begin
      res(integer32(t.dg(1))) := t.cf;
      cont := true;
    end Scan_Coefficient;
    procedure Scan_Coefficients is new Visiting_Iterator(Scan_Coefficient);

  begin
    Scan_Coefficients(p);
   -- for i in 0..d-1 loop  -- why leave the leading coefficient !??
    for i in 0..integer32(d) loop
      res(i) := res(i)/res(integer32(d));
    end loop;
    return res;
  end Coefficient_Vector;

  function Coefficient_Vector
              ( d : natural32; p : DoblDobl_Complex_Polynomials.Poly )
              return DoblDobl_Complex_Vectors.Vector is

  -- DECRIPTION :
  --   Returns the coefficient vector of p, eventually divided 
  --   by its leading coefficient if p is not monic.

  -- REQUIRED : d is the degree of p.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Polynomials;

    res : DoblDobl_Complex_Vectors.Vector(0..integer32(d))
        := (0..integer32(d) => Create(integer(0)));

    procedure Scan_Coefficient ( t : in Term; cont : out boolean ) is
    begin
      res(integer32(t.dg(1))) := t.cf;
      cont := true;
    end Scan_Coefficient;
    procedure Scan_Coefficients is new Visiting_Iterator(Scan_Coefficient);

  begin
    Scan_Coefficients(p);
   -- for i in 0..d-1 loop  -- why leave the leading coefficient !??
    for i in 0..integer32(d) loop
      res(i) := res(i)/res(integer32(d));
    end loop;
    return res;
  end Coefficient_Vector;

  function Coefficient_Vector
              ( d : natural32; p : QuadDobl_Complex_Polynomials.Poly )
              return QuadDobl_Complex_Vectors.Vector is

  -- DECRIPTION :
  --   Returns the coefficient vector of p, eventually divided 
  --   by its leading coefficient if p is not monic.

  -- REQUIRED : d is the degree of p.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Polynomials;

    res : QuadDobl_Complex_Vectors.Vector(0..integer32(d))
        := (0..integer32(d) => Create(integer(0)));

    procedure Scan_Coefficient ( t : in Term; cont : out boolean ) is
    begin
      res(integer32(t.dg(1))) := t.cf;
      cont := true;
    end Scan_Coefficient;
    procedure Scan_Coefficients is new Visiting_Iterator(Scan_Coefficient);

  begin
    Scan_Coefficients(p);
   -- for i in 0..d-1 loop  -- why leave the leading coefficient !??
    for i in 0..integer32(d) loop
      res(i) := res(i)/res(integer32(d));
    end loop;
    return res;
  end Coefficient_Vector;

  function Coefficient_Vector
              ( d : natural32; p : Multprec_Complex_Polynomials.Poly )
              return Multprec_Complex_Vectors.Vector is

  -- DECRIPTION :
  --   Returns the coefficient vector of p, eventually divided 
  --   by its leading coefficient if p is not monic.

  -- REQUIRED : d is the degree of p.

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Polynomials;

    res : Multprec_Complex_Vectors.Vector(0..integer32(d))
        := (0..integer32(d) => Create(integer(0)));

    procedure Scan_Coefficient ( t : in Term; cont : out boolean ) is
    begin
      Copy(t.cf,res(integer32(t.dg(1))));
      cont := true;
    end Scan_Coefficient;
    procedure Scan_Coefficients is new Visiting_Iterator(Scan_Coefficient);

  begin
    Scan_Coefficients(p);
   -- for i in 0..d-1 loop  -- why leave the leading coefficient !??
    for i in 0..integer32(d) loop
      Div(res(i),res(integer32(d)));
    end loop;
    return res;
  end Coefficient_Vector;

  function Coefficient_Vector
              ( mind,maxd : integer32;
                p : Standard_Complex_Laurentials.Poly )
              return Standard_Complex_Vectors.Vector is

  -- DECRIPTION :
  --   Returns the coefficient vector of p, eventually divided 
  --   by its leading coefficient if p is not monic.

  -- REQUIRED : mind is the minimal and maxd the maximal degree of p.

    use Standard_Complex_Numbers;
    use Standard_Complex_Laurentials;

    res : Standard_Complex_Vectors.Vector(mind..maxd)
        := (mind..maxd => Create(0.0));

    procedure Scan_Coefficient ( t : in Term; cont : out boolean ) is
    begin
      res(t.dg(1)) := t.cf;
      cont := true;
    end Scan_Coefficient;
    procedure Scan_Coefficients is new Visiting_Iterator(Scan_Coefficient);

  begin
    Scan_Coefficients(p);
   -- for i in mind..maxd-1 loop  -- same mistake as above !?
    for i in mind..maxd-1 loop
      res(i) := res(i)/res(maxd);
    end loop;
    return res;
  end Coefficient_Vector;

  function Coefficient_Vector
              ( mind,maxd : integer32;
                p : Multprec_Complex_Laurentials.Poly )
              return Multprec_Complex_Vectors.Vector is

  -- DECRIPTION :
  --   Returns the coefficient vector of p, eventually divided 
  --   by its leading coefficient if p is not monic.

  -- REQUIRED : mind is the minimal and maxd the maximal degree of p.

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Laurentials;

    res : Multprec_Complex_Vectors.Vector(mind..maxd)
        := (mind..maxd => Create(integer(0)));

    procedure Scan_Coefficient ( t : in Term; cont : out boolean ) is
    begin
      Copy(t.cf,res(t.dg(1)));
      cont := true;
    end Scan_Coefficient;
    procedure Scan_Coefficients is new Visiting_Iterator(Scan_Coefficient);

  begin
    Scan_Coefficients(p);
   -- for i in mind..maxd-1 loop  -- same mistake as above !?
    for i in mind..maxd-1 loop
      Div(res(i),res(maxd));
    end loop;
    return res;
  end Coefficient_Vector;

  function Shift ( v : Standard_Complex_Vectors.Vector )
                 return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   If v'first < 0, then all entries will be shifted so the
  --   vector on return starts at 0.

    res : Standard_Complex_Vectors.Vector(0..v'last-v'first);
    ind : integer32 := 0;

  begin
    for i in v'range loop
      res(ind) := v(i);
      ind := ind + 1;
    end loop;
    return res;
  end Shift;

  function Shift ( v : Multprec_Complex_Vectors.Vector )
                 return Multprec_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   If v'first < 0, then all entries will be shifted so the
  --   vector on return starts at 0.

    res : Multprec_Complex_Vectors.Vector(0..v'last-v'first);
    ind : integer32 := 0;

  begin
    for i in v'range loop
      Multprec_Complex_Numbers.Copy(v(i),res(ind));
      ind := ind + 1;
    end loop;
    return res;
  end Shift;

  procedure Call_Durand_Kerner
              ( file : in file_type; d : in integer32;
                p : in Standard_Complex_Polynomials.Poly;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    cff : Standard_Complex_Vectors.Vector(0..d);

  begin
    tstart(timer);
    cff := Coefficient_Vector(natural32(d),p);
    new_line(file);
    put_line(file,"The coefficient vector :");
    put_line(file,cff);
    if d = 1 
     then sols := Create_Solution_List(-cff(0));
     else Standard_Find_Roots(file,d,cff,sols);
    end if;
    tstop(timer);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    new_line(file);
    print_times(file,timer,"application of the method of Weierstrass");
  end Call_Durand_Kerner;

  procedure Call_Durand_Kerner
              ( file : in file_type; d : in integer32;
                p : in DoblDobl_Complex_Polynomials.Poly;
                sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    timer : Timing_Widget;
    cff : DoblDobl_Complex_Vectors.Vector(0..d);

  begin
    tstart(timer);
    cff := Coefficient_Vector(natural32(d),p);
    new_line(file);
    put_line(file,"The coefficient vector :");
    put_line(file,cff);
    if d = 1 
     then sols := Create_Solution_List(-cff(0));
     else DoblDobl_Find_Roots(file,d,cff,sols);
    end if;
    tstop(timer);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    new_line(file);
    print_times(file,timer,"application of the method of Weierstrass");
  end Call_Durand_Kerner;

  procedure Call_Durand_Kerner
              ( file : in file_type; d : in integer32;
                p : in QuadDobl_Complex_Polynomials.Poly;
                sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    timer : Timing_Widget;
    cff : QuadDobl_Complex_Vectors.Vector(0..d);

  begin
    tstart(timer);
    cff := Coefficient_Vector(natural32(d),p);
    new_line(file);
    put_line(file,"The coefficient vector :");
    put_line(file,cff);
    if d = 1 
     then sols := Create_Solution_List(-cff(0));
     else QuadDobl_Find_Roots(file,d,cff,sols);
    end if;
    tstop(timer);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    new_line(file);
    print_times(file,timer,"application of the method of Weierstrass");
  end Call_Durand_Kerner;

  procedure Call_Durand_Kerner
              ( file : in file_type; mind,maxd : in integer32;
                p : in Standard_Complex_Laurentials.Poly;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    cff : Standard_Complex_Vectors.Vector(mind..maxd);

  begin
    tstart(timer);
    cff := Coefficient_Vector(mind,maxd,p);
    new_line(file);
    put_line(file,"The coefficient vector :");
    put_line(file,cff);
    if mind = 0
     then Standard_Find_Roots(file,maxd,cff,sols);
     else Standard_Find_Roots(file,maxd-mind,Shift(cff),sols);
    end if;
    tstop(timer);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    new_line(file);
    print_times(file,timer,"application of the method of Weierstrass");
  end Call_Durand_Kerner;

  procedure Call_Durand_Kerner
              ( file : in file_type; mind,maxd : in integer32;
                p : in Multprec_Complex_Laurentials.Poly;
                size : in natural32;
                sols : out Multprec_Complex_Solutions.Solution_List ) is

    use Multprec_Complex_Solutions;

    timer : Timing_Widget;
    cff : Multprec_Complex_Vectors.Vector(mind..maxd);

  begin
    tstart(timer);
    cff := Coefficient_Vector(mind,maxd,p);
    new_line(file);
    put_line(file,"The coefficient vector :");
    put_line(file,cff);
    if mind = 0
     then Multprec_Find_Roots(file,maxd,size,cff,sols);
     else Multprec_Find_Roots(file,maxd-mind,size,Shift(cff),sols);
    end if;
    tstop(timer);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    new_line(file);
    print_times(file,timer,"application of the method of Weierstrass");
  end Call_Durand_Kerner;

  procedure Black_Box_Durand_Kerner
              ( p : in Standard_Complex_Polynomials.Poly;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Numbers;

    d : constant integer32 := Standard_Complex_Polynomials.Degree(p);
    cff : constant Standard_Complex_Vectors.Vector(0..d)
        := Coefficient_Vector(natural32(d),p);

  begin
    if d = 1 then
      sols := Create_Solution_List(-cff(0));
    elsif d > 1 then
      Standard_Find_Roots(d,cff,sols);
    end if;
  end Black_Box_Durand_Kerner;

  procedure Black_Box_Durand_Kerner
              ( p : in DoblDobl_Complex_Polynomials.Poly;
                sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Numbers;

    d : constant integer32 := DoblDobl_Complex_Polynomials.Degree(p);
    cff : constant DoblDobl_Complex_Vectors.Vector(0..d)
        := Coefficient_Vector(natural32(d),p);

  begin
    if d = 1 then
      sols := Create_Solution_List(-cff(0));
    elsif d > 1 then
      DoblDobl_Find_Roots(d,cff,sols);
    end if;
  end Black_Box_Durand_Kerner;

  procedure Black_Box_Durand_Kerner
              ( p : in QuadDobl_Complex_Polynomials.Poly;
                sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Numbers;

    d : constant integer32 := QuadDobl_Complex_Polynomials.Degree(p);
    cff : constant QuadDobl_Complex_Vectors.Vector(0..d)
        := Coefficient_Vector(natural32(d),p);

  begin
    if d = 1 then
      sols := Create_Solution_List(-cff(0));
    elsif d > 1 then
      QuadDobl_Find_Roots(d,cff,sols);
    end if;
  end Black_Box_Durand_Kerner;

  procedure Black_Box_Durand_Kerner
              ( p : in Multprec_Complex_Polynomials.Poly;
                size : in natural32;
                sols : out Multprec_Complex_Solutions.Solution_List ) is

    use Multprec_Complex_Numbers;

    d : constant integer32 := Multprec_Complex_Polynomials.Degree(p);
    cff : Multprec_Complex_Vectors.Vector(0..d)
        := Coefficient_Vector(natural32(d),p);

  begin
    Multprec_Complex_Vector_Tools.Set_Size(cff,size);
    if d = 1 then
      sols := Create_Solution_List(-cff(0));
    elsif d > 1 then
      Multprec_Find_Roots(d,size,cff,sols);
    end if;
    Multprec_Complex_Vectors.Clear(cff);
  end Black_Box_Durand_Kerner;

  procedure Black_Box_Durand_Kerner
              ( p : in Standard_Complex_Laurentials.Poly;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Laurentials;
    mindeg : constant integer32 := Minimal_Degree(p,1);
    maxdeg : constant integer32 := Maximal_Degree(p,1);
    c : constant Standard_Complex_Vectors.Vector(mindeg..maxdeg)
      := Coefficient_Vector(mindeg,maxdeg,p);

  begin
    if mindeg = 0
     then Standard_Find_Roots(maxdeg,c,sols);
     else Standard_Find_Roots(maxdeg-mindeg,Shift(c),sols);
    end if;
  end Black_Box_Durand_Kerner;

  procedure Black_Box_Durand_Kerner
              ( p : in Multprec_Complex_Laurentials.Poly;
                size : in natural32;
                sols : out Multprec_Complex_Solutions.Solution_List ) is

    use Multprec_Complex_Laurentials;
    mindeg : constant integer32 := Minimal_Degree(p,1);
    maxdeg : constant integer32 := Maximal_Degree(p,1);
    c : constant Multprec_Complex_Vectors.Vector(mindeg..maxdeg)
      := Coefficient_Vector(mindeg,maxdeg,p);

  begin
    if mindeg = 0
     then Multprec_Find_Roots(maxdeg,size,c,sols);
     else Multprec_Find_Roots(maxdeg-mindeg,size,Shift(c),sols);
    end if;
  end Black_Box_Durand_Kerner;

  procedure Black_Box_Durand_Kerner
              ( file : in file_type;
                p : in Standard_Complex_Polynomials.Poly;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    d : constant integer32 := Standard_Complex_Polynomials.Degree(p);

  begin
    put_line(file,"1 1");
    put(file,p); new_line(file);
    if d = 0 then
      new_line(file);
      put_line(file,"There are no roots to compute.");
      new_line(file);
    else
      Call_Durand_Kerner(file,d,p,sols);
    end if;
  end Black_Box_Durand_Kerner;

  procedure Black_Box_Durand_Kerner
              ( file : in file_type;
                p : in DoblDobl_Complex_Polynomials.Poly;
                sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    d : constant integer32 := DoblDobl_Complex_Polynomials.Degree(p);

  begin
    put_line(file,"1 1");
    put(file,p); new_line(file);
    if d = 0 then
      new_line(file);
      put_line(file,"There are no roots to compute.");
      new_line(file);
    else
      Call_Durand_Kerner(file,d,p,sols);
    end if;
  end Black_Box_Durand_Kerner;

  procedure Black_Box_Durand_Kerner
              ( file : in file_type;
                p : in QuadDobl_Complex_Polynomials.Poly;
                sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    d : constant integer32 := QuadDobl_Complex_Polynomials.Degree(p);

  begin
    put_line(file,"1 1");
    put(file,p); new_line(file);
    if d = 0 then
      new_line(file);
      put_line(file,"There are no roots to compute.");
      new_line(file);
    else
      Call_Durand_Kerner(file,d,p,sols);
    end if;
  end Black_Box_Durand_Kerner;

  procedure Black_Box_Durand_Kerner
              ( file : in file_type;
                p : in Standard_Complex_Laurentials.Poly;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Laurentials;
    mind : constant integer32 := Minimal_Degree(p,1);
    maxd : constant integer32 := Maximal_Degree(p,1);

  begin
    put_line(file,"1 1");
    put(file,p); new_line(file);
    if mind = maxd then
      new_line(file);
      put_line(file,"There are no roots to compute.");
      new_line(file);
    else
      Call_Durand_Kerner(file,mind,maxd,p,sols);
    end if;
  end Black_Box_Durand_Kerner;

  procedure Black_Box_Durand_Kerner
              ( file : in file_type;
                p : in Multprec_Complex_Laurentials.Poly;
                size : in natural32;
                sols : out Multprec_Complex_Solutions.Solution_List ) is

    use Multprec_Complex_Laurentials;
    mind : constant integer32 := Minimal_Degree(p,1);
    maxd : constant integer32 := Maximal_Degree(p,1);

  begin
    put_line(file,"1 1");
    put(file,p); new_line(file);
    if mind = maxd then
      new_line(file);
      put_line(file,"There are no roots to compute.");
      new_line(file);
    else
      Call_Durand_Kerner(file,mind,maxd,p,size,sols);
    end if;
  end Black_Box_Durand_Kerner;

end Black_Box_Univariate_Solvers;
