with text_io;                            use text_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with DoblDobl_Complex_VecVecs_io;        use DoblDobl_Complex_VecVecs_io;
with QuadDobl_Complex_VecVecs_io;        use QUadDobl_Complex_VecVecs_io;
with Standard_Series_Matrix_Solvers;
with DoblDobl_Series_Matrix_Solvers;
with QuadDobl_Series_Matrix_Solvers;

package body Newton_Convolutions is

  function Series_Coefficients
             ( v : Standard_Complex_Vectors.Vector;
               d : integer32 )
             return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(v'range);

  begin
    for k in res'range loop
      declare
        cff : Standard_Complex_Vectors.Vector(0..d);
      begin
        cff(0) := v(k);
        cff(1..d) := (1..d => Standard_Complex_Numbers.Create(0.0));
        res(k) := new Standard_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end Series_Coefficients;

  function Series_Coefficients
             ( v : DoblDobl_Complex_Vectors.Vector;
               d : integer32 )
             return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(v'range);
    zero : constant double_double := create(0.0);

  begin
    for k in res'range loop
      declare
        cff : DoblDobl_Complex_Vectors.Vector(0..d);
      begin
        cff(0) := v(k);
        cff(1..d) := (1..d => DoblDobl_Complex_Numbers.Create(zero));
        res(k) := new DoblDobl_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end Series_Coefficients;

  function Series_Coefficients
             ( v : QuadDobl_Complex_Vectors.Vector;
               d : integer32 )
             return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(v'range);
    zero : constant quad_double := create(0.0);

  begin
    for k in res'range loop
      declare
        cff : QuadDobl_Complex_Vectors.Vector(0..d);
      begin
        cff(0) := v(k);
        cff(1..d) := (1..d => QuadDobl_Complex_Numbers.Create(zero));
        res(k) := new QuadDobl_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end Series_Coefficients;

  procedure Minus ( v : in Standard_Complex_VecVecs.VecVec ) is

    cf : Standard_Complex_Vectors.Link_to_Vector;

  begin
    for i in v'range loop
      cf := v(i);
      for j in cf'range loop
        Standard_Complex_Numbers.Min(cf(j));
      end loop;
    end loop;
  end Minus;

  procedure Minus ( v : in DoblDobl_Complex_VecVecs.VecVec ) is

    cf : DoblDobl_Complex_Vectors.Link_to_Vector;

  begin
    for i in v'range loop
      cf := v(i);
      for j in cf'range loop
        DoblDobl_Complex_Numbers.Min(cf(j));
      end loop;
    end loop;
  end Minus;

  procedure Minus ( v : in QuadDobl_Complex_VecVecs.VecVec ) is

    cf : QuadDobl_Complex_Vectors.Link_to_Vector;

  begin
    for i in v'range loop
      cf := v(i);
      for j in cf'range loop
        QuadDobl_Complex_Numbers.Min(cf(j));
      end loop;
    end loop;
  end Minus;

  procedure Update ( x,y : in Standard_Complex_VecVecs.VecVec ) is

    xcf,ycf : Standard_Complex_Vectors.Link_to_Vector;

  begin
    for i in x'range loop
      xcf := x(i);
      ycf := y(i);
      for j in xcf'range loop
        Standard_Complex_Numbers.Add(xcf(j),ycf(j));
      end loop;
    end loop;
  end Update;

  procedure Update ( x,y : in DoblDobl_Complex_VecVecs.VecVec ) is

    xcf,ycf : DoblDobl_Complex_Vectors.Link_to_Vector;

  begin
    for i in x'range loop
      xcf := x(i);
      ycf := y(i);
      for j in xcf'range loop
        DoblDobl_Complex_Numbers.Add(xcf(j),ycf(j));
      end loop;
    end loop;
  end Update;

  procedure Update ( x,y : in QuadDobl_Complex_VecVecs.VecVec ) is

    xcf,ycf : QuadDobl_Complex_Vectors.Link_to_Vector;

  begin
    for i in x'range loop
      xcf := x(i);
      ycf := y(i);
      for j in xcf'range loop
        QuadDobl_Complex_Numbers.Add(xcf(j),ycf(j));
      end loop;
    end loop;
  end Update;

  procedure Standard_Newton_Step
              ( s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                info : out integer32;
                ipvt : in out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector ) is
  begin
    put_line("scf :"); put_line(scf);
    Standard_Speelpenning_Convolutions.Compute(s.pwt,s.mxe,scf);
    Standard_Speelpenning_Convolutions.EvalDiff(s,scf);
    put_line("vy :"); put_line(s.vy);
    Minus(s.vy);
    Standard_Series_Matrix_Solvers.Solve_by_lufac(s.vm,s.vy,ipvt,info,wrk);
    put_line("dx :"); put_line(s.vy);
    Standard_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    Update(scf,s.yv);
  end Standard_Newton_Step;

  procedure DoblDobl_Newton_Step
              ( s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                info : out integer32;
                ipvt : in out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector ) is
  begin
    put_line("scf :"); put_line(scf);
    DoblDobl_Speelpenning_Convolutions.Compute(s.pwt,s.mxe,scf);
    DoblDobl_Speelpenning_Convolutions.EvalDiff(s,scf);
    put_line("vy :"); put_line(s.vy);
    Minus(s.vy);
    DoblDobl_Series_Matrix_Solvers.Solve_by_lufac(s.vm,s.vy,ipvt,info,wrk);
    put_line("dx :"); put_line(s.vy);
    DoblDobl_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    Update(scf,s.yv);
  end DoblDobl_Newton_Step;

  procedure QuadDobl_Newton_Step
              ( s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in QuadDobl_Complex_VecVecs.VecVec;
                info : out integer32;
                ipvt : in out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector ) is
  begin
    put_line("scf :"); put_line(scf);
    QuadDobl_Speelpenning_Convolutions.Compute(s.pwt,s.mxe,scf);
    QuadDobl_Speelpenning_Convolutions.EvalDiff(s,scf);
    put_line("vy :"); put_line(s.vy);
    Minus(s.vy);
    QuadDobl_Series_Matrix_Solvers.Solve_by_lufac(s.vm,s.vy,ipvt,info,wrk);
    put_line("dx :"); put_line(s.vy);
    QuadDobl_Speelpenning_Convolutions.Delinearize(s.vy,s.yv);
    Update(scf,s.yv);
  end QuadDobl_Newton_Step;

end Newton_Convolutions;
