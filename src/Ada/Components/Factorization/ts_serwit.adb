with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with Standard_Complex_Poly_Systems;
with Standard_Series_Matrix_Solvers;
with Standard_Speelpenning_Convolutions;
with System_Convolution_Circuits;        use System_Convolution_Circuits;
with Standard_Complex_Solutions;
with Witness_Sets_io;

procedure ts_serwit is

-- DESCRIPTION :
--   Develops series of a curve defined by a witness set.

  function Series_Coefficients
             ( v : Standard_Complex_Vectors.Vector;
               d : integer32 )
             return Standard_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns the coefficients of a vector of series,
  --   with leading coefficients the components of a solution v.
  --   The series are of degree d.

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

  procedure Add_Continuation_Parameter
              ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit ) is

  -- DESCRIPTION :
  --   To every coefficient adds 't' as the linear coefficient
  --   in all power series coefficients of the circuit.

    cf : Standard_Complex_Vectors.Link_to_Vector;

    use Standard_Complex_Vectors;

  begin
    for k in c.cff'range loop
      cf := c.cff(k);
      cf(1) := Standard_Complex_Numbers.Create(1.0);
    end loop;
    if c.cst /= null then
      c.cst(1) := Standard_Complex_Numbers.Create(1.0);
    end if;
  end Add_Continuation_Parameter;

  procedure Add_Continuation_Parameter
              ( c : in Standard_Speelpenning_Convolutions.Circuits ) is

  -- DESCRIPTION :
  --   To every coefficient adds 't' as the linear coefficient
  --   in all power series coefficients of the circuits.

  begin
    for k in c'range loop
      Add_Continuation_Parameter(c(k));
    end loop;
  end Add_Continuation_Parameter;

  procedure Minus ( v : in Standard_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Flips the sign of all coefficients in v.

    cf : Standard_Complex_Vectors.Link_to_Vector;

  begin
    for i in v'range loop
      cf := v(i);
      for j in cf'range loop
        Standard_Complex_Numbers.Min(cf(j));
      end loop;
    end loop;
  end Minus;

  procedure Update ( x,y : in Standard_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Adds to all coefficients of x the corresponding coefficient of y.

  -- REQUIRED : x'range = y'range, and for all k in x'range
  --   x(k)'range = y(k)'range.

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

  procedure Standard_Newton_Step
              ( s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                info : out integer32;
                ipvt : in out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Applies one Newton step on the convolution circuits c,
  --   departing from the series coefficients in s.

  -- ON ENTRY :
  --   s        system of convolution circuits;
  --   scf      vector of coefficients of power series;
  --   ipvt     vector for the pivoting information in the LU factorization;
  --   wrk      work space for the matrix series solver.

  -- ON RETURN :
  --   info     info from the LU factorization;
  --   ipvt     pivoting of the LU factorization on the lead matrix.

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

  procedure Standard_Newton_Steps
              ( s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                dim : in integer32 ) is

  -- DESCRIPTION :
  --   Applies several Newton steps on the system of convolution circuits s,
  --   departing from the series coefficients in scf.

    info : integer32;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    wrk : Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector(1..dim);
    ans : character;

  begin
    loop
      Standard_Newton_Step(s,scf,info,ipvt,wrk);
      put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    Standard_Complex_Vectors.Clear(wrk);
  end Standard_Newton_Steps;

  procedure Standard_Test ( degree : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a witness set.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    genpts : Standard_Complex_Solutions.Solution_List;
    deg,dim : natural32;

    use Standard_Speelpenning_Convolutions;

  begin
    Witness_Sets_io.Standard_Read_Embedding(lp,genpts,dim);
    deg := Standard_Complex_Solutions.Length_Of(genpts);
    new_line;
    put("Read witness set of dimension "); 
    put(dim,1); put(" and degree "); put(deg,1); put_line(".");
    declare
      c : constant Circuits(lp'range)
        := Make_Convolution_Circuits(lp.all,natural32(degree));
      s : Link_to_System := Create(c,lp'last,degree);
      sol : constant Standard_Complex_Solutions.Link_to_Solution
          := Standard_Complex_Solutions.Head_Of(genpts);
      scf : constant Standard_Complex_VecVecs.VecVec(1..sol.n)
          := Series_Coefficients(sol.v,degree);
    begin
      Add_Continuation_Parameter(c);
      put_line("The coefficients of the convolution circuits :");
      for i in c'range loop
        put_line(c(i).cff);
      end loop;
      Standard_Newton_Steps(s,scf,lp'last);
      Clear(s);
    end;
  end Standard_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the degree of the series.
  --   Launches the test.

    deg : integer32 := 0;

  begin
    new_line;
    put_line("Developing series of a curve defined by a witness set ...");
    new_line;
    put("Give the degree of the power series : "); get(deg); skip_line;
    Standard_Test(deg);
  end Main;

begin
  Main;
end ts_serwit;
