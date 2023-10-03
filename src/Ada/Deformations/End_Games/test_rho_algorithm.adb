with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Mathematical_Functions;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Double_Rho_Algorithm;
with Double_Complex_Rho_Algorithm;
with DoblDobl_Complex_Rho_Algorithm;

package body Test_Rho_Algorithm is

  procedure Double_Example
              ( nbr : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   The vector nbr on return contains the sum of 1/k^2,
  --   for k from 1 to m, which converges to pi^2/6.
  --   This is example 3 of the Brezinski paper and
  --   the example of the 1955 Wynn paper on the rho algorithm.
  --   The epsilon algorithm does not work on this example.

  -- REQUIRED : nbr'first = 0, which corresponds to m = 1.

    cff : double_float;
    kp1 : integer32;

  begin
    nbr(0) := 1.0;
    for k in 1..nbr'last loop
      kp1 := k+1;
      cff := 1.0/double_float(kp1*kp1);
      nbr(k) := nbr(k-1) + cff;
    end loop;
  end Double_Example;

  procedure Double_Complex_Example
              ( nbr : out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Defines the coefficients of the Taylor series for sqrt(I - t).
  --   The rho algorithm is correct for sqrt(1 - t)
  --   with only four terms. 
  --   The coefficients were computedwith ts_fabryhom.

  -- REQUIRED : nbr'first = 0 and nbr'last >= 4.

    c0 : constant Complex_Number
       := create(-7.07106781186548E-01, -7.07106781186548E-01);
    c1 : constant Complex_Number
       := create( 3.53553390593274E-01, -3.53553390593274E-01);
    c2 : constant Complex_Number
       := create(-8.83883476483185E-02, -8.83883476483185E-02);
    c3 : constant Complex_Number
       := create(-4.41941738241592E-02,  4.41941738241592E-02);
    c4 : constant Complex_Number
       := create( 2.76213586400995E-02,  2.76213586400995E-02);
    c5 : constant Complex_Number
       := create( 1.93349510480697E-02, -1.93349510480697E-02);
    c6 : constant Complex_Number
       := create(-1.45012132860523E-02, -1.45012132860523E-02);
    c7 : constant Complex_Number
       := create(-1.13938104390411E-02,  1.13938104390411E-02);
    c8 : constant Complex_Number
       := create( 9.25747098172088E-03,  9.25747098172088E-03);

  begin
    nbr(0) := c0/c1;
    nbr(1) := c1/c2;
    nbr(2) := c2/c3;
    nbr(3) := c3/c4;
    nbr(4) := c4/c5;
    if nbr'last >= 5 then
      nbr(5) := c5/c6;
    end if;
    if nbr'last >= 6 then
      nbr(6) := c6/c6;
    end if;
    if nbr'last >= 7 then
      nbr(7) := c7/c8;
    end if;
  end Double_Complex_Example;

  procedure Double_Run_Rho
              ( tab : in Standard_Floating_VecVecs.VecVec;
                idx : in out Standard_Integer_Vectors.Vector;
                seq : in Standard_Floating_Vectors.Vector;
                exa : in double_float := 0.0;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Runs the rho algorithm in double precision
  --   on the sequence seq, defining the extrapolation table tab,
  --   and computing the error using exa.

  -- REQUIRED : tab'last = idx'last = seq'last.

    dim : integer32 := 0;
    err : double_float;
 
  begin
    for k in seq'range loop
      if verbose then
        put("adding "); put(seq(k)); put_line(" ...");
      end if;
      Double_Rho_Algorithm.Extrapolate(tab,dim,idx,seq(k),verbose);
    end loop;
    put("computed "); put(dim,1); put_line(" columns :");
    for col in 0..dim-1 loop
      if col mod 2 = 0 then
        put(col,6); put(" : "); put(tab(col)(idx(col)));
        if exa /= 0.0 then
          err := abs(tab(col)(idx(col)) - exa);
          put("  error : "); put(err,3);
        end if;
        new_line;
      end if;
    end loop;
  end Double_Run_Rho;

  procedure Double_Complex_Run_Rho
              ( tab : in Standard_Complex_VecVecs.VecVec;
                idx : in out Standard_Integer_Vectors.Vector;
                seq : in Standard_Complex_Vectors.Vector;
                exa : in Complex_Number;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Runs the rho algorithm in double precision
  --   on the sequence seq, defining the extrapolation table tab,
  --   and computing the error using exa.

  -- REQUIRED : tab'last = idx'last = seq'last.

    dim : integer32 := 0;
    err : double_float;
    zero : constant Complex_Number := create(0.0);
 
  begin
    for k in seq'range loop
      if verbose then
        put("adding "); put(seq(k)); put_line(" ...");
      end if;
      Double_Complex_Rho_Algorithm.Extrapolate(tab,dim,idx,seq(k),verbose);
    end loop;
    put("computed "); put(dim,1); put_line(" columns :");
    for col in 0..dim-1 loop
      if col mod 2 = 0 then
        put(col,6); put(" : "); put(tab(col)(idx(col)));
        if exa /= zero then
          err := AbsVal(tab(col)(idx(col)) - exa);
          put("  error : "); put(err,3);
        end if;
        new_line;
      end if;
    end loop;
  end Double_Complex_Run_Rho;

  procedure Double_Main is

    dim : constant integer32 := 11;

  begin
   -- put("Give the dimension : "); get(dim);
    declare
      table : Standard_Floating_VecVecs.VecVec(0..dim);
      tblix : Standard_Integer_Vectors.Vector(0..dim);
      seq : Standard_Floating_Vectors.Vector(0..dim);
      exa : constant double_float
          := Standard_Mathematical_Functions.Pi**2/6.0;
    begin
      put_line("running an example ...");
      Double_Example(seq);
      Double_Rho_Algorithm.Allocate(table,tblix,dim);
      Double_Run_Rho(table,tblix,seq,exa,true);
      put("pi^2/6 : "); put(exa); new_line;
    end;
  end Double_Main;

  procedure Double_Complex_Main is

    dim : constant integer32 := 4;

  begin
   -- put("Give the dimension : "); get(dim);
    declare
      table : Standard_Complex_VecVecs.VecVec(0..dim);
      tblix : Standard_Integer_Vectors.Vector(0..dim);
      seq : Standard_Complex_Vectors.Vector(0..dim);
      exa : constant Complex_Number := create(0.0,1.0); 
    begin
      put_line("running a complex example ...");
      Double_Complex_Example(seq);
      Double_Complex_Rho_Algorithm.Allocate(table,tblix,dim);
      Double_Complex_Run_Rho(table,tblix,seq,exa,true);
    end;
  end Double_Complex_Main;

  procedure Main is
  begin
    Double_Main;
    Double_Complex_Main;
  end Main;

end Test_Rho_Algorithm;
