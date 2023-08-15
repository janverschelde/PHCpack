with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Mathematical_Functions;
with Double_Double_Constants;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Double_Double_Vectors;
with Double_Double_VecVecs;
with Double_Rho_Algorithm;
with Double_Double_Rho_Algorithm;

package body Test_Rho_Algorithm is

  procedure Double_Example3
              ( nbr : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   The vector nbr on return contains the sum of 1/k^2,
  --   for k from 1 to m, which converges to pi^2/6.
  --   This is example 3 of the Brezinski paper and
  --   the example in Table 3 of the 1955 Wynn paper on the rho algorithm.
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
  end Double_Example3;

  procedure Double_Double_Example3
              ( nbr : out Double_Double_Vectors.Vector ) is

  -- DESCRIPTION :
  --   The vector nbr on return contains the sum of 1/k^2,
  --   for k from 1 to m, which converges to pi^2/6.
  --   This is example 3 of the Brezinski paper and
  --   the example in Table 3 of the 1955 Wynn paper on the rho algorithm.
  --   The epsilon algorithm does not work on this example.

  -- REQUIRED : nbr'first = 0, which corresponds to m = 1.

    one : constant double_double := Double_Double_Numbers.Create(1.0);
    cff,den : double_double;
    kp1 : integer32;

  begin
    nbr(0) := one;
    for k in 1..nbr'last loop
      kp1 := k+1;
      den := Double_Double_Numbers.Create(kp1*kp1);
      cff := one/den;
      nbr(k) := nbr(k-1) + cff;
    end loop;
  end Double_Double_Example3;

  procedure Double_Example4
              ( nbr : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   The vector nbr on return contains the sum of 1/(2*k*(2*k-1)),
  --   for k from 1 to m, which converges to ln(2).
  --   This is example in Table 4 of the 1955 Wynn paper
  --   on the rho algorithm.

  -- REQUIRED : nbr'first = 0, which corresponds to m = 1.

    cff : double_float;
    den : integer32;

  begin
    nbr(0) := 0.5;
    for k in 1..nbr'last loop
      den := 2*(k+1)*(2*k+1); -- 2*k*(2*k-1);
      cff := 1.0/double_float(den);
      nbr(k) := nbr(k-1) + cff;
    end loop;
  end Double_Example4;

  procedure Double_Double_Example4
              ( nbr : out Double_Double_Vectors.Vector ) is

  -- DESCRIPTION :
  --   The vector nbr on return contains the sum of 1/(2*k*(2*k-1)),
  --   for k from 1 to m, which converges to ln(2).
  --   This is example in Table 4 of the 1955 Wynn paper
  --   on the rho algorithm.

  -- REQUIRED : nbr'first = 0, which corresponds to m = 1.

    one : constant double_double := Double_Double_Numbers.create(1.0);
    cff,dd_den : double_double;
    den : integer32;

  begin
    nbr(0) := Double_Double_Numbers.create(0.5);
    for k in 1..nbr'last loop
      den := 2*(k+1)*(2*k+1); -- 2*k*(2*k-1);
      dd_den := Double_Double_Numbers.create(den);
      cff := one/dd_den;
      nbr(k) := nbr(k-1) + cff;
    end loop;
  end Double_Double_Example4;

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

  procedure Double_Double_Run_Rho
              ( tab : in Double_Double_VecVecs.VecVec;
                idx : in out Standard_Integer_Vectors.Vector;
                seq : in Double_Double_Vectors.Vector;
                exa : in double_double;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Runs the rho algorithm in double precision
  --   on the sequence seq, defining the extrapolation table tab,
  --   and computing the error using exa.

  -- REQUIRED : tab'last = idx'last = seq'last.

    zero : constant double_double := Double_Double_Numbers.create(0.0);
    dim : integer32 := 0;
    err : double_double;
 
  begin
    for k in seq'range loop
      if verbose then
        put("adding "); put(seq(k)); put_line(" ...");
      end if;
      Double_Double_Rho_Algorithm.Extrapolate(tab,dim,idx,seq(k),verbose);
    end loop;
    put("computed "); put(dim,1); put_line(" columns :");
    for col in 0..dim-1 loop
      if col mod 2 = 0 then
        put(col,6); put(" : "); put(tab(col)(idx(col)));
        if exa /= zero then
          err := abs(tab(col)(idx(col)) - exa);
          put("  error : "); put(err,3);
        end if;
        new_line;
      end if;
    end loop;
  end Double_Double_Run_Rho;

  procedure Double_Main is

  -- DESCRIPTION :
  --   Runs tests in double precision.

    dim : constant integer32 := 11;

  begin
   -- put("Give the dimension : "); get(dim);
    declare
      table : Standard_Floating_VecVecs.VecVec(0..dim);
      tblix : Standard_Integer_Vectors.Vector(0..dim);
      seq3 : Standard_Floating_Vectors.Vector(0..dim);
      exa3 : constant double_float
           := Standard_Mathematical_Functions.Pi**2/6.0;
      seq4 : Standard_Floating_Vectors.Vector(0..dim);
      exa4 : constant double_float
           := Standard_Mathematical_Functions.LN(2.0);
    begin
      put_line("running the Table 3 example ...");
      Double_Example3(seq3);
      Double_Rho_Algorithm.Allocate(table,tblix,dim);
      Double_Run_Rho(table,tblix,seq3,exa3,true);
      put("pi^2/6 : "); put(exa3); new_line;
      put_line("running the Table 4 example ...");
      Double_Example4(seq4);
      Double_Rho_Algorithm.Allocate(table,tblix,dim);
      Double_Run_Rho(table,tblix,seq4,exa4,true);
      put(" ln(2) : "); put(exa4); new_line;
    end;
  end Double_Main;

  procedure Double_Double_Main is

  -- DESCRIPTION :
  --   Runs tests in double double precision.

    dim : constant integer32 := 29;

  begin
   -- put("Give the dimension : "); get(dim);
    declare
      table : Double_Double_VecVecs.VecVec(0..dim);
      tblix : Standard_Integer_Vectors.Vector(0..dim);
      seq3 : Double_Double_Vectors.Vector(0..dim);
      six : constant double_double := Double_Double_Numbers.create(6.0);
      ddpi : constant double_double := Double_Double_Constants.pi;
      exa3 : constant double_double := (ddpi*ddpi)/six;
      seq4 : Double_Double_Vectors.Vector(0..dim);
      two : constant double_double := Double_Double_Numbers.create(2.0);
      exa4 : constant double_double := Double_Double_Numbers.log(two);
    begin
      put_line("running the Table 3 example ...");
      Double_Double_Example3(seq3);
      Double_Double_Rho_Algorithm.Allocate(table,tblix,dim);
      Double_Double_Run_Rho(table,tblix,seq3,exa3,false);
      put("pi^2/6 : "); put(exa3); new_line;
      put_line("running the Table 4 example ...");
      Double_Double_Example4(seq4);
      Double_Double_Rho_Algorithm.Allocate(table,tblix,dim);
      Double_Double_Run_Rho(table,tblix,seq4,exa4,false);
      put(" ln(2) : "); put(exa4); new_line;
    end;
  end Double_Double_Main;

  procedure Main is
  begin
    Double_Main;
    Double_Double_Main;
  end Main;

end Test_Rho_Algorithm;
