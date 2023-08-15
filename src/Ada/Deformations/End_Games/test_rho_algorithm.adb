with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Mathematical_Functions;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Double_Rho_Algorithm;

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

  procedure Main is

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
  end Main;

end Test_Rho_Algorithm;
