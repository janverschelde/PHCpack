with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Mathematical_Functions;
with Double_Double_Constants;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Constants;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Octo_Double_Constants;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Octo_Double_Numbers_io;             use Octo_Double_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Double_Double_Vectors;
with Double_Double_VecVecs;
with Quad_Double_Vectors;
with Quad_Double_VecVecs;
with Octo_Double_Vectors;
with Octo_Double_VecVecs;
with Double_Theta_Algorithm;
with Double_Double_Theta_Algorithm;
with Quad_Double_Theta_Algorithm;
with Octo_Double_Theta_Algorithm;

package body Test_Theta_Algorithm is

  procedure Double_Example2
              ( nbr : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   The vector nbr on return contains the sum of (-1)**(k+1)/k,
  --   for k from 1 to m, which converges to log(2).
  --   This is example 2 of the Brezinski paper.
  --   The epsilon algorithm works well on this example.

  -- REQUIRED : nbr'first = 0, which corresponds to m = 1.

    cff : double_float;

  begin
    nbr(0) := 1.0;
    for k in 1..nbr'last loop
      cff := 1.0/double_float(k+1);
      if k mod 2 = 1
       then cff := -cff;
      end if;
      nbr(k) := nbr(k-1) + cff;
    end loop;
  end Double_Example2;

  procedure Double_Double_Example2
              ( nbr : out Double_Double_Vectors.Vector ) is

  -- DESCRIPTION :
  --   The vector nbr on return contains the sum of (-1)**(k+1)/k,
  --   for k from 1 to m, which converges to log(2).
  --   This is example 2 of the Brezinski paper.
  --   The epsilon algorithm works well on this example.

  -- REQUIRED : nbr'first = 0, which corresponds to m = 1.

    one : constant double_double := Double_Double_Numbers.create(1.0);
    cff,kp1 : double_double;

  begin
    nbr(0) := Double_Double_Numbers.create(1.0);
    for k in 1..nbr'last loop
      kp1 := Double_Double_Numbers.create(k+1);
      cff := one/kp1;
      if k mod 2 = 1
       then cff := -cff;
      end if;
      nbr(k) := nbr(k-1) + cff;
    end loop;
  end Double_Double_Example2;

  procedure Quad_Double_Example2
              ( nbr : out Quad_Double_Vectors.Vector ) is

  -- DESCRIPTION :
  --   The vector nbr on return contains the sum of (-1)**(k+1)/k,
  --   for k from 1 to m, which converges to log(2).
  --   This is example 2 of the Brezinski paper.
  --   The epsilon algorithm works well on this example.

  -- REQUIRED : nbr'first = 0, which corresponds to m = 1.

    one : constant quad_double := Quad_Double_Numbers.create(1.0);
    cff,kp1 : quad_double;

  begin
    nbr(0) := Quad_Double_Numbers.create(1.0);
    for k in 1..nbr'last loop
      kp1 := Quad_Double_Numbers.create(k+1);
      cff := one/kp1;
      if k mod 2 = 1
       then cff := -cff;
      end if;
      nbr(k) := nbr(k-1) + cff;
    end loop;
  end Quad_Double_Example2;

  procedure Octo_Double_Example2
              ( nbr : out Octo_Double_Vectors.Vector ) is

  -- DESCRIPTION :
  --   The vector nbr on return contains the sum of (-1)**(k+1)/k,
  --   for k from 1 to m, which converges to log(2).
  --   This is example 2 of the Brezinski paper.
  --   The epsilon algorithm works well on this example.

  -- REQUIRED : nbr'first = 0, which corresponds to m = 1.

    one : constant octo_double := Octo_Double_Numbers.create(1.0);
    cff,kp1 : octo_double;

  begin
    nbr(0) := Octo_Double_Numbers.create(1.0);
    for k in 1..nbr'last loop
      kp1 := Octo_Double_Numbers.create(k+1);
      cff := one/kp1;
      if k mod 2 = 1
       then cff := -cff;
      end if;
      nbr(k) := nbr(k-1) + cff;
    end loop;
  end Octo_Double_Example2;

  procedure Double_Example3
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
  end Double_Example3;

  procedure Double_Double_Example3
              ( nbr : out Double_Double_Vectors.Vector ) is

  -- DESCRIPTION :
  --   The vector nbr on return contains the sum of 1/k^2,
  --   for k from 1 to m, which converges to pi^2/6.
  --   This is example 3 of the Brezinski paper and
  --   the example of the 1955 Wynn paper on the rho algorithm.
  --   The epsilon algorithm does not work on this example.

  -- REQUIRED : nbr'first = 0, which corresponds to m = 1.

    one : constant double_double := Double_Double_Numbers.create(1.0);
    cff,den : double_double;
    kp1 : integer32;

  begin
    nbr(0) := Double_Double_Numbers.create(1.0);
    for k in 1..nbr'last loop
      kp1 := k+1;
      den := Double_Double_Numbers.create(kp1*kp1);
      cff := one/den;
      nbr(k) := nbr(k-1) + cff;
    end loop;
  end Double_Double_Example3;

  procedure Quad_Double_Example3
              ( nbr : out Quad_Double_Vectors.Vector ) is

  -- DESCRIPTION :
  --   The vector nbr on return contains the sum of 1/k^2,
  --   for k from 1 to m, which converges to pi^2/6.
  --   This is example 3 of the Brezinski paper and
  --   the example of the 1955 Wynn paper on the rho algorithm.
  --   The epsilon algorithm does not work on this example.

  -- REQUIRED : nbr'first = 0, which corresponds to m = 1.

    one : constant quad_double := Quad_Double_Numbers.create(1.0);
    cff,den : quad_double;
    kp1 : integer32;

  begin
    nbr(0) := Quad_Double_Numbers.create(1.0);
    for k in 1..nbr'last loop
      kp1 := k+1;
      den := Quad_Double_Numbers.create(kp1*kp1);
      cff := one/den;
      nbr(k) := nbr(k-1) + cff;
    end loop;
  end Quad_Double_Example3;

  procedure Octo_Double_Example3
              ( nbr : out Octo_Double_Vectors.Vector ) is

  -- DESCRIPTION :
  --   The vector nbr on return contains the sum of 1/k^2,
  --   for k from 1 to m, which converges to pi^2/6.
  --   This is example 3 of the Brezinski paper and
  --   the example of the 1955 Wynn paper on the rho algorithm.
  --   The epsilon algorithm does not work on this example.

  -- REQUIRED : nbr'first = 0, which corresponds to m = 1.

    one : constant octo_double := Octo_Double_Numbers.create(1.0);
    cff,den : octo_double;
    kp1 : integer32;

  begin
    nbr(0) := Octo_Double_Numbers.create(1.0);
    for k in 1..nbr'last loop
      kp1 := k+1;
      den := Octo_Double_Numbers.create(kp1*kp1);
      cff := one/den;
      nbr(k) := nbr(k-1) + cff;
    end loop;
  end Octo_Double_Example3;

  procedure Double_Run_Theta
              ( tab : in Standard_Floating_VecVecs.VecVec;
                idx : in out Standard_Integer_Vectors.Vector;
                seq : in Standard_Floating_Vectors.Vector;
                exa : in double_float := 0.0;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Runs the theta algorithm in double precision
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
      Double_Theta_Algorithm.Extrapolate(tab,dim,idx,seq(k),verbose);
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
  end Double_Run_Theta;

  procedure Double_Double_Run_Theta
              ( tab : in Double_Double_VecVecs.VecVec;
                idx : in out Standard_Integer_Vectors.Vector;
                seq : in Double_Double_Vectors.Vector;
                exa : in double_double;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Runs the theta algorithm in double double precision
  --   on the sequence seq, defining the extrapolation table tab,
  --   and computing the error using exa.

  -- REQUIRED : tab'last = idx'last = seq'last.

    dim : integer32 := 0;
    err : double_double;
    zero : constant double_double := Double_Double_Numbers.create(0.0);
 
  begin
    for k in seq'range loop
      if verbose then
        put("adding "); put(seq(k)); put_line(" ...");
      end if;
      Double_Double_Theta_Algorithm.Extrapolate(tab,dim,idx,seq(k),verbose);
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
  end Double_Double_Run_Theta;

  procedure Quad_Double_Run_Theta
              ( tab : in Quad_Double_VecVecs.VecVec;
                idx : in out Standard_Integer_Vectors.Vector;
                seq : in Quad_Double_Vectors.Vector;
                exa : in quad_double;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Runs the theta algorithm in quad double precision
  --   on the sequence seq, defining the extrapolation table tab,
  --   and computing the error using exa.

  -- REQUIRED : tab'last = idx'last = seq'last.

    dim : integer32 := 0;
    err : quad_double;
    zero : constant quad_double := Quad_Double_Numbers.create(0.0);
 
  begin
    for k in seq'range loop
      if verbose then
        put("adding "); put(seq(k)); put_line(" ...");
      end if;
      Quad_Double_Theta_Algorithm.Extrapolate(tab,dim,idx,seq(k),verbose);
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
  end Quad_Double_Run_Theta;

  procedure Octo_Double_Run_Theta
              ( tab : in Octo_Double_VecVecs.VecVec;
                idx : in out Standard_Integer_Vectors.Vector;
                seq : in Octo_Double_Vectors.Vector;
                exa : in octo_double;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Runs the theta algorithm in octo double precision
  --   on the sequence seq, defining the extrapolation table tab,
  --   and computing the error using exa.

  -- REQUIRED : tab'last = idx'last = seq'last.

    dim : integer32 := 0;
    err : octo_double;
    zero : constant octo_double := Octo_Double_Numbers.create(0.0);
 
  begin
    for k in seq'range loop
      if verbose then
        put("adding "); put(seq(k)); put_line(" ...");
      end if;
      Octo_Double_Theta_Algorithm.Extrapolate(tab,dim,idx,seq(k),verbose);
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
  end Octo_Double_Run_Theta;

  procedure Double_Main is

    dim : constant integer32 := 17;

  begin
   -- put("Give the dimension : "); get(dim);
    declare
      table : Standard_Floating_VecVecs.VecVec(0..dim);
      tblix : Standard_Integer_Vectors.Vector(0..dim);
      seq2 : Standard_Floating_Vectors.Vector(0..dim);
      exa2 : constant double_float
           := Standard_Mathematical_Functions.LN(2.0);
      seq3 : Standard_Floating_Vectors.Vector(0..dim);
      exa3 : constant double_float
           := Standard_Mathematical_Functions.Pi**2/6.0;
    begin
      put_line("running example 2 ...");
      Double_Example2(seq2);
      Double_Theta_Algorithm.Allocate(table,tblix,dim);
      Double_Run_Theta(table,tblix,seq2,exa2,false);
      put("log(2) : "); put(exa2); new_line;
      put_line("running example 3 ...");
      Double_Example3(seq3);
      Double_Theta_Algorithm.Allocate(table,tblix,dim);
      Double_Run_Theta(table,tblix,seq3,exa3,false);
      put("pi^2/6 : "); put(exa3); new_line;
    end;
  end Double_Main;

  procedure Double_Double_Main is

    dim : constant integer32 := 21;

  begin
   -- put("Give the dimension : "); get(dim);
    declare
      table : Double_Double_VecVecs.VecVec(0..dim);
      tblix : Standard_Integer_Vectors.Vector(0..dim);
      seq2 : Double_Double_Vectors.Vector(0..dim);
      two : constant double_double := Double_Double_Numbers.create(2.0);
      exa2 : constant double_double := Double_Double_Numbers.log(two);
      seq3 : Double_Double_Vectors.Vector(0..dim);
      six : constant double_double := Double_Double_Numbers.create(6.0);
      ddpi : constant double_double := Double_Double_Constants.pi;
      exa3 : constant double_double := (ddpi*ddpi)/six;
    begin
      put_line("running example 2 ...");
      Double_Double_Example2(seq2);
      Double_Double_Theta_Algorithm.Allocate(table,tblix,dim);
      Double_Double_Run_Theta(table,tblix,seq2,exa2,false);
      put("log(2) : "); put(exa2); new_line;
      put_line("running example 3 ...");
      Double_Double_Example3(seq3);
      Double_Double_Theta_Algorithm.Allocate(table,tblix,dim);
      Double_Double_Run_Theta(table,tblix,seq3,exa3,false);
      put("pi^2/6 : "); put(exa3); new_line;
    end;
  end Double_Double_Main;

  procedure Quad_Double_Main is

    dim : constant integer32 := 25;

  begin
   -- put("Give the dimension : "); get(dim);
    declare
      table : Quad_Double_VecVecs.VecVec(0..dim);
      tblix : Standard_Integer_Vectors.Vector(0..dim);
      seq2 : Quad_Double_Vectors.Vector(0..dim);
      two : constant quad_double := Quad_Double_Numbers.create(2.0);
      exa2 : constant quad_double := Quad_Double_Numbers.log(two);
      seq3 : Quad_Double_Vectors.Vector(0..dim);
      six : constant quad_double := Quad_Double_Numbers.create(6.0);
      qdpi : constant quad_double := Quad_Double_Constants.pi;
      exa3 : constant quad_double := (qdpi*qdpi)/six;
    begin
      put_line("running example 2 ...");
      Quad_Double_Example2(seq2);
      Quad_Double_Theta_Algorithm.Allocate(table,tblix,dim);
      Quad_Double_Run_Theta(table,tblix,seq2,exa2,false);
      put("log(2) : "); put(exa2); new_line;
      put_line("running example 3 ...");
      Quad_Double_Example3(seq3);
      Quad_Double_Theta_Algorithm.Allocate(table,tblix,dim);
      Quad_Double_Run_Theta(table,tblix,seq3,exa3,false);
      put("pi^2/6 : "); put(exa3); new_line;
    end;
  end Quad_Double_Main;

  procedure Octo_Double_Main is

    dim : constant integer32 := 33;

  begin
   -- put("Give the dimension : "); get(dim);
    declare
      table : Octo_Double_VecVecs.VecVec(0..dim);
      tblix : Standard_Integer_Vectors.Vector(0..dim);
      seq2 : Octo_Double_Vectors.Vector(0..dim);
      two : constant octo_double := Octo_Double_Numbers.create(2.0);
      exa2 : constant octo_double := Octo_Double_Numbers.log(two);
      seq3 : Octo_Double_Vectors.Vector(0..dim);
      six : constant octo_double := Octo_Double_Numbers.create(6.0);
      qdpi : constant octo_double := Octo_Double_Constants.pi;
      exa3 : constant octo_double := (qdpi*qdpi)/six;
    begin
      put_line("running example 2 ...");
      Octo_Double_Example2(seq2);
      Octo_Double_Theta_Algorithm.Allocate(table,tblix,dim);
      Octo_Double_Run_Theta(table,tblix,seq2,exa2,false);
      put("log(2) : "); put(exa2); new_line;
      put_line("running example 3 ...");
      Octo_Double_Example3(seq3);
      Octo_Double_Theta_Algorithm.Allocate(table,tblix,dim);
      Octo_Double_Run_Theta(table,tblix,seq3,exa3,false);
      put("pi^2/6 : "); put(exa3); new_line;
    end;
  end Octo_Double_Main;

  procedure Main is
  begin
    Double_Main;
    Double_Double_Main;
    Quad_Double_Main;
    Octo_Double_Main;
  end Main;

end Test_Theta_Algorithm;
