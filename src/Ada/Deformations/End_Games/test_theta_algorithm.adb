with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Double_Theta_Algorithm;

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

  begin
    nbr(0) := 1.0;
    for k in 1..nbr'last loop
      cff := 1.0/double_float(k*k);
      nbr(k) := nbr(k-1) + cff;
    end loop;
  end Double_Example3;

  procedure Double_Run_Theta
              ( tab : in Standard_Floating_VecVecs.VecVec;
                idx : in out Standard_Integer_Vectors.Vector;
                seq : in Standard_Floating_Vectors.Vector;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Runs the theta algorithm on the sequence seq.

  -- REQUIRED : tab'last = idx'last.

    dim : integer32 := 0;
 
  begin
    if verbose then
      put("adding "); put(seq(0)); put_line(" ...");
    end if;
    Double_Theta_Algorithm.Initialize(tab,dim,idx,seq(0),verbose);
    if verbose then
      put("adding "); put(seq(1)); put_line(" ...");
    end if;
    Double_Theta_Algorithm.Initialize(tab,dim,idx,seq(1),verbose);
    for k in 2..seq'last loop
      if verbose then
        put("adding "); put(seq(k)); put_line(" ...");
      end if;
      Double_Theta_Algorithm.Extrapolate(tab,dim,idx,seq(k),verbose);
    end loop;
    put("computed "); put(dim,1); put_line(" columns :");
    for col in 0..dim-1 loop
      if col mod 2 = 0 then
        put(tab(col)(idx(col))); new_line;
      end if;
    end loop;
  end Double_Run_Theta;

  procedure Main is

    dim : integer32 := 0;

  begin
    put("Give the dimension : "); get(dim);
    declare
      table : Standard_Floating_VecVecs.VecVec(0..dim);
      tblix : Standard_Integer_Vectors.Vector(0..dim);
      seq2 : Standard_Floating_Vectors.Vector(0..dim);
      seq3 : Standard_Floating_Vectors.Vector(0..dim);
    begin
      put_line("running example 2 ...");
      Double_Example2(seq2);
      Double_Theta_Algorithm.Allocate(table,tblix,dim);
      Double_Run_Theta(table,tblix,seq2);
     -- put_line("running example 3 ...");
     -- Double_Example3(seq3);
     -- Double_Theta_Algorithm.Allocate(table,tblix,dim);
     -- Double_Run_Theta(table,tblix,seq3);
    end;
  end Main;

end Test_Theta_Algorithm;
