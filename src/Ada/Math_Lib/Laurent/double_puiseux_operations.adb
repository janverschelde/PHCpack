with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Standard_Floating_Matrices_io;     use Standard_Floating_Matrices_io;
with Double_Weighted_Assignment;

package body Double_Puiseux_Operations is

  procedure Leading_Powers
              ( dim : in integer32; tol : in double_float;
                A : in Standard_Floating_Matrices.Matrix;
                b : in Standard_Floating_Vectors.Vector;
                wrkm : in Standard_Integer_VecVecs.VecVec;
                d : out Standard_Floating_Vectors.Vector;
                idxone,idxtwo : out Standard_Integer_Vectors.Vector;
                fail : out boolean; vrblvl : in integer32 := 0 ) is

    c : Standard_Floating_Vectors.Vector(0..dim);
    abc : Standard_Floating_Matrices.Matrix(1..dim,1..dim+1);
    mv : Standard_Floating_Vectors.Vector(1..dim);
    idxmv : Standard_Integer_Vectors.Vector(1..2*dim);

  begin
    if vrblvl > 0
     then put_line("-> in Double_Puiseux_Operations.leading_powers ...");
    end if;
    for i in wrkm'range loop
      for j in wrkm(i)'range loop
        wrkm(i)(j) := 0;
      end loop;
    end loop;
    Double_Weighted_Assignment.cramer_vector(A,b,c,wrkm,vrblvl-1);
    abc := Double_Weighted_Assignment.Abc_Matrix(A,b,c);
    if vrblvl > 0 then
      put_line("-> checking if minimum is attained twice ...");
      put_line("A : "); put(A,1,A'last(1),3);
      put_line("b : "); put(b,3); new_line;
      put_line("c : "); put(c,3); new_line;
      put_line("Cramer vector added to A | b :");
      for i in abc'range(1) loop
        for j in abc'first(2)..abc'last(2)-1 loop
          put(abc(i,j),3);
        end loop;
        put(" | ");
        put(abc(i,abc'last(2)),3);
        new_line;
      end loop;
    end if;
    Double_Weighted_Assignment.Abc_ArgMin(abc,tol,mv,idxmv,fail,vrblvl-1);
    if fail and (vrblvl > 0)
     then put_line("Failed to reach minimum exactly twice!");
    end if;
    for i in 1..dim loop
      idxone(i) := idxmv(2*i-1);
      idxtwo(i) := idxmv(2*i);
    end loop;
    if vrblvl > 0 then
      put("1st index :"); put(idxone); new_line;
      put("2nd index :"); put(idxtwo); new_line;
    end if;
    for i in d'range loop
      d(i) := c(i) - c(0);
    end loop;
  end Leading_Powers;

  procedure Check_Correctness
              ( dim : in integer32; tol : in double_float;
                x,d : in Standard_Floating_Vectors.Vector;
                idx1,idx2 : in Standard_Integer_Vectors.Vector;
                correct : out Boolean_Vectors.Vector;
                cd : in out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    cnt : integer32 := 0;
    err : double_float;

  begin
    if vrblvl > 0
     then put_line("-> in Double_Puiseux_Operations.check_correctness ...");
    end if;
    if vrblvl = 0 then
      for i in idx2'range loop
        if idx2(i) = dim+1 then
          err := abs(x(idx1(i)) - d(idx1(i)));
          if err < tol then
            correct(idx1(i)) := true;
            cd(idx1(i)) := d(idx1(i));
          end if;
        end if;
      end loop;
    else
      for i in idx2'range loop
        put("-> value for "); put(idx1(i),1);
        if idx2(i) = dim+1
         then put(" is correct? ");
         else put(" not correct ");
        end if;
        err := abs(x(idx1(i)) - d(idx1(i)));
        put(":"); put(err,3);
        if idx2(i) = dim+1 then
          if err < 1.0e-12 then
            put_line(", okay"); cnt := cnt + 1;
            correct(idx1(i)) := true;
            cd(idx1(i)) := d(idx1(i));
          else
            put_line(", bug!");
          end if;
        else
          if err < 1.0e-12 
           then put_line(", minimum attained 3 times?");
           else put_line(", okay");
          end if;
        end if;
      end loop;
      put("Found "); put(cnt,1);
      if cnt = 1
       then put_line(" correct value.");
       else put_line(" correct values.");
      end if;
      for i in correct'range loop
        if correct(i) then
          put("-> value "); put(i,1); put_line(" is correct :");
          put("d("); put(i,1); put(") :"); put(cd(i)); new_line;
          put("x("); put(i,1); put(") :"); put(x(i));
          err := cd(i) - x(i);
          put(", err :"); put(err,3); new_line;
        end if;
      end loop;
    end if;
  end Check_Correctness;

end Double_Puiseux_Operations;
