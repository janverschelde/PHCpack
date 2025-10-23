with Ada.Text_IO;                       use Ada.Text_IO;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Standard_Random_Numbers;
with Standard_Random_Vectors;
with Standard_Random_Matrices;
with Standard_Floating_Matrices_io;     use Standard_Floating_Matrices_io;
with Double_Puiseux_Operations;

package body Test_Leading_Powers is

  function Row_Min_Plus
             ( A : Matrix; x : Standard_Floating_Vectors.Vector )
             return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(A'range(1));
    val : double_float;

  begin
    for i in res'range loop
      res(i) := A(i,A'first(2)) + x(x'first);
      for j in A'first(2)+1..A'last(2) loop
        val := A(i,j) + x(j);
        if val < res(i)
         then res(i) := val;
        end if;
      end loop;
    end loop;
    return res;
  end Row_Min_Plus;

  function Row_Min_Plus
             ( A : Matrix; x : Standard_Floating_Vectors.Vector;
               skip : Boolean_Vectors.Vector )
             return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(A'range(1));
    val : double_float;

  begin
    for i in res'range loop
      res(i) := 1.0E+99;
      for j in A'range(2) loop
        if not skip(j) then
          val := A(i,j) + x(j);
          if val < res(i)
           then res(i) := val;
          end if;
        end if;
      end loop;
    end loop;
    return res;
  end Row_Min_Plus;

  function id ( n : integer32 ) return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..n);

  begin
    for i in 1..n loop
      res(i) := i;
    end loop;
    return res;
  end id;

  procedure Shuffle ( p : in out Standard_Integer_Vectors.Vector ) is

    dim : constant integer32 := p'last;
    j,v : integer32;

  begin
    for i in p'range loop
      j := Standard_Random_Numbers.Random(1,dim);
      if j /= i then -- swap
        v := p(i);
        p(i) := p(j);
        p(j) := v;
      end if;
    end loop;
  end Shuffle;

  function Random_Permutation
             ( n : integer32 ) return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..n) := id(n);

  begin
    Shuffle(res);
    return res;
  end Random_Permutation;

  procedure Random_Leading_Input 
              ( dim : in integer32;
                p : out Standard_Integer_Vectors.Vector;
                A : out Matrix;
                x,b : out Standard_Floating_Vectors.Vector ) is

    n : constant natural32 := natural32(dim);

  begin
    A := Standard_Random_Matrices.Random_Matrix(n,n);
    x := Standard_Random_Vectors.Random_Vector(1,dim);
    p := Random_Permutation(dim);
    put("A random permutation :"); put(p); new_line;
    for i in 1..dim loop
      x(i) := abs(x(i))/double_float(2*dim);
      for j in A'range(2) loop
        if j = p(i) -- permutation determines on which row minimum occurs
         then A(i,j) := abs(A(i,j));
         else A(i,j) := 1.0 + abs(A(i,j));
        end if;
      end loop;
    end loop;
    put("A random "); put(dim,1);
    put_line("-dimensional matrix :"); put(A,3);
    put("A random "); put(dim,1);
    put_line("-dimensional vector :"); put(x,3); new_line;
    b := Row_Min_Plus(A,x);
    put_line("The right hand side vector : "); put(b,3); new_line;
  end Random_Leading_Input;

  procedure Random_General_Input 
              ( dim : in integer32;
                A : out Matrix;
                x,b : out Standard_Floating_Vectors.Vector ) is

    n : constant natural32 := natural32(dim);

  begin
    A := Standard_Random_Matrices.Random_Matrix(n,n);
    x := Standard_Random_Vectors.Random_Vector(1,dim);
    for i in 1..dim loop
      x(i) := abs(x(i));
      for j in A'range(2) loop
        A(i,j) := abs(A(i,j));
      end loop;
    end loop;
    put("A random "); put(dim,1);
    put_line("-dimensional matrix :"); put(A,3);
    put("A random "); put(dim,1);
    put_line("-dimensional vector :"); put(x,3); new_line;
    b := Row_Min_Plus(A,x);
    put_line("The right hand side vector : "); put(b,3); new_line;
  end Random_General_Input;

  procedure Check_Differences
              ( dim : in integer32; A : in Matrix;
                b,x,d : in Standard_Floating_Vectors.Vector ) is

    r : Standard_Floating_Vectors.Vector(1..dim);
    xmd,bmr : double_float;

  begin
    put_line("The leading degrees d : "); put(d,3); new_line;
    put_line("The original x : "); put(x,3); new_line;
    xmd := 0.0;
    for i in x'range loop
      xmd := xmd + abs(x(i) - d(i));
    end loop;
    put("-> difference |x-d| :"); put(xmd,3); new_line;
    r := Row_Min_Plus(A,d);
    put_line("Comparing right hand sides :");
    bmr := 0.0;
    for i in r'range loop
      put("b("); put(i,1); put(") :"); put(b(i));
      put(" - r("); put(i,1); put(") :"); put(r(i));
      put(" = "); put(b(i)-r(i),2,3,3); new_line;
      bmr := bmr + abs(b(i)-r(i));
    end loop;
    put("-> difference |b-r| :"); put(bmr,3); new_line;
    if xmd + bmr < 1.0e-10
     then put_line("Test succeeded.");
     else put_line("Differences are too high, test failed.");
    end if;
  end Check_Differences;

  procedure Test_Leading_Random ( dim : in integer32 ) is

    A : Matrix(1..dim,1..dim);
    x,b,d : Standard_Floating_Vectors.Vector(1..dim);
    p,idx1,idx2 : Standard_Integer_Vectors.Vector(1..dim);
    m : Standard_Integer_VecVecs.VecVec(0..dim);
    fail : boolean;

  begin
    put_line("-> generating random data ...");
    Random_Leading_Input(dim,p,A,x,b);
    put_line("-> computing the tropical Cramer vector ...");
    for i in m'range loop
      m(i) := new Standard_Integer_Vectors.Vector'(1..dim => 0);
    end loop;
    Double_Puiseux_Operations.Leading_Powers
      (dim,1.0E-12,A,b,m,d,idx1,idx2,fail,2);
    if fail then
      put_line("Unexpected failure reported!");
    else
      fail := false;
      for i in idx2'range loop
        fail := fail or (idx2(i) /= dim+1);
      end loop;
      if not fail then
        put("All second indices equal ");
        put(dim+1,1); put_line(", okay.");
      else
        put("All second indices are not equal to ");
        put(dim+1,1); put_line(", bug!");
      end if;
      put_line("-> checking differences ...");
      Check_Differences(dim,A,b,x,d);
    end if;
    for i in m'range loop
      Standard_Integer_Vectors.Clear(m(i));
    end loop;
  end Test_Leading_Random;

  procedure Test_General_Random ( dim : in integer32 ) is

    A : Matrix(1..dim,1..dim);
    x,b,d,cd : Standard_Floating_Vectors.Vector(1..dim);
    idx1,idx2 : Standard_Integer_Vectors.Vector(1..dim);
    m : Standard_Integer_VecVecs.VecVec(0..dim);
    fail,done : boolean;
    correct : Boolean_Vectors.Vector(1..dim) := (1..dim => false);

  begin
    put_line("-> generating random data ...");
    Random_General_Input(dim,A,x,b);
    cd := (1..dim => -1.0);
    for i in 1..dim loop
      put("*** running step "); put(i,1); put_line(" ***");
      put_line("-> computing the tropical Cramer vector ...");
      for i in m'range loop
        m(i) := new Standard_Integer_Vectors.Vector'(1..dim => 0);
      end loop;
      Double_Puiseux_Operations.Leading_Powers
        (dim,1.0E-12,A,b,m,d,idx1,idx2,fail,2);
      if not fail then
        put_line("No failure reported, unexpected.");
        put_line("-> checking differences ...");
        Check_Differences(dim,A,b,x,d);
      else
        put_line("Failure reported, as expected.");
        Double_Puiseux_Operations.Check_Correctness
          (dim,1.0E-12,x,d,idx1,idx2,correct,cd,1);
        done := true;
        for i in correct'range loop
          done := done and correct(i);
        end loop;
        if done then
          fail := false;
          put("A step "); put(i,1);
          put_line(", all values are correct, done!");
        else
          b := Row_Min_Plus(A,x,correct);
        end if;
      end if;
      exit when not fail;
    end loop;
    for i in m'range loop
      Standard_Integer_Vectors.Clear(m(i));
    end loop;
  end Test_General_Random;

  procedure Main is

    dim : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Testing leading powers computation ...");
    put("-> Test specific random data ? (y/n) ");
    Ask_Yes_or_No(ans);
    put("-> Give the dimension : "); get(dim);
    if ans = 'y'
     then Test_Leading_Random(dim);
     else Test_General_Random(dim);
    end if;
  end Main;

end Test_Leading_Powers;
