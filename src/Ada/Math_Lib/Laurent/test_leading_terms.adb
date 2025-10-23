with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Boolean_Vectors;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Matrices_io;     use Standard_Floating_Matrices_io;
with Standard_Random_Vectors;
with Standard_Random_Matrices;
with Double_Puiseux_Operations;

package body Test_Leading_Terms is

  procedure Series_Product
              ( A : in Standard_Floating_Matrices.Matrix;
                x : in Standard_Floating_Vectors.Vector;
                cA : in Standard_Complex_Matrices.Matrix;
                cx : in Standard_Complex_Vectors.Vector;
                B : out Standard_Floating_Matrices.Matrix;
                cB : out Standard_Complex_Matrices.Matrix ) is

    minb : double_float;
    minidx : integer32;
    nbr : Complex_Number;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        B(i,j) := A(i,j) + x(j);
        cB(i,j) := cA(i,j)*cx(j);
      end loop;
    end loop;
    for i in B'range(1) loop                -- sort i-th row
      for j in B'first(2)..B'last(2)-1 loop -- select minimum
        minb := B(i,j); minidx := j;
        for k in j+1..B'last(2) loop
          if B(i,k) < minb
           then minb := B(i,k); minidx := k;
          end if;
        end loop;
        if minidx /= j then -- swap if minimum not at position j
          B(i,minidx) := B(i,j); B(i,j) := minb;
          nbr := cB(i,minidx);  
          cB(i,minidx) := cB(i,j); cB(i,j) := nbr;
        end if;
      end loop;
    end loop;
  end Series_Product;

  procedure Random_Input
              ( dim : in integer32;
                A,B : out Standard_Floating_Matrices.Matrix;
                x : out Standard_Floating_Vectors.Vector;
                cA,cB : out Standard_Complex_Matrices.Matrix;
                cx : out Standard_Complex_Vectors.Vector ) is

    n : constant natural32 := natural32(dim);

  begin
    A := Standard_Random_Matrices.Random_Matrix(n,n);
    cA := Standard_Random_Matrices.Random_Matrix(n,n);
    x := Standard_Random_Vectors.Random_Vector(1,dim);
    cx := Standard_Random_Vectors.Random_Vector(1,dim);
    for i in 1..dim loop
      x(i) := abs(x(i));
      for j in A'range(2) loop
        A(i,j) := abs(A(i,j));
      end loop;
    end loop;
    put("A random "); put(dim,1);
    put_line("-dimensional power matrix :"); put(A,3);
    put("A random "); put(dim,1);
    put_line("-dimensional power vector :"); put(x,3); new_line;
    Series_Product(A,x,cA,cx,B,cB);
    put_line("Power matrix after product :"); put(B,3);
  end Random_Input;

  procedure Test_Random_Input ( dim : in integer32 ) is

    rA,rB : Standard_Floating_Matrices.Matrix(1..dim,1..dim);
    cA,cB : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    rx,vb,ry,rz : Standard_Floating_Vectors.Vector(1..dim);
    cx,cy : Standard_Complex_Vectors.Vector(1..dim);
    m : Standard_Integer_VecVecs.VecVec(0..dim);
    idx1,idx2 : Standard_Integer_Vectors.Vector(1..dim);
    fail : boolean;
    prev,next : Boolean_Vectors.Vector(1..dim) := (1..dim => false);
    rowidx : integer32;
    err : double_float;

  begin
    put_line("-> generating random data ...");
    Random_Input(dim,rA,rB,rx,cA,cB,cx);
    put_line("-> computing the tropical Cramer vector ...");
    for i in m'range loop
      m(i) := new Standard_Integer_Vectors.Vector'(1..dim => 0);
    end loop;
    for i in vb'range loop
      vb(i) := rB(i,rB'first(2));
    end loop;
    Double_Puiseux_Operations.Leading_Powers
      (dim,1.0E-12,rA,vb,m,ry,idx1,idx2,fail,2);
    Double_Puiseux_Operations.Check_Correctness
      (dim,1.0E-12,rx,ry,idx1,idx2,next,rz,1);
    for i in next'range loop
      if not prev(i) and next(i) then -- found new correct power
        rowidx := 0;
        for j in idx1'range loop -- look for row index
          if idx1(j) = i
           then rowidx := j; exit;
          end if;
        end loop;
        cy(i) := cB(rowidx,cB'first(2))/cA(rowidx,i);
        put("cy("); put(i,1); put(") :"); put(cy(i)); new_line;
        put("cx("); put(i,1); put(") :"); put(cx(i));
        err := AbsVal(cy(i) - cx(i));
        put(", err :"); put(err,3); new_line;
      end if;
    end loop;
    prev := next; -- for the next round
  end Test_Random_Input;

  procedure Main is

    dim : integer32 := 0;

  begin
    new_line;
    put_line("Testing leading terms computation ...");
    put("-> Give the dimension : "); get(dim);
    Test_Random_Input(dim);
  end Main;

end Test_Leading_Terms;
