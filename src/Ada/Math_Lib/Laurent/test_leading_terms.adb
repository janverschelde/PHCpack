with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Standard_Integer_VecVecs;
with Standard_Floating_Matrices_io;     use Standard_Floating_Matrices_io;
with Standard_Random_Numbers;
with Standard_Random_Vectors;
with Standard_Random_Matrices;
with Double_Puiseux_Operations;

package body Test_Leading_Terms is

  procedure Sort ( B : in out Standard_Floating_Matrices.Matrix;
                   cB : in out Standard_Complex_Matrices.Matrix;
                   nbrcols : in integer32 ) is

    minval : double_float;
    minidx : integer32;
    nbr : Complex_Number;

  begin
    for i in B'range(1) loop              -- sort i-th row
      for j in B'first(2)..nbrcols-1 loop -- select minimum
        minval := B(i,j);
        minidx := j;
        for k in j+1..nbrcols loop
          if B(i,k) < minval then
            minval := B(i,k);
            minidx := k;
          end if;
        end loop;
        if minidx /= j then -- swap if minimum not at position j
          B(i,minidx) := B(i,j);
          B(i,j) := minval;
          nbr := cB(i,minidx);  
          cB(i,minidx) := cB(i,j);
          cB(i,j) := nbr;
        end if;
      end loop;
    end loop;
  end Sort;

  procedure Series_Product
              ( A : in Standard_Floating_Matrices.Matrix;
                x : in Standard_Floating_Vectors.Vector;
                cA : in Standard_Complex_Matrices.Matrix;
                cx : in Standard_Complex_Vectors.Vector;
                B : out Standard_Floating_Matrices.Matrix;
                cB : out Standard_Complex_Matrices.Matrix ) is
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        B(i,j) := A(i,j) + x(j);
        cB(i,j) := cA(i,j)*cx(j);
      end loop;
    end loop;
    Sort(B,cB,B'last(2));
  end Series_Product;

  procedure Series_Product
              ( A : in Standard_Floating_Matrices.Matrix;
                x : in Standard_Floating_Vectors.Vector;
                cA : in Standard_Complex_Matrices.Matrix;
                cx : in Standard_Complex_Vectors.Vector;
                skip : in Boolean_Vectors.Vector;
                B : out Standard_Floating_Matrices.Matrix;
                cB : out Standard_Complex_Matrices.Matrix;
                nbrcols : out integer32 ) is

    colidx : integer32 := 0;

  begin
    for j in A'range(2) loop
      if not skip(j) then
        colidx := colidx + 1;
        for i in A'range(1) loop
          B(i,colidx) := A(i,j) + x(j);
          cB(i,colidx) := cA(i,j)*cx(j);
        end loop;
      end if;
    end loop;
    nbrcols := colidx;
    Sort(B,cB,nbrcols);
  end Series_Product;

  procedure Series_Product
              ( A,X : in Standard_Floating_Matrices.Matrix;
                cA,cX : in Standard_Complex_Matrices.Matrix;
                B : out Standard_Floating_Matrices.Matrix;
                cB : out Standard_Complex_Matrices.Matrix ) is

    offset : integer32 := 0; -- column offset for B

  begin
    for k in cX'range(2) loop -- run over the columns of X
      for i in A'range(1) loop
        for j in A'range(2) loop
          B(i,j+offset) := A(i,j) + X(j,k);
          cB(i,j+offset) := cA(i,j)*cX(j,k);
        end loop;
      end loop;
      offset := offset + A'last(1);
    end loop;
    Sort(B,cB,B'last(2));
  end Series_Product;

  procedure Series_Product
              ( A,X : in Standard_Floating_Matrices.Matrix;
                cA,cX : in Standard_Complex_Matrices.Matrix;
                skip : in Standard_Integer_Vectors.Vector;
                B : out Standard_Floating_Matrices.Matrix;
                cB : out Standard_Complex_Matrices.Matrix;
                nbrcols : out integer32 ) is

    colidx : integer32 := 0;
    xidx : integer32;

  begin
    for j in X'range(1) loop -- run over the rows of X
      xidx := skip(j) + 1;
      if xidx <= X'last(2) then -- otherwise skip element
        colidx := colidx + 1;
        for i in A'range(1) loop
          B(i,colidx) := A(i,j) + X(j,xidx);
          cB(i,colidx) := cA(i,j)*cX(j,xidx);
        end loop;
      end if;
    end loop;
    nbrcols := colidx;
    Sort(B,cB,nbrcols);
  end Series_Product;

  procedure Random_Vector
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
  end Random_Vector;

  procedure Random_Series
              ( dim,nbr : in integer32;
                X : out Standard_Floating_Matrices.Matrix;
                cX : out Standard_Complex_Matrices.Matrix ) is

    nrows : constant natural32 := natural32(dim);
    ncols : constant natural32 := natural32(nbr);
    rnd : double_float;

  begin
    cX := Standard_Random_Matrices.Random_Matrix(nrows,ncols);
    for i in 1..dim loop
      rnd := Standard_Random_Numbers.Random;
      X(i,1) := abs(rnd);
      for j in 2..nbr loop
        rnd := Standard_Random_Numbers.Random;
        X(i,j) := X(i,j-1) + 0.1 + abs(rnd);
      end loop;
    end loop;
  end Random_Series;

  procedure Random_Series_System
              ( dim,nbr : in integer32;
                A,X,B : out Standard_Floating_Matrices.Matrix;
                cA,cX,cB : out Standard_Complex_Matrices.Matrix ) is

    n : constant natural32 := natural32(dim);

  begin
    A := Standard_Random_Matrices.Random_Matrix(n,n);
    cA := Standard_Random_Matrices.Random_Matrix(n,n);
    for i in A'range(1) loop
      for j in A'range(2) loop
        A(i,j) := abs(A(i,j));
      end loop;
    end loop;
    put("A random "); put(dim,1);
    put_line("-dimensional power matrix :"); put(A,3);
    Random_Series(dim,nbr,X,cX);
    put("A random "); put(dim,1);
    put_line("-dimensional series exponents :"); put(X,3);
    Series_Product(A,X,cA,cX,B,cB);
    put_line("Right hand side exponents :"); put(B,3);
  end Random_Series_System;

  procedure Next_Coefficients
              ( cy : in out Standard_Complex_Vectors.Vector;
                cA,cB : in Standard_Complex_Matrices.Matrix;
                idx1 : in Standard_Integer_Vectors.Vector;
                prev,next : in Boolean_Vectors.Vector;
                cx : in Standard_Complex_Vectors.Vector ) is

    rowidx : integer32;
    err : double_float;

  begin
    for i in next'range loop
      if not prev(i) and next(i) then -- found new correct power
        rowidx := 0;
        for j in idx1'range loop -- look for row index
          if idx1(j) = i
           then rowidx := j; exit;
          end if;
        end loop;
        cy(i) := cB(rowidx,cB'first(2))/cA(rowidx,i);
        put("cy("); put(i,1); put(") : "); put(cy(i)); new_line;
        put("cx("); put(i,1); put(") : "); put(cx(i));
        err := AbsVal(cy(i) - cx(i));
        put(", err :"); put(err,3); new_line;
      end if;
    end loop;
  end Next_Coefficients;

  procedure Next_Series_Coefficients
              ( cy : in out Standard_Complex_Vectors.Vector;
                cA,cB : in Standard_Complex_Matrices.Matrix;
                idx1 : in Standard_Integer_Vectors.Vector;
                next : in Boolean_Vectors.Vector;
                correct : in out Standard_Integer_Vectors.Vector;
                cx : in Standard_Complex_Vectors.Vector;
                cZ : in out Standard_Complex_Matrices.Matrix ) is

    rowidx : integer32;
    err : double_float;

  begin
    for i in next'range loop
      if next(i) then -- found new correct power
        rowidx := 0;
        for j in idx1'range loop -- look for row index
          if idx1(j) = i
           then rowidx := j; exit;
          end if;
        end loop;
        cy(i) := cB(rowidx,cB'first(2))/cA(rowidx,i);
        put("cy("); put(i,1); put(") : "); put(cy(i)); new_line;
        put("cx("); put(i,1); put(") : "); put(cx(i));
        err := AbsVal(cy(i) - cx(i));
        put(", err :"); put(err,3); new_line;
        correct(i) := correct(i) + 1;
        cZ(i,correct(i)) := cy(i);
      end if;
    end loop;
  end Next_Series_Coefficients;

  procedure Coefficient_Check
              ( tol : in double_float;
                cX,cY : in Standard_Complex_Vectors.Vector ) is

    err : double_float;
    sumerr : double_float := 0.0;

  begin
    for i in cY'range loop
      put("cY("); put(i,1); put(") : "); put(cy(i)); new_line;
      put("cX("); put(i,1); put(") : "); put(cx(i));
      err := AbsVal(cY(i) - cX(i));
      put(", err :"); put(err,3); new_line;
      sumerr := sumerr + err;
    end loop;
    put("Sum of errors :"); put(sumerr,3);
    if sumerr < tol
     then put_line(", okay.");
     else put_line(", failure!");
    end if;
  end Coefficient_Check;

  procedure Coefficient_Check
              ( tol : in double_float;
                cX,cY : in Standard_Complex_Matrices.matrix ) is

    err : double_float;
    sumerr : double_float := 0.0;

  begin
    for i in cY'range(1) loop
      for j in cY'range(2) loop
        put("cY("); put(i,1); put(","); put(j,1); put(") : ");
        put(cY(i,j)); new_line;
        put("cX("); put(i,1); put(","); put(j,1); put(") : ");
        put(cX(i,j));
        err := AbsVal(cY(i,j) - cX(i,j));
        put(", err :"); put(err,3); new_line;
        sumerr := sumerr + err;
      end loop;
    end loop;
    put("Sum of errors :"); put(sumerr,3);
    if sumerr < tol
     then put_line(", okay.");
     else put_line(", failure!");
    end if;
  end Coefficient_Check;

  procedure Power_Check
              ( tol : in double_float;
                eX,eY : in Standard_Floating_Matrices.matrix ) is

    err : double_float;
    sumerr : double_float := 0.0;

  begin
    for i in eY'range(1) loop
      for j in eY'range(2) loop
        put("eY("); put(i,1); put(","); put(j,1); put(") : "); put(eY(i,j));
        new_line;
        put("eX("); put(i,1); put(","); put(j,1); put(") : "); put(eX(i,j));
        err := abs(eY(i,j) - eX(i,j));
        put(", err :"); put(err,3); new_line;
        sumerr := sumerr + err;
      end loop;
    end loop;
    put("Sum of errors :"); put(sumerr,3);
    if sumerr < tol
     then put_line(", okay.");
     else put_line(", failure!");
    end if;
  end Power_Check;

  procedure Test_Random_Vector ( dim : in integer32 ) is

    rA,rB : Standard_Floating_Matrices.Matrix(1..dim,1..dim);
    cA,cB : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    rx,vb,ry,rz : Standard_Floating_Vectors.Vector(1..dim);
    cx,cy : Standard_Complex_Vectors.Vector(1..dim);
    m : Standard_Integer_VecVecs.VecVec(0..dim);
    idx1,idx2 : Standard_Integer_Vectors.Vector(1..dim);
    fail,done : boolean;
    prev,next : Boolean_Vectors.Vector(1..dim) := (1..dim => false);
    nbrcols : integer32;
    tol : constant double_float := 1.0E-12;

  begin
    put_line("-> generating random data ...");
    Random_Vector(dim,rA,rB,rx,cA,cB,cx);
    for i in m'range loop
      m(i) := new Standard_Integer_Vectors.Vector'(1..dim => 0);
    end loop;
    for step in 1..dim loop
      put("*** running step "); put(step,1); put_line(" ***");
      for i in vb'range loop
        vb(i) := rB(i,rB'first(2));
      end loop;
      put_line("-> computing the tropical Cramer vector ...");
      Double_Puiseux_Operations.Leading_Powers
        (dim,tol,rA,vb,m,ry,idx1,idx2,fail,2);
      Double_Puiseux_Operations.Check_Correctness
        (dim,tol,rx,ry,idx1,idx2,next,rz,1);
      Next_Coefficients(cy,cA,cB,idx1,prev,next,cx);
      done := true;
      for i in next'range loop
        done := done and next(i);
      end loop;
      if done then
        put("At step "); put(step,1);
        put_line(", all values are correct, done!"); exit;
      end if;
      prev := next; -- for the next round
      Series_Product(rA,rx,cA,cx,next,rB,cB,nbrcols);
    end loop;
    put_line("-> checking final coefficients :");
    Coefficient_Check(tol,cx,cy);
  end Test_Random_Vector;

  procedure Test_Random_Series ( dim,nbr : in integer32 ) is

    prd : constant integer32 := dim*nbr;
    eA : Standard_Floating_Matrices.Matrix(1..dim,1..dim);
    eX,eZ : Standard_Floating_Matrices.Matrix(1..dim,1..nbr);
    eB : Standard_Floating_Matrices.Matrix(1..dim,1..prd);
    cA : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    cX,cZ : Standard_Complex_Matrices.Matrix(1..dim,1..nbr);
    cB : Standard_Complex_Matrices.Matrix(1..dim,1..prd);
    vB,ry,rX,rz : Standard_Floating_Vectors.Vector(1..dim);
    cy,vX : Standard_Complex_Vectors.Vector(1..dim);
    m : Standard_Integer_VecVecs.VecVec(0..dim);
    idx1,idx2,correct : Standard_Integer_Vectors.Vector(1..dim);
    next : Boolean_Vectors.Vector(1..dim) := (1..dim => false);
    fail,done : boolean;
    tol : constant double_float := 1.0E-12;
    nbrcols : integer32;

  begin
    put_line("-> generating random data ...");
    Random_Series_System(dim,nbr,eA,eX,eB,cA,cX,cB);
    for i in m'range loop
      m(i) := new Standard_Integer_Vectors.Vector'(1..dim => 0);
    end loop;
    correct := (1..dim => 0); -- indices in X that are correct
    for step in 1..prd loop
      put("*** running step "); put(step,1); put_line(" ***");
      for i in 1..dim loop
        vB(i) := eB(i,eB'first(2));
        if correct(i) < eX'last(2) then
          rX(i) := eX(i,correct(i)+1);
          vX(i) := cX(i,correct(i)+1);
        else
          rX(i) := 1.0E+99;
        end if;
      end loop;
      put_line("-> computing the tropical Cramer vector ...");
      Double_Puiseux_Operations.Leading_Powers
        (dim,tol,eA,vB,m,rY,idx1,idx2,fail,2);
      Double_Puiseux_Operations.Check_Correctness
        (dim,tol,rX,rY,idx1,idx2,next,rz,1);
      Next_Series_Coefficients(cy,cA,cB,idx1,next,correct,vX,cZ);
      for i in next'range loop
        if next(i)
         then eZ(i,correct(i)) := rY(i);
        end if;
      end loop;
      put("correct indices :"); put(correct); new_line;
      done := true;
      for i in next'range loop
        done := done and (correct(i) = eX'last(2));
      end loop;
      if done then
        put("At step "); put(step,1);
        put_line(", all values are correct, done!"); exit;
      end if;
      next := (1..dim => false);
      Series_Product(eA,eX,cA,cX,correct,eB,cB,nbrcols);
    end loop;
    put_line("-> checking final powers :");
    Power_Check(tol,eX,eZ);
    put_line("-> checking final coefficients :");
    Coefficient_Check(tol,cX,cZ);
  end Test_Random_Series;

  procedure Main is

    dim,nbr : integer32 := 0;

  begin
    new_line;
    put_line("Testing leading terms computation ...");
    put("-> Give the dimension : "); get(dim);
    put("-> Give the number of terms : "); get(nbr);
    if nbr = 1
     then Test_Random_Vector(dim);
     else Test_Random_Series(dim,nbr);
    end if;
  end Main;

end Test_Leading_Terms;
