with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Numbers_IO;      use Standard_Floating_Numbers_IO;
with Standard_Complex_Numbers_IO;       use Standard_Complex_Numbers_IO;
with Standard_Integer_Vectors_IO;       use Standard_Integer_Vectors_IO;
with Standard_Floating_Vectors_IO;      use Standard_Floating_Vectors_IO;
with Standard_Floating_Matrices_IO;     use Standard_Floating_Matrices_IO;
with Double_Weighted_Assignment;

package body Double_Puiseux_Operations is

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
                cB : out Standard_Complex_Matrices.Matrix;
                tosort : in boolean := false ) is
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        B(i,j) := A(i,j) + x(j);
        cB(i,j) := cA(i,j)*cx(j);
      end loop;
    end loop;
    if tosort
     then Sort(B,cB,B'last(2));
    end if;
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
                cB : out Standard_Complex_Matrices.Matrix;
                tosort : in boolean := false ) is

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
    if tosort
     then Sort(B,cB,B'last(2));
    end if;
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

  procedure Assign_Correctness
              ( dim : in integer32;
                d : in Standard_Floating_Vectors.Vector;
                idx1,idx2 : in Standard_Integer_Vectors.Vector;
                correct : out Boolean_Vectors.Vector;
                cd : in out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    index1 : integer32 := 0;

  begin
    if vrblvl > 0
     then put_line("-> in Double_Puiseux_Operations.assign_correctness ...");
    end if;
    for i in idx2'range loop
      if idx2(i) = dim+1 then
        index1 := idx1(i);
        correct(index1) := true;
        cd(index1) := d(index1);
        if vrblvl > 0 then
          put("-> value "); put(index1,1); put_line(" is correct :");
          put("d("); put(index1,1); put(") :"); put(cd(index1)); new_line;
        end if;
      end if;
    end loop;
  end Assign_Correctness;

  procedure Next_Coefficients
              ( cy : in out Standard_Complex_Vectors.Vector;
                cA : in Standard_Complex_Matrices.Matrix;
                cb : in Standard_Complex_Vectors.Vector;
                idx1 : in Standard_Integer_Vectors.Vector;
                prev,next : in Boolean_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    rowidx : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in Double_Puiseux_Operations.next_coefficients ...");
    end if;
    for i in next'range loop
      if not prev(i) and next(i) then -- found new correct power
        rowidx := 0;
        for j in idx1'range loop -- look for row index
          if idx1(j) = i
           then rowidx := j; exit;
          end if;
        end loop;
        cy(i) := cb(rowidx)/cA(rowidx,i);
        if vrblvl > 0
         then put("cy("); put(i,1); put(") : "); put(cy(i)); new_line;
        end if;
      end if;
    end loop;
  end Next_Coefficients;

  function Leading_Right_Power
             ( rB : Standard_Floating_Matrices.Matrix; rowidx : integer32;
               skipcols : Boolean_Vectors.Vector; vrblvl : integer32 := 0 )
             return double_float is

    res : double_float := -1.0;

  begin
    if vrblvl > 0 then
      put("-> in Double_Puiseux_Operations.");
      put_line("leading_right_power ...");
    end if;
    for k in rB'range(2) loop           -- run over all columns
      if not skipcols(k) then           -- do we skip column k?
        if res = -1.0 then
          res := rB(rowidx,k);          -- first valid column
        elsif rB(rowidx,k) < res then
          res := rB(rowidx,k);          -- new minimal value
        end if;
      end if;
    end loop;
    return res;
  end Leading_Right_Power;

  procedure Leading_Right_Term
             ( rB : in Standard_Floating_Matrices.Matrix;
               cB : in Standard_Complex_Matrices.Matrix;
               rowidx : in integer32;
               skipcols : in Boolean_Vectors.Vector;
               minpow : out double_float; mincff : out Complex_Number;
               vrblvl : in integer32 := 0 ) is

    minidx : integer32 := -1;

  begin
    if vrblvl > 0 then
      put("-> in Double_Puiseux_Operations.");
      put_line("leading_right_term ...");
    end if;
    minpow := -1.0;
    for k in rB'range(2) loop           -- run over all columns
      if not skipcols(k) then           -- do we skip column k?
        if minidx = -1 then
          minidx := k;                  -- first valid column
          minpow := rB(rowidx,minidx);
        elsif rB(rowidx,k) < minpow then
          minidx := k;                  -- new minimal value
          minpow := rB(rowidx,minidx);
        end if;
      end if;
    end loop;
    if minidx = -1 then
      if vrblvl > 0
       then put_line("No leading power found in the right hand side!");
      end if;
    else
      mincff := cB(rowidx,minidx);
    end if;
  end Leading_Right_Term;

  function Leading_Right_Powers
             ( rB : Standard_Floating_Matrices.Matrix;
               skipcols : Boolean_Vectors.Vector; vrblvl : integer32 := 0 )
             return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(rB'range(1));

  begin
    if vrblvl > 0 then
      put("-> in Double_Puiseux_Operations.");
      put_line("leading_right_powers ...");
    end if;
    for i in rB'range(1) loop
      res(i) := Leading_Right_Power(rB,i,skipcols,vrblvl-1);
    end loop;
    return res;
  end Leading_Right_Powers;

  procedure Leading_Right_Terms
             ( rB : in Standard_Floating_Matrices.Matrix;
               cB : in Standard_Complex_Matrices.Matrix;
               skipcols : in Boolean_Vectors.Vector;
               minpow : out Standard_Floating_Vectors.Vector;
               mincff : out Standard_Complex_Vectors.Vector;
               vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put("-> in Double_Puiseux_Operations.");
      put_line("leading_right_terms ...");
    end if;
    for i in rB'range(1) loop
      Leading_Right_Term(rB,cB,i,skipcols,minpow(i),mincff(i),vrblvl-1);
    end loop;
  end Leading_Right_Terms;

  procedure Leading_Solver
              ( dim : in integer32; tol : in double_float;
                rA,rB : in Standard_Floating_Matrices.Matrix;
                cA,cB : in Standard_Complex_Matrices.Matrix;
                ry : out Standard_Floating_Vectors.Vector;
                cy : out Standard_Complex_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    mrb,dy : Standard_Floating_Vectors.Vector(1..dim);
    mcb : Standard_Complex_Vectors.Vector(1..dim);
    wrkcrm : Standard_Integer_VecVecs.VecVec(0..dim);
    idx1,idx2 : Standard_Integer_Vectors.Vector(1..dim);
    fail,done : boolean;
    prev,next : Boolean_Vectors.Vector(1..dim) := (1..dim => false);

  begin
    if vrblvl > 0
     then put_line("-> in Double_Puiseux_Operations.leading_solver ...");
    end if;
    ry := (1..dim => 0.0);
    cy := (1..dim => Standard_Complex_Numbers.create(0.0));
    for i in wrkcrm'range loop
      wrkcrm(i) := new Standard_Integer_Vectors.Vector'(1..dim => 0);
    end loop;
    for step in 1..dim loop
      if vrblvl > 0
       then put("*** running step "); put(step,1); put_line(" ***");
      end if;
      Leading_Right_Terms(rB,cB,next,mrb,mcb,vrblvl-1);
      if vrblvl > 0
       then put_line("-> computing the tropical Cramer vector ...");
      end if;
      Leading_Powers(dim,tol,rA,mrb,wrkcrm,dy,idx1,idx2,fail,vrblvl-1);
      Assign_Correctness(dim,dy,idx1,idx2,next,ry,vrblvl-1);
      Next_Coefficients(cy,cA,mcb,idx1,prev,next,vrblvl-1);
      done := true;
      for i in next'range loop
        done := done and next(i);
      end loop;
      if done then
        if vrblvl > 0 then
          put("At step "); put(step,1);
          put_line(", all values are correct, done!"); exit;
        end if;
      end if;
      prev := next; -- for the next round
    end loop;
  end Leading_Solver;

end Double_Puiseux_Operations;
