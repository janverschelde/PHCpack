with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Double_Laurent_Series;

package body Double_Linear_Laurent_Solvers is

  procedure Allocate_Series_Coefficients
              ( dim,deg : in integer32;
                cff : out Standard_Complex_VecVecs.Link_to_VecVec ) is

    res : Standard_Complex_VecVecs.VecVec(1..dim);
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);

  begin
    for i in 1..dim loop
      declare
        val : constant Standard_Complex_Vectors.Vector(0..deg)
            := (0..deg => zero);
      begin
        res(i) := new Standard_Complex_Vectors.Vector'(val);
      end;
    end loop;
    cff := new Standard_Complex_VecVecs.VecVec'(res);
  end Allocate_Series_Coefficients;

  procedure Write ( e : in Standard_Integer_Matrices.Matrix;
                    c : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                    s : in string := "A" ) is
  begin
    for i in e'range(1) loop
      for j in e'range(2) loop
        put(s & "("); put(i,1); put(","); put(j,1); put_line(") :");
        Double_Laurent_Series.Write(e(i,j),c(i)(j).all);
      end loop;
    end loop;
  end Write;

  procedure Write ( e : in Standard_Integer_Vectors.Vector;
                    c : in Standard_Complex_VecVecs.Link_to_VecVec;
                    s : in string := "v" ) is
  begin
    Write(standard_output,e,c,s);
  end Write;

  procedure Write ( file : in file_type;
                    e : in Standard_Integer_Vectors.Vector;
                    c : in Standard_Complex_VecVecs.Link_to_VecVec;
                    s : in string := "v" ) is
  begin
    for i in e'range loop
      put(file,s & "("); put(file,i,1); put_line(file,") :");
      Double_Laurent_Series.Write(file,e(i),c(i).all);
    end loop;
  end Write;

  procedure Matrix_Vector_Product
              ( d : in integer32;
                eA : in Standard_Integer_Matrices.Matrix;
                cA : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                ex : in Standard_Integer_Vectors.Vector;
                cx : in Standard_Complex_VecVecs.Link_to_VecVec;
                ey : out Standard_Integer_Vectors.Vector;
                cy : out Standard_Complex_VecVecs.Link_to_VecVec ) is

    ze,ewrk : integer32;
    zc,cwrk : Standard_Complex_Vectors.Vector(0..d);

  begin
    for i in eA'range(1) loop
      declare
        cAi : constant Standard_Complex_VecVecs.Link_to_VecVec := cA(i);
        cyi : constant Standard_Complex_Vectors.Link_to_Vector := cy(i);
      begin -- initialize with first product instead of with zero
        Double_Laurent_Series.Multiply
          (d,eA(i,eA'first(2)),ex(ex'first),cAi(cAi'first).all,
           cx(cx'first).all,ey(i),cyi.all);
        for j in eA'first(2)+1..eA'last(2) loop
          Double_Laurent_Series.Multiply
            (d,eA(i,j),ex(j),cAi(j).all,cx(j).all,ze,zc);
          Double_Laurent_Series.Add(d,ey(i),ze,cyi.all,zc,ewrk,cwrk);
          ey(i) := ewrk;
          for k in 0..d loop
            cyi(k) := cwrk(k);
          end loop;
        end loop;
      end;
    end loop;
  end Matrix_Vector_Product;

  procedure Forward_Substitution
              ( d : in integer32;
                eL : in Standard_Integer_Matrices.Matrix;
                cL : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                eb : in Standard_Integer_Vectors.Vector;
                cb : in Standard_Complex_VecVecs.Link_to_VecVec;
                ex : out Standard_Integer_Vectors.Vector;
                cx : out Standard_Complex_VecVecs.Link_to_VecVec ) is

    ze,ewrk : integer32;
    zc,cwrk : Standard_Complex_Vectors.Vector(0..d);

  begin
    for i in eb'range loop
      ex(i) := eb(i);
      declare
        cbi : constant Standard_Complex_Vectors.Link_to_Vector := cb(i);
        cLi : constant Standard_Complex_VecVecs.Link_to_VecVec := cL(i);
        cxi : constant Standard_Complex_Vectors.Link_to_Vector := cx(i);
      begin
        for k in 0..d loop
          cxi(k) := cbi(k);
        end loop;
        for j in ex'first..(i-1) loop
          Double_Laurent_Series.Multiply
            (d,eL(i,j),ex(j),cLi(j).all,cx(j).all,ze,zc);
         -- put("L("); put(i,1); put(","); put(j,1);
         -- put(")*x("); put(j,1); put_line(") :");
         -- Double_Laurent_Series.Write(ze,zc);
          Double_Laurent_Series.Subtract(d,ex(i),ze,cxi.all,zc,ewrk,cwrk);
          ex(i) := ewrk;
          for k in 0..d loop
            cxi(k) := cwrk(k);
          end loop;
         -- put("x("); put(i,1); put_line(") after the update :");
         -- Double_Laurent_Series.Write(ex(i),cxi.all);
        end loop;
      end;
    end loop;
  end Forward_Substitution;

  procedure Backward_Substitution
              ( d : in integer32;
                eU : in Standard_Integer_Matrices.Matrix;
                cU : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                eb : in Standard_Integer_Vectors.Vector;
                cb : in Standard_Complex_VecVecs.Link_to_VecVec;
                ex : out Standard_Integer_Vectors.Vector;
                cx : out Standard_Complex_VecVecs.Link_to_VecVec ) is

    ze,ewrk : integer32;
    zc,cwrk : Standard_Complex_Vectors.Vector(0..d);
    cbi,cxi : Standard_Complex_Vectors.Link_to_Vector;
    cUi : Standard_Complex_VecVecs.Link_to_VecVec;

  begin
    cbi := cb(cb'last); cxi := cx(cx'last); cUi := cU(cU'last);
    Double_Laurent_Series.Divide
      (d,eb(eb'last),eU(eU'last(1),eU'last(2)),cbi.all,cUi(cUi'last).all,
       ex(ex'last),cxi.all,cwrk);
    for i in reverse eb'first..eb'last-1 loop
      ex(i) := eb(i);
      cbi := cb(i); cxi := cx(i); cUi := cU(i);
      for k in 0..d loop
        cxi(k) := cbi(k);
      end loop;
      for j in (i+1)..ex'last loop
        Double_Laurent_Series.Multiply
          (d,eU(i,j),ex(j),cUi(j).all,cx(j).all,ze,zc);
       -- put("U("); put(i,1); put(","); put(j,1);
       -- put(")*x("); put(j,1); put_line(") :");
       -- Double_Laurent_Series.Write(ze,zc);
        Double_Laurent_Series.Subtract(d,ex(i),ze,cxi.all,zc,ewrk,cwrk);
        ex(i) := ewrk;
        for k in 0..d loop
          cxi(k) := cwrk(k);
        end loop;
       -- put("x("); put(i,1); put_line(") after the update :");
       -- Double_Laurent_Series.Write(ex(i),cxi.all);
      end loop;
      Double_Laurent_Series.Divide
        (d,ex(i),eU(i,i),cxi.all,cUi(i).all,ze,zc,cwrk);
      ex(i) := ze;
      for k in 0..d loop
        cxi(k) := zc(k);
      end loop;
    end loop;
  end Backward_Substitution;

  function Pivot_Row
              ( nrows,row : in integer32;
                lead : in Standard_Integer_Matrices.Matrix;
                column : in Standard_Complex_VecVecs.Link_to_VecVec )
              return integer32 is

    res : integer32 := row;
    leadval : integer32 := lead(row,row);
    cff : Standard_Complex_Vectors.Link_to_Vector := column(row);
    themax : double_float := Standard_Complex_Numbers.AbsVal(cff(0));
    valmax : double_float;

  begin
    for i in (row+1)..nrows loop
      cff := column(i);
      if lead(i,row) > leadval then
        null; -- higher leading exponents cannot be pivots
      else
        valmax := Standard_Complex_Numbers.AbsVal(cff(0));
        if lead(i,row) < leadval then -- lower exponents must be pivots
          themax := valmax;
          leadval := lead(i,row);
          res := i;
        else -- if lead(i,row) = leadval, then compare coefficients
          if valmax > themax then
            themax := valmax;
            res := i;
          end if;
        end if;
      end if;
    end loop;
    return res;
  end Pivot_Row;

  procedure Swap_Rows
              ( ncols,row,pivrow : in integer32;
                lead : in out Standard_Integer_Matrices.Matrix;
                cffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                pivots : in out Standard_Integer_Vectors.Vector ) is

    itmp : integer32;
    vtmp : Standard_Complex_VecVecs.Link_to_VecVec;
  
  begin
    itmp := pivots(row);
    pivots(row) := pivots(pivrow);
    pivots(pivrow) := itmp;
    vtmp := cffs(row);
    cffs(row) := cffs(pivrow);
    cffs(pivrow) := vtmp;
    for j in 1..ncols loop
      itmp := lead(row,j);
      lead(row,j) := lead(pivrow,j);
      lead(pivrow,j) := itmp;
    end loop;
  end Swap_Rows;

  procedure LU_Factorization
              ( nrows,ncols,deg : in integer32;
                Alead : in out Standard_Integer_Matrices.Matrix;
                Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                pivots : out Standard_Integer_Vectors.Vector ) is

    irow,jrow,Acol : Standard_Complex_VecVecs.Link_to_VecVec;
    icff,jcff : Standard_Complex_Vectors.Link_to_Vector;
    ze,ewrk,eprd,idx : integer32;
    zc,cwrk,cprd : Standard_Complex_Vectors.Vector(0..deg);

  begin
    for i in pivots'range loop
      pivots(i) := i;
    end loop;
    for j in 1..ncols loop
      Acol := Acffs(j);
      idx := Pivot_Row(nrows,j,Alead,Acol);
      if idx /= j
       then Swap_Rows(ncols,j,idx,Alead,Acffs,pivots);
      end if;
      for i in (j+1)..nrows loop
        irow := Acffs(i); icff := irow(j); -- icff is A(i,j)
        jrow := Acffs(j); jcff := jrow(j); -- jcff is A(j,j)
        Double_Laurent_Series.Divide
          (deg,Alead(i,j),Alead(j,j),icff.all,jcff.all,ze,zc,cwrk);
        Alead(i,j) := ze;    -- leading exponent of A(i,j)/A(j,j)
        for k in 0..deg loop -- zc has coefficients of A(i,j)/A(j,j)
          icff(k) := zc(k);  -- A(i,j) := A(i,j)/A(j,j)
        end loop;
        ewrk := 0;
        for k in 0..deg loop
          cwrk(k) := Standard_Complex_Numbers.Create(0.0);
        end loop;
        for k in (j+1)..nrows loop -- A(i,k) := A(i,k) - A(i,j)*A(j,k)
          irow := Acffs(i); icff := irow(j); -- icff is A(i,j)
          jrow := Acffs(j); jcff := jrow(k); -- jcff is A(j,k)
          Double_Laurent_Series.Multiply
            (deg,Alead(i,j),Alead(j,k),icff.all,jcff.all,eprd,cprd);
         -- eprd is the leading exponent of A(i,j)*A(j,k)
         -- cprd has the coefficients of A(i,j)*A(j,k)
          icff := irow(k); -- icff is A(i,k)
          ewrk := Alead(i,k);
          for L in 0..deg loop
            cwrk(L) := icff(L);
          end loop;
          Double_Laurent_Series.Subtract(deg,ewrk,eprd,cwrk,cprd,ze,zc);
          Alead(i,k) := ze;
          for L in 0..deg loop
            icff(L) := zc(L);
          end loop;
        end loop;
      end loop;
    end loop;
  end LU_Factorization;

end Double_Linear_Laurent_Solvers;
