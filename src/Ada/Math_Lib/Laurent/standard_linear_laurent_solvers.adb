with Standard_Complex_Vectors;
with Standard_Laurent_Series;

package body Standard_Linear_Laurent_Solvers is

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
        Standard_Laurent_Series.Multiply
          (d,eA(i,eA'first(2)),ex(ex'first),cAi(cAi'first).all,
           cx(cx'first).all,ey(i),cyi.all);
        for j in eA'first(2)+1..eA'last(2) loop
          Standard_Laurent_Series.Multiply
            (d,eA(i,j),ex(j),cAi(j).all,cx(j).all,ze,zc);
          Standard_Laurent_Series.Add(d,ey(i),ze,cyi.all,zc,ewrk,cwrk);
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
          Standard_Laurent_Series.Multiply
            (d,eL(i,j),ex(j),cLi(j).all,cx(j).all,ze,zc);
         -- put("L("); put(i,1); put(","); put(j,1);
         -- put(")*x("); put(j,1); put_line(") :");
         -- Standard_Laurent_Series.Write(ze,zc);
          Standard_Laurent_Series.Subtract(d,ex(i),ze,cxi.all,zc,ewrk,cwrk);
          ex(i) := ewrk;
          for k in 0..d loop
            cxi(k) := cwrk(k);
          end loop;
         -- put("x("); put(i,1); put_line(") after the update :");
         -- Standard_Laurent_Series.Write(ex(i),cxi.all);
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
    Standard_Laurent_Series.Divide
      (d,eb(eb'last),eU(eU'last(1),eU'last(2)),cbi.all,cUi(cUi'last).all,
       ex(ex'last),cxi.all,cwrk);
    for i in reverse eb'first..eb'last-1 loop
      ex(i) := eb(i);
      cbi := cb(i); cxi := cx(i); cUi := cU(i);
      for k in 0..d loop
        cxi(k) := cbi(k);
      end loop;
      for j in (i+1)..ex'last loop
        Standard_Laurent_Series.Multiply
          (d,eU(i,j),ex(j),cUi(j).all,cx(j).all,ze,zc);
       -- put("U("); put(i,1); put(","); put(j,1);
       -- put(")*x("); put(j,1); put_line(") :");
       -- Standard_Laurent_Series.Write(ze,zc);
        Standard_Laurent_Series.Subtract(d,ex(i),ze,cxi.all,zc,ewrk,cwrk);
        ex(i) := ewrk;
        for k in 0..d loop
          cxi(k) := cwrk(k);
        end loop;
       -- put("x("); put(i,1); put_line(") after the update :");
       -- Standard_Laurent_Series.Write(ex(i),cxi.all);
      end loop;
      Standard_Laurent_Series.Divide
        (d,ex(i),eU(i,i),cxi.all,cUi(i).all,ze,zc,cwrk);
      ex(i) := ze;
      for k in 0..d loop
        cxi(k) := zc(k);
      end loop;
    end loop;
  end Backward_Substitution;

end Standard_Linear_Laurent_Solvers;
