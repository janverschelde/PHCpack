with Standard_Complex_Matrices;
with Standard_Vector_Splitters;
with Standard_Matrix_Splitters;
with Standard_Inlined_Linear_Solvers;

package body Standard_Inlined_Linearization is

  procedure Row_Matrix_Multiply
              ( rArows,iArows : in Standard_Floating_VecVecs.Link_to_VecVec;
                rx,ix,ry,iy : in Standard_Floating_Vectors.Link_to_Vector ) is

    rlnk,ilnk : Standard_Floating_Vectors.Link_to_Vector;
    pr,pi,qr,qi : double_float;

  begin
    for k in rArows'range loop
      rlnk := rArows(k); ilnk := iArows(k);
      ry(k) := 0.0; iy(k) := 0.0;
      for j in rx'range loop
        pr := rlnk(j); pi := ilnk(j);
        qr := rx(j);   qi := ix(j);
        ry(k) := ry(k) + pr*qr - pi*qi;
        iy(k) := iy(k) + pr*qi + pi*qr;
      end loop;
    end loop;
  end Row_Matrix_Multiply;

  procedure Inlined_Solve_Head_by_lufac
              ( dim : in integer32; 
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rb,ib : in Standard_Floating_Vectors.Link_to_Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 ) is
  begin
    Standard_Inlined_Linear_Solvers.lufac(rc,ic,dim,ipvt,info);
    if info = 0 then
      Standard_Inlined_Linear_Solvers.lusolve(rc,ic,dim,ipvt,rb,ib);
    end if;
  end Inlined_Solve_Head_by_lufac;

  procedure Inlined_Solve_Head_by_lufco
              ( dim : in integer32; 
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rb,ib : in Standard_Floating_Vectors.Link_to_Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_float;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    Standard_Inlined_Linear_Solvers.lufco(rc,ic,dim,ipvt,ry,iy,rcond);
    if 1.0 + rcond /= 1.0 then
      Standard_Inlined_Linear_Solvers.lusolve(rc,ic,dim,ipvt,rb,ib);
    end if;
  end Inlined_Solve_Head_by_lufco;

  procedure Inlined_Solve_Tail_by_lusolve
              ( deg,dim : in integer32; 
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector;
                idx : in integer32 := 1 ) is

    rlnk,ilnk : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for k in idx..deg loop                         -- loop to compute x(k)
      Row_Matrix_Multiply(rv(k),iv(k),rb(0),ib(0),ry,iy); -- y = A(k)*x(0)
      rlnk := rb(k); ilnk := ib(k);
      for j in rlnk'range loop          -- compute b(k) = b(k) - A(k)*x(0)
        rlnk(j) := rlnk(j) - ry(j);
        ilnk(j) := ilnk(j) - iy(j);
      end loop;
      for i in 1..(k-1) loop
        Row_Matrix_Multiply(rv(k-i),iv(k-i),rb(i),ib(i),ry,iy);
        for j in rlnk'range loop         -- subtract A(k-1)*x(k) from b(k)
          rlnk(j) := rlnk(j) - ry(j);
          ilnk(j) := ilnk(j) - iy(j);
        end loop;
      end loop;
      Standard_Inlined_Linear_Solvers.lusolve(rc,ic,dim,ipvt,rb(k),ib(k));
    end loop;
  end Inlined_Solve_Tail_by_lusolve;

  procedure Inlined_Solve_by_lufac
              ( dim : in integer32; 
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    Inlined_Solve_by_lufac(rb'last,dim,rc,ic,rv,iv,rb,ib,ipvt,info,ry,iy);
  end Inlined_Solve_by_lufac;

  procedure Inlined_Solve_by_lufac
              ( deg,dim : in integer32; 
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector ) is

    rlnk,ilnk : Standard_Floating_Vectors.Link_to_Vector;

  begin
    Standard_Inlined_Linear_Solvers.lufac(rc,ic,dim,ipvt,info);
    if info = 0 then
      Standard_Inlined_Linear_Solvers.lusolve(rc,ic,dim,ipvt,rb(0),ib(0));
      for k in 1..deg loop                           -- loop to compute x(k)
        Row_Matrix_Multiply(rv(k),iv(k),rb(0),ib(0),ry,iy); -- y = A(k)*x(0)
        rlnk := rb(k); ilnk := ib(k);
        for j in rlnk'range loop          -- compute b(k) = b(k) - A(k)*x(0)
          rlnk(j) := rlnk(j) - ry(j);
          ilnk(j) := ilnk(j) - iy(j);
        end loop;
        for i in 1..(k-1) loop
          Row_Matrix_Multiply(rv(k-i),iv(k-i),rb(i),ib(i),ry,iy);
          for j in rlnk'range loop         -- subtract A(k-1)*x(k) from b(k)
            rlnk(j) := rlnk(j) - ry(j);
            ilnk(j) := ilnk(j) - iy(j);
          end loop;
        end loop;
        Standard_Inlined_Linear_Solvers.lusolve(rc,ic,dim,ipvt,rb(k),ib(k));
      end loop;
    end if;
  end Inlined_Solve_by_lufac;

  procedure Inlined_Solve_by_lufco
              ( dim : in integer32; 
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_float;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    Inlined_Solve_by_lufco(rb'last,dim,rc,ic,rv,iv,rb,ib,ipvt,rcond,ry,iy);
  end Inlined_Solve_by_lufco;

  procedure Inlined_Solve_by_lufco
              ( deg,dim : in integer32; 
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_float;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector ) is

    rlnk,ilnk : Standard_Floating_Vectors.Link_to_Vector;

  begin
    Standard_Inlined_Linear_Solvers.lufco(rc,ic,dim,ipvt,ry,iy,rcond);
    if 1.0 + rcond /= 1.0 then
      Standard_Inlined_Linear_Solvers.lusolve(rc,ic,dim,ipvt,rb(0),ib(0));
      for k in 1..deg loop                           -- loop to compute x(k)
        Row_Matrix_Multiply(rv(k),iv(k),rb(0),ib(0),ry,iy); -- y = A(k)*x(0)
        rlnk := rb(k); ilnk := ib(k);
        for j in rlnk'range loop          -- compute b(k) = b(k) - A(k)*x(0)
          rlnk(j) := rlnk(j) - ry(j);
          ilnk(j) := ilnk(j) - iy(j);
        end loop;
        for i in 1..(k-1) loop
          Row_Matrix_Multiply(rv(k-i),iv(k-i),rb(i),ib(i),ry,iy);
          for j in rlnk'range loop         -- subtract A(k-1)*x(k) from b(k)
            rlnk(j) := rlnk(j) - ry(j);
            ilnk(j) := ilnk(j) - iy(j);
          end loop;
        end loop;
        Standard_Inlined_Linear_Solvers.lusolve(rc,ic,dim,ipvt,rb(k),ib(k));
      end loop;
    end if;
  end Inlined_Solve_by_lufco;

  procedure Inlined_Solve_by_lufac
              ( A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 ) is

    a0lu : constant Standard_Complex_Matrices.Link_to_Matrix := A(0);
    dim : constant integer32 := a0lu'last(1);
    deg : constant integer32 := A'last;
    rcols,icols,rb,ib : Standard_Floating_VecVecs.Link_to_VecVec;
    ry,iy : Standard_Floating_Vectors.Link_to_Vector;
    rv,iv : Standard_Floating_VecVecVecs.Link_to_VecVecVec;

  begin
    Standard_Floating_VecVecVecs.Allocate(rv,1,deg,1,dim,1,dim);
    Standard_Floating_VecVecVecs.Allocate(iv,1,deg,1,dim,1,dim);
    Standard_Matrix_Splitters.Split_Rows(A,rv,iv);
    rcols := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
    icols := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
    rb := Standard_Vector_Splitters.Allocate(deg,dim,0,1);
    ib := Standard_Vector_Splitters.Allocate(deg,dim,0,1);
    ry := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);
    iy := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);
    Standard_Matrix_Splitters.Complex_Parts(a0lu.all,rcols,icols);
    Standard_Vector_Splitters.Complex_Parts(b,rb,ib);
    Inlined_Solve_by_lufac(dim,rcols,icols,rv,iv,rb,ib,ipvt,info,ry,iy);
    Standard_Vector_Splitters.Complex_Merge(rb,ib,b);
    Standard_Floating_VecVecs.Deep_Clear(rcols);
    Standard_Floating_VecVecs.Deep_Clear(icols);
    Standard_Floating_VecVecs.Deep_Clear(rb);
    Standard_Floating_VecVecs.Deep_Clear(ib);
    Standard_Floating_Vectors.Clear(ry);
    Standard_Floating_Vectors.Clear(iy);
    Standard_Floating_VecVecVecs.Clear(rv);
    Standard_Floating_VecVecVecs.Clear(iv);
  end Inlined_Solve_by_lufac;

  procedure Inlined_Solve_by_lufco
              ( A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_float ) is

    a0lu : constant Standard_Complex_Matrices.Link_to_Matrix := A(0);
    dim : constant integer32 := a0lu'last(1);
    deg : constant integer32 := A'last;
    rcols,icols,rb,ib : Standard_Floating_VecVecs.Link_to_VecVec;
    ry,iy : Standard_Floating_Vectors.Link_to_Vector;
    rv,iv : Standard_Floating_VecVecVecs.Link_to_VecVecVec;

  begin
    Standard_Floating_VecVecVecs.Allocate(rv,1,deg,1,dim,1,dim);
    Standard_Floating_VecVecVecs.Allocate(iv,1,deg,1,dim,1,dim);
    Standard_Matrix_Splitters.Split_Rows(A,rv,iv);
    rcols := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
    icols := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
    rb := Standard_Vector_Splitters.Allocate(deg,dim,0,1);
    ib := Standard_Vector_Splitters.Allocate(deg,dim,0,1);
    ry := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);
    iy := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);
    Standard_Matrix_Splitters.Complex_Parts(a0lu.all,rcols,icols);
    Standard_Vector_Splitters.Complex_Parts(b,rb,ib);
    Inlined_Solve_by_lufco(dim,rcols,icols,rv,iv,rb,ib,ipvt,rcond,ry,iy);
    Standard_Vector_Splitters.Complex_Merge(rb,ib,b);
    Standard_Floating_VecVecs.Deep_Clear(rcols);
    Standard_Floating_VecVecs.Deep_Clear(icols);
    Standard_Floating_VecVecs.Deep_Clear(rb);
    Standard_Floating_VecVecs.Deep_Clear(ib);
    Standard_Floating_Vectors.Clear(ry);
    Standard_Floating_Vectors.Clear(iy);
    Standard_Floating_VecVecVecs.Clear(rv);
    Standard_Floating_VecVecVecs.Clear(iv);
  end Inlined_Solve_by_lufco;

end Standard_Inlined_Linearization;
