with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Complex_Vector_Norms;
with Double_Laurent_Series;
with Double_Linear_Laurent_Solvers;     use Double_Linear_Laurent_Solvers;
with Test_Double_Lseries_Matrices;

package body Double_Lseries_Newton_Steps is

  procedure Make_Series
              ( sol : in Standard_Complex_Vectors.Vector;
                deg : in integer32;
                lead : out Standard_Integer_Vectors.Vector;
                cffs : out Standard_Complex_VecVecs.Link_to_VecVec ) is

    xcffs : Standard_Complex_VecVecs.VecVec(sol'range);

  begin
    lead := (sol'range => 0);
    for i in sol'range loop
      declare
        cff : Standard_Complex_Vectors.Vector(0..deg);
      begin
        cff(0) := sol(i);
        for k in 1..deg loop
          cff(k) := Standard_Complex_Numbers.Create(0.0);
        end loop;
        xcffs(i) := new Standard_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    cffs := new Standard_Complex_VecVecs.VecVec'(xcffs);
  end Make_Series;

  procedure Make_Series
              ( sol : in Laur_Sys; deg : in integer32;
                lead : out Standard_Integer_Vectors.Vector;
                cffs : out Standard_Complex_VecVecs.Link_to_VecVec ) is

    xcffs : Standard_Complex_VecVecs.VecVec(sol'range);
    dg : constant Standard_Integer_Vectors.Vector(1..1) := (1..1 => 0);
    ldg : Standard_Complex_Laurentials.Degrees
        := new Standard_Integer_Vectors.Vector'(dg);

  begin
    lead := (sol'range => 0);
    for i in sol'range loop
      put("Number of variables : ");
      put(integer32(Number_of_Unknowns(sol(i))),1); new_line;
      declare
        cff : Standard_Complex_Vectors.Vector(0..deg);
      begin
        for k in 0..deg loop
          ldg(1) := k;
          cff(k) := Standard_Complex_Laurentials.Coeff(sol(i),ldg);
        end loop;
        Double_Laurent_Series.Normalize(deg,lead(i),cff);
        xcffs(i) := new Standard_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    cffs := new Standard_Complex_VecVecs.VecVec'(xcffs);
    Standard_Complex_Laurentials.Clear(ldg);
  end Make_Series;

  procedure Set_Leading_Exponents 
              ( lead : in out Standard_Integer_Vectors.Vector ) is

    ans : character;

  begin
    put("Leading exponents : "); put(lead); new_line;
    for i in lead'range loop
      new_line;
      loop
        put("-> the leading degree of series "); put(i,1);
        put(" : "); put(lead(i)); new_line;
        put("   Change the leading degree ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
        put("Give the new leading degree : "); get(lead(i));
      end loop;
    end loop;
  end Set_Leading_Exponents;

  function Max_Norm ( v : Standard_Complex_VecVecs.Link_to_VecVec )
                    return double_float is

    res : double_float
        := Standard_Complex_Vector_Norms.Max_Norm(v(v'first).all);
    val : double_float;

  begin
    for k in v'first+1..v'last loop
      val := Standard_Complex_Vector_Norms.Max_Norm(v(k).all);
      if val > res
       then res := val;
      end if;
    end loop;
    return res;
  end Max_Norm;

  procedure Newton_Step
              ( deg : in integer32;
                p : in Table_Vector; jp : in Table_Vector_Array;
                xlead : in out Standard_Integer_Vectors.Vector;
                xcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                ylead : in out Standard_Integer_Vectors.Vector;
                ycffs : in out Standard_Complex_VecVecs.Link_to_VecVec;
                Alead : in out Standard_Integer_Matrices.Matrix;
                Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                Blead : in out Standard_Integer_Matrices.Matrix;
                Bcffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                dxlead : in out Standard_Integer_Vectors.Vector;
                dxcffs : in out Standard_Complex_VecVecs.Link_to_VecVec;
                rlead : in out Standard_Integer_Vectors.Vector;
                rcffs : in out Standard_Complex_VecVecs.Link_to_VecVec;
                dxnrm,pxnrm : out double_float;
                verbose : in boolean := true ) is

    pivots : Standard_Integer_Vectors.Vector(xlead'range);
    acff,bcff : Standard_Complex_Vectors.Link_to_Vector;

    use Standard_Complex_Numbers;

  begin
    if verbose
     then put_line("Evaluating the table vector ...");
    end if;
    Eval(deg,p,xlead,xcffs,ylead,ycffs.all,verbose);
    if verbose
     then Double_Linear_Laurent_Solvers.Write(ylead,ycffs,"y");
    end if;
    if verbose
     then put_line("Evaluating the table vector array ...");
    end if;
    Eval(deg,jp,xlead,xcffs,Alead,Acffs,verbose);
    if verbose then
      Test_Double_Lseries_Matrices.Copy
        (p.nbt,p.nbt,deg,Alead,Acffs,Blead,Bcffs);
      Double_Linear_Laurent_Solvers.Write(Blead,Bcffs,"B");
    end if;
    LU_Factorization(p.nbt,p.nbt,deg,Alead,Acffs,pivots);
    if verbose
     then put("The pivots : "); put(pivots); new_line;
    end if;
    for i in pivots'range loop
      dxlead(i) := ylead(i);
      acff := ycffs(pivots(i));
      bcff := dxcffs(i);
      for k in 0..deg loop
        bcff(k) := -acff(k);
      end loop;
    end loop;
    if verbose
     then Double_Linear_Laurent_Solvers.Write(dxlead,dxcffs,"b");
    end if;
    Forward_Substitution(deg,Alead,Acffs,dxlead,dxcffs,ylead,ycffs);
    Backward_Substitution(deg,Alead,Acffs,ylead,ycffs,dxlead,dxcffs);
    dxnrm := Max_Norm(dxcffs);
    pxnrm := Max_Norm(rcffs);
    if verbose then
      Double_Linear_Laurent_Solvers.Write(dxlead,dxcffs,"dx");
      Matrix_Vector_Product(deg,Blead,Bcffs,dxlead,dxcffs,rlead,rcffs);
      Double_Linear_Laurent_Solvers.Write(rlead,rcffs,"r");
      put("Maximum |dx| :"); put(dxnrm,3);
      put("  residual :"); put(pxnrm,3); new_line;
    end if;
    for i in dxlead'range loop
      Double_Laurent_Series.Add(deg,
        xlead(i),dxlead(i),xcffs(i).all,dxcffs(i).all,ylead(i),ycffs(i).all);
      xlead(i) := ylead(i);
      acff := ycffs(i);
      bcff := xcffs(i);
      for k in 0..deg loop
        bcff(k) := acff(k);
      end loop;
    end loop;
  end Newton_Step;

  procedure Run_Newton_Steps
              ( deg : in integer32;
                tv : in Table_Vector; tva : in Table_Vector_Array;
                xlead : in out Standard_Integer_Vectors.Vector;
                xcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                dxnrm,pxnrm : out double_float;
                numit : out integer32; maxit : in integer32 := 4;
                dxtol : in double_float := 1.0E-8;
                pxtol : in double_float := 1.0E-8;
                verbose : in boolean := true ) is
             
    neq : constant integer32 := tva'last;
    nvr : constant integer32 := xlead'last;
    ylead,dxlead,rlead : Standard_Integer_Vectors.Vector(1..nvr);
    ycffs,dxcffs,rcffs : Standard_Complex_VecVecs.Link_to_VecVec;
    Alead,Blead : Standard_Integer_Matrices.Matrix(1..neq,1..nvr);
    Acffs,Bcffs : Standard_Complex_VecVecVecs.Link_to_VecVecVec;
    stepcnt : integer32 := 1;

  begin
    Double_Linear_Laurent_Solvers.Allocate_Series_Coefficients(nvr,deg,ycffs);
    Double_Linear_Laurent_Solvers.Allocate_Series_Coefficients(nvr,deg,dxcffs);
    Double_Linear_Laurent_Solvers.Allocate_Series_Coefficients(nvr,deg,rcffs);
    Standard_Complex_VecVecVecs.Allocate(Acffs,1,neq,1,nvr,0,deg);
    Standard_Complex_VecVecVecs.Allocate(Bcffs,1,neq,1,nvr,0,deg);
    loop
      if verbose
       then put("Step "); put(stepcnt,1); put_line(" ...");
      end if;
      Double_Lseries_Newton_Steps.Newton_Step
        (deg,tv,tva,xlead,xcffs,ylead,ycffs,Alead,Acffs,Blead,Bcffs,
         dxlead,dxcffs,rlead,rcffs,dxnrm,pxnrm,verbose);
      if verbose
       then Double_Linear_Laurent_Solvers.Write(xlead,xcffs,"x");
      end if;
      exit when (stepcnt = maxit);
      exit when ((dxnrm < dxtol) and (pxnrm < pxtol));
      stepcnt := stepcnt + 1;
    end loop;
    numit := stepcnt;
    Standard_Complex_VecVecs.Deep_Clear(ycffs);
    Standard_Complex_VecVecs.Deep_Clear(dxcffs);
    Standard_Complex_VecVecs.Deep_Clear(rcffs);
    Standard_Complex_VecVecVecs.Clear(Acffs);
    Standard_Complex_VecVecVecs.Clear(Bcffs);
  end Run_Newton_Steps;

end Double_Lseries_Newton_Steps;
