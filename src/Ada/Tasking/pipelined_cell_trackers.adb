with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Double_Double_Vectors;
with DoblDobl_Complex_Vectors;
with Quad_Double_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_Linear_Solvers;
With DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Linear_Solvers;
with Standard_Integer_Matrices;
with Standard_Integer64_Matrices;
with Standard_Integer_Linear_Solvers;
with Standard_Integer64_Linear_Solvers;
with Supports_of_Polynomial_Systems;
with Standard_Radial_Solvers;
with Standard_Binomial_Solvers;
with Standard_Binomial_Systems;
with DoblDobl_Radial_Solvers;
with DoblDobl_Binomial_Solvers;
with DoblDobl_Binomial_Systems;
with QuadDobl_Radial_Solvers;
with QuadDobl_Binomial_Systems;
with QuadDobl_Binomial_Solvers;
with Polyhedral_Start_Systems;          use Polyhedral_Start_Systems;
with Single_Polyhedral_Trackers;        use Single_Polyhedral_Trackers;

package body Pipelined_Cell_Trackers is

  procedure Standard_Track_Cell
              ( sem : in out Semaphore.Lock; idtask,n,r : in integer32; 
                mix : in Standard_Integer_Vectors.Vector;
                mic : in Mixed_Cell;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                cff : in Standard_Complex_VecVecs.VecVec;
                dpw : in Standard_Floating_VecVecs.Link_to_VecVec;
                cft : in Standard_Complex_VecVecs.Link_to_VecVec;
                epv : in Exponent_Vectors.Exponent_Vectors_Array;
                hom : in Standard_Complex_Laur_SysFun.Eval_Coeff_Laur_Sys;
                ejf : in Standard_Complex_Laur_JacoMats.Eval_Coeff_Jaco_Mat;
                jmf : in Standard_Complex_Laur_JacoMats.Mult_Factors;
                q : in Standard_Complex_Laur_Systems.Laur_Sys;
                tmv : in out natural32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                last : in out Standard_Complex_Solutions.Solution_List ) is

    use Supports_of_Polynomial_Systems;

    s_c : Standard_Complex_Vectors.Vector(1..2*n);    -- for fully mixed
    CC : Standard_Complex_Matrices.Matrix(1..n,1..n); -- for semi mixed
    piv : Standard_Integer_Vectors.Vector(1..n);
    pdetU : natural32;
    info : integer32;
    A,M,U : Standard_Integer_Matrices.Matrix(q'range,q'range);
    b,wrk : Standard_Complex_Vectors.Vector(q'range);
    brd,logbrd,logx,e10x : Standard_Floating_Vectors.Vector(b'range);
    bsc : Standard_Complex_Vectors.Vector(b'range);
    ls : Standard_Complex_Solutions.Link_to_Solution;
    ptr : Standard_Complex_Solutions.Solution_List;

  begin
    if r = n then
      Select_Coefficients(cff,epv,mic.pts.all,s_c);
      Fully_Mixed_To_Binomial_Format(s_c,mic.pts.all,A,b);
    else
      Select_Subsystem_to_Matrix_Format(cff,epv,mix,mic.pts.all,A,CC,b);
      Standard_Complex_Linear_Solvers.lufac(CC,n,piv,info);
      Standard_Complex_Linear_Solvers.lusolve(CC,n,piv,b);
    end if;
    U := A;
    Standard_Integer_Linear_Solvers.Upper_Triangulate(M,U);
    pdetU := Volume_of_Diagonal(U);
    Semaphore.Request(sem);
    ptr := last;
    Allocate(sols,last,n,integer32(pdetU));
    if Standard_Complex_Solutions.Is_Null(ptr)
     then ptr := sols;
     else ptr := Standard_Complex_Solutions.Tail_Of(ptr);
    end if;
    brd := Standard_Radial_Solvers.Radii(b);
    bsc := Standard_Radial_Solvers.Scale(b,brd);
    Standard_Binomial_Solvers.Solve_Upper_Square(U,bsc,ptr);
    Semaphore.Release(sem);
    logbrd := Standard_Radial_Solvers.Log10(brd);
    logx := Standard_Radial_Solvers.Radial_Upper_Solve(U,logbrd);
    logx := Standard_Radial_Solvers.Multiply(M,logx);
    e10x := Standard_Radial_Solvers.Exp10(logx);
    for i in 1..pdetU loop
      ls := Standard_Complex_Solutions.Head_Of(ptr);
      Standard_Binomial_Systems.Eval(M,ls.v,wrk);
      ls.v := wrk;
      Standard_Radial_Solvers.Multiply(ls.v,e10x);
      Track_Path(mix,lif,mic.nor,cff,dpw,cft,epv,hom,ejf,jmf,ls);
      ptr := Standard_Complex_Solutions.Tail_Of(ptr);
    end loop;
    tmv := tmv + pdetU;
  end Standard_Track_Cell;

  procedure DoblDobl_Track_Cell
              ( sem : in out Semaphore.Lock; idtask,n,r : in integer32; 
                mix : in Standard_Integer_Vectors.Vector;
                mic : in Mixed_Cell;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                cff : in DoblDobl_Complex_VecVecs.VecVec;
                dpw : in Standard_Floating_VecVecs.Link_to_VecVec;
                cft : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                epv : in Exponent_Vectors.Exponent_Vectors_Array;
                hom : in DoblDobl_Complex_Laur_SysFun.Eval_Coeff_Laur_Sys;
                ejf : in DoblDobl_Complex_Laur_JacoMats.Eval_Coeff_Jaco_Mat;
                jmf : in DoblDobl_Complex_Laur_JacoMats.Mult_Factors;
                q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                tmv : in out natural32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                last : in out DoblDobl_Complex_Solutions.Solution_List ) is

    use Supports_of_Polynomial_Systems;

    s_c : DoblDobl_Complex_Vectors.Vector(1..2*n);    -- for fully mixed
    CC : DoblDobl_Complex_Matrices.Matrix(1..n,1..n); -- for semi mixed
    piv : Standard_Integer_Vectors.Vector(1..n);
    pdetU : natural64;
    info : integer32;
    A,M,U : Standard_Integer64_Matrices.Matrix(q'range,q'range);
    b,wrk : DoblDobl_Complex_Vectors.Vector(q'range);
    brd,logbrd,logx,e10x : Double_Double_Vectors.Vector(b'range);
    bsc : DoblDobl_Complex_Vectors.Vector(b'range);
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    ptr : DoblDobl_Complex_Solutions.Solution_List;

  begin
    if r = n then
      Select_Coefficients(cff,epv,mic.pts.all,s_c);
      Fully_Mixed_To_Binomial_Format(s_c,mic.pts.all,A,b);
    else
      Select_Subsystem_to_Matrix_Format(cff,epv,mix,mic.pts.all,A,CC,b);
      DoblDobl_Complex_Linear_Solvers.lufac(CC,n,piv,info);
      DoblDobl_Complex_Linear_Solvers.lusolve(CC,n,piv,b);
    end if;
    U := A;
    Standard_Integer64_Linear_Solvers.Upper_Triangulate(M,U);
    pdetU := Volume_of_Diagonal(U);
    Semaphore.Request(sem);
    ptr := last;
    Allocate(sols,last,n,integer32(pdetU));
    if DoblDobl_Complex_Solutions.Is_Null(ptr)
     then ptr := sols;
     else ptr := DoblDobl_Complex_Solutions.Tail_Of(ptr);
    end if;
    brd := DoblDobl_Radial_Solvers.Radii(b);
    bsc := DoblDobl_Radial_Solvers.Scale(b,brd);
    DoblDobl_Binomial_Solvers.Solve_Upper_Square(U,bsc,ptr);
    Semaphore.Release(sem);
    logbrd := DoblDobl_Radial_Solvers.Log10(brd);
    logx := DoblDobl_Radial_Solvers.Radial_Upper_Solve(U,logbrd);
    logx := DoblDobl_Radial_Solvers.Multiply(M,logx);
    e10x := DoblDobl_Radial_Solvers.Exp10(logx);
    for i in 1..pdetU loop
      ls := DoblDobl_Complex_Solutions.Head_Of(ptr);
      DoblDobl_Binomial_Systems.Eval(M,ls.v,wrk);
      ls.v := wrk;
      DoblDobl_Radial_Solvers.Multiply(ls.v,e10x);
      Track_Path(mix,lif,mic.nor,cff,dpw,cft,epv,hom,ejf,jmf,ls);
      ptr := DoblDobl_Complex_Solutions.Tail_Of(ptr);
    end loop;
    tmv := tmv + natural32(pdetU);
  end DoblDobl_Track_Cell;

  procedure QuadDobl_Track_Cell
              ( sem : in out Semaphore.Lock; idtask,n,r : in integer32; 
                mix : in Standard_Integer_Vectors.Vector;
                mic : in Mixed_Cell;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                cff : in QuadDobl_Complex_VecVecs.VecVec;
                dpw : in Standard_Floating_VecVecs.Link_to_VecVec;
                cft : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                epv : in Exponent_Vectors.Exponent_Vectors_Array;
                hom : in QuadDobl_Complex_Laur_SysFun.Eval_Coeff_Laur_Sys;
                ejf : in QuadDobl_Complex_Laur_JacoMats.Eval_Coeff_Jaco_Mat;
                jmf : in QuadDobl_Complex_Laur_JacoMats.Mult_Factors;
                q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                tmv : in out natural32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                last : in out QuadDobl_Complex_Solutions.Solution_List ) is

    use Supports_of_Polynomial_Systems;

    s_c : QuadDobl_Complex_Vectors.Vector(1..2*n);    -- for fully mixed
    CC : QuadDobl_Complex_Matrices.Matrix(1..n,1..n); -- for semi mixed
    piv : Standard_Integer_Vectors.Vector(1..n);
    pdetU : natural64;
    info : integer32;
    A,M,U : Standard_Integer64_Matrices.Matrix(q'range,q'range);
    b,wrk : QuadDobl_Complex_Vectors.Vector(q'range);
    brd,logbrd,logx,e10x : Quad_Double_Vectors.Vector(b'range);
    bsc : QuadDobl_Complex_Vectors.Vector(b'range);
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    ptr : QuadDobl_Complex_Solutions.Solution_List;

  begin
    if r = n then
      Select_Coefficients(cff,epv,mic.pts.all,s_c);
      Fully_Mixed_To_Binomial_Format(s_c,mic.pts.all,A,b);
    else
      Select_Subsystem_to_Matrix_Format(cff,epv,mix,mic.pts.all,A,CC,b);
      QuadDobl_Complex_Linear_Solvers.lufac(CC,n,piv,info);
      QuadDobl_Complex_Linear_Solvers.lusolve(CC,n,piv,b);
    end if;
    U := A;
    Standard_Integer64_Linear_Solvers.Upper_Triangulate(M,U);
    pdetU := Volume_of_Diagonal(U);
    Semaphore.Request(sem);
    ptr := last;
    Allocate(sols,last,n,integer32(pdetU));
    if QuadDobl_Complex_Solutions.Is_Null(ptr)
     then ptr := sols;
     else ptr := QuadDobl_Complex_Solutions.Tail_Of(ptr);
    end if;
    brd := QuadDobl_Radial_Solvers.Radii(b);
    bsc := QuadDobl_Radial_Solvers.Scale(b,brd);
    QuadDobl_Binomial_Solvers.Solve_Upper_Square(U,bsc,ptr);
    Semaphore.Release(sem);
    logbrd := QuadDobl_Radial_Solvers.Log10(brd);
    logx := QuadDobl_Radial_Solvers.Radial_Upper_Solve(U,logbrd);
    logx := QuadDobl_Radial_Solvers.Multiply(M,logx);
    e10x := QuadDobl_Radial_Solvers.Exp10(logx);
    for i in 1..pdetU loop
      ls := QuadDobl_Complex_Solutions.Head_Of(ptr);
      QuadDobl_Binomial_Systems.Eval(M,ls.v,wrk);
      ls.v := wrk;
      QuadDobl_Radial_Solvers.Multiply(ls.v,e10x);
      Track_Path(mix,lif,mic.nor,cff,dpw,cft,epv,hom,ejf,jmf,ls);
      ptr := QuadDobl_Complex_Solutions.Tail_Of(ptr);
    end loop;
    tmv := tmv + natural32(pdetU);
  end QuadDobl_Track_Cell;

end Pipelined_Cell_Trackers;
