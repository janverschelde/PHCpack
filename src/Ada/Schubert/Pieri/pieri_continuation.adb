with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Complex_VecVecs;           use Standard_Complex_VecVecs;
with Standard_Natural_Matrices_io;       use Standard_Natural_Matrices_io;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Homotopy;
with Standard_IncFix_Continuation;       use Standard_IncFix_Continuation;
with Standard_Root_Refiners;             use Standard_Root_Refiners;
with Determinantal_Systems;              use Determinantal_Systems;
with Plane_Representations;              use Plane_Representations;
with Curves_into_Grassmannian;           use Curves_into_Grassmannian;
with Curves_into_Grassmannian_io;        use Curves_into_Grassmannian_io;
with Verification_with_Determinants;     use Verification_with_Determinants;

package body Pieri_Continuation is

-- AUXILIARY FORMAT CONVERTORS :

  function Create ( v : Standard_Complex_Vectors.Vector;
                    conpar : integer32 ) return Solution is

  -- DESCRIPTION :
  --   Returns a solution that contains a solution with vector v in it.
  --   The v(conpar) is the continuation parameter.

    sol : Solution(v'last-1);

  begin
    sol.t := v(conpar);
    sol.m := 1;
    sol.v(v'first..conpar-1) := v(v'first..conpar-1);
    sol.v(conpar..sol.v'last) := v(conpar+1..v'last);
    sol.err := 0.0;
    sol.rco := 0.0;
    sol.res := 0.0;
    return sol;
  end Create;

  function Create ( v : Standard_Complex_Vectors.Vector )
                  return Solution_List is

  -- DESCRIPTION :
  --   Returns a solution list that contains a solution with vector v in it.

  -- NOTE : only needed in the original Pieri homotopy algorithm.

    res : Solution_List;

  begin
    Add(res,Create(v,v'last));
    return res;
  end Create;

  function Create ( v : Standard_Complex_VecVecs.VecVec;
                    conpar : integer32 ) return Solution_List is

  -- DESCRIPTION :
  --   Returns a solution list that contains the vectors in v.

  -- REQUIRED : v is not empty.
  --   The value for the continuation parameter is in v(conpar).

    res,res_last : Solution_List;

  begin
    for i in v'range loop
      Append(res,res_last,Create(v(i).all,conpar));
    end loop;
    return res;
  end Create;

  function Convert ( sols : Solution_List )
                   return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Converts the solutions in sols to the vector representation.

  -- NOTE : only needed in the original Pieri homotopy algorithm.

    ls : constant Link_to_Solution := Head_Of(sols);
    res : Standard_Complex_Vectors.Vector(1..ls.n+1);

  begin
    res(ls.v'range) := ls.v;
    res(res'last) := ls.t;
    return res;
  end Convert;

  function Convert ( sol : Solution; conpar : integer32 )
                   return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Converts the solution to the vector representation, where the
  --   continuation parameter is put into the entry indexed by conpar.

    res : Standard_Complex_Vectors.Vector(1..sol.n+1);

  begin
    res(sol.v'first..conpar-1) := sol.v(sol.v'first..conpar-1);
    res(conpar) := sol.t;
    res(conpar+1..sol.n+1) := sol.v(conpar..sol.v'last);
    return res;
  end Convert;

  function Convert ( sols : Solution_List; conpar : integer32 )
                   return Standard_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Converts the solutions in sols to the vector representation.

    res : Standard_Complex_VecVecs.VecVec(1..integer32(Length_Of(sols)));
    tmp : Solution_List := sols;

  begin
    for i in res'range loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
        v : Standard_Complex_Vectors.Vector(ls.v'first..ls.v'last+1);
      begin
        v := Convert(ls.all,conpar);
        res(i) := new Standard_Complex_Vectors.Vector'(v);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Convert;

  function Square ( n : integer32; p : Poly_Sys ) return Poly_Sys is

  -- DESCRIPTION :
  --   Returns a n-by-n system, by adding random linear combinations of 
  --   the polynomials p(i), i > n, to the first n polynomials.

    res : Poly_Sys(1..n);
    acc : Poly;

  begin
    for i in res'range loop
      Copy(p(i),res(i));
      for j in n+1..p'last loop
        acc := Random1*p(j);
        Add(res(i),acc);
        Clear(acc);
      end loop;
    end loop;
    return res;
  end Square;

  procedure Refine_Roots
                ( p : in Poly_Sys; sols : in out Solution_List ) is

    epsxa : constant double_float := 10.0**(-8);
    epsfa : constant double_float := 10.0**(-8);
    tolsing : constant double_float := 10.0**(-8);
    max : constant natural32 := 3;
    numit : natural32 := 0;
    deflate : boolean := false;

  begin
    Silent_Root_Refiner(p,sols,epsxa,epsfa,tolsing,numit,max,deflate);
  end Refine_Roots;

  procedure Refine_Roots
                ( file : in file_type; p : in Poly_Sys;
                  sols : in out Solution_List; report : in boolean ) is

    epsxa : constant double_float := 10.0**(-8);
    epsfa : constant double_float := 10.0**(-8);
    tolsing : constant double_float := 10.0**(-8);
    max : constant natural32 := 3;
    numit : natural32 := 0;
    deflate : boolean := false;

  begin
    if report
     then Reporting_Root_Refiner
            (file,p,sols,epsxa,epsfa,tolsing,numit,max,deflate,false);
     else Silent_Root_Refiner(p,sols,epsxa,epsfa,tolsing,numit,max,deflate);
    end if;
  end Refine_Roots;

  procedure Determinantal_Polynomial_Continuation
                ( file : in file_type; lochom : in Poly_Sys;  -- out added
                  locsol : in out Standard_Complex_Vectors.Vector;
                  report,outlog : in boolean ) is

  -- DESCRIPTION :
  --   Given a homotopy and solution with as last variable the continuation
  --   parameter the polynomial continuation will launched.
  --   This homotopy uses the determinantal expansions defined in the paper.

  -- NOTE : 
  --   This routine is only used in the original Pieri homotopy algorithm.

    sols : Solution_List := Create(locsol);
    thesys,squsys : Poly_Sys(locsol'first..locsol'last-1);
    square_case : boolean;
   -- ran : Complex_Number;

  begin
   -- for i in lochom'range loop           -- added to deal with real input
   --   ran := Random1;
   --   Mul(lochom(i),ran);
   -- end loop;
    if lochom'length + 1 = locsol'length then
      square_case := true;
      Standard_Homotopy.Create(lochom,locsol'last);
    else
      square_case := false;
      squsys := Square(locsol'last-1,lochom);
      if outlog then
        put_line(file,"The homotopy as square system : ");
        put_line(file,squsys);
      end if;
      Standard_Homotopy.Create(squsys,locsol'last);
    end if;
    declare
      procedure Sil_Cont is
        new Silent_Continue(Max_Norm,Standard_Homotopy.Eval,
                            Standard_Homotopy.Diff,Standard_Homotopy.Diff);
      procedure Rep_Cont is
        new Reporting_Continue(Max_Norm,Standard_Homotopy.Eval,
                               Standard_Homotopy.Diff,Standard_Homotopy.Diff);
    begin
      if report
       then Rep_Cont(file,sols,false,target=>Create(1.0));
       else Sil_Cont(sols,false,target=>Create(1.0));
      end if;
    end;
    if square_case then
      thesys := Eval(lochom,Create(1.0),locsol'last);
      Refine_Roots(file,thesys,sols,report);
    else
      thesys := Eval(squsys,Create(1.0),locsol'last);
      Refine_Roots(file,thesys,sols,report);
      Clear(squsys);
    end if;
    Clear(thesys);
    locsol := Convert(sols);
    Clear(sols);
    Standard_Homotopy.Clear;
  end Determinantal_Polynomial_Continuation;

  procedure Determinantal_Polynomial_Continuation
                ( lochom : in Poly_Sys; conpar : in integer32;
                  locsols : in out Standard_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Given a homotopy and solution with as continuation parameter
  --   the variable indexed by conpar, the continuation will be launched.
  --   This homotopy uses the determinantal expansions defined in the paper.

    sols : Solution_List := Create(locsols,conpar);
    firstsol : constant Standard_Complex_Vectors.Vector
             := locsols(locsols'first).all;
    thesys,squsys : Poly_Sys(firstsol'first..firstsol'last-1);
    square_case : boolean;

  begin
    if lochom'length + 1 = firstsol'length then
      square_case := true;
      Standard_Homotopy.Create(lochom,conpar);
    else
      square_case := false;
      squsys := Square(firstsol'last-1,lochom);
      Standard_Homotopy.Create(squsys,conpar);
    end if;
    declare
      procedure Sil_Cont is
        new Silent_Continue(Max_Norm,Standard_Homotopy.Eval,
                            Standard_Homotopy.Diff,Standard_Homotopy.Diff);
    begin
      Sil_Cont(sols,false,target=>Create(1.0));
    end;
    if square_case then
      thesys := Eval(lochom,Create(1.0),conpar);
      Refine_Roots(thesys,sols);
    else
      thesys := Eval(squsys,Create(1.0),conpar);
      Refine_Roots(thesys,sols);
      Clear(squsys);
    end if;
    Clear(thesys);
    locsols := Convert(sols,conpar);
    Clear(sols);
    Standard_Homotopy.Clear;
  end Determinantal_Polynomial_Continuation;

  procedure Determinantal_Polynomial_Continuation
                ( file : in file_type; lochom : in Poly_Sys;
                  conpar : in integer32;
                  locsols : in out Standard_Complex_VecVecs.VecVec;
                  report,outlog : in boolean ) is

  -- DESCRIPTION :
  --   Given a homotopy and solution with as continuation parameter
  --   the variable indexed by conpar, the continuation will be launched.
  --   This homotopy uses the determinantal expansions defined in the paper.

    sols : Solution_List := Create(locsols,conpar);
    firstsol : constant Standard_Complex_Vectors.Vector
             := locsols(locsols'first).all;
    thesys,squsys : Poly_Sys(firstsol'first..firstsol'last-1);
    square_case : boolean;

  begin
    if lochom'length + 1 = firstsol'length then
      square_case := true;
      Standard_Homotopy.Create(lochom,conpar);
    else
      square_case := false;
      squsys := Square(firstsol'last-1,lochom);
      if outlog then
        put_line(file,"The homotopy as square system : ");
        put_line(file,squsys);
      end if;
      Standard_Homotopy.Create(squsys,conpar);
    end if;
    declare
      procedure Sil_Cont is
        new Silent_Continue(Max_Norm,Standard_Homotopy.Eval,
                            Standard_Homotopy.Diff,Standard_Homotopy.Diff);
      procedure Rep_Cont is
        new Reporting_Continue(Max_Norm,Standard_Homotopy.Eval,
                               Standard_Homotopy.Diff,Standard_Homotopy.Diff);
    begin
      if report
       then Rep_Cont(file,sols,false,target=>Create(1.0));
       else Sil_Cont(sols,false,target=>Create(1.0));
      end if;
    end;
    if square_case then
      thesys := Eval(lochom,Create(1.0),conpar);
      Refine_Roots(file,thesys,sols,report);
    else
      thesys := Eval(squsys,Create(1.0),conpar);
      Refine_Roots(file,thesys,sols,report);
      Clear(squsys);
    end if;
    Clear(thesys);
    locsols := Convert(sols,conpar);
    Clear(sols);
    Standard_Homotopy.Clear;
  end Determinantal_Polynomial_Continuation;

-- TARGET ROUTINES :

  procedure Trace_Paths ( file : in file_type; homsys : in Poly_Sys;
                          locmap : in Standard_Natural_Matrices.Matrix;
                          report,outlog : in boolean;
                          plane : in out Standard_Complex_Matrices.Matrix ) is

  -- NOTE :
  --   This routine is only used in the original Pieri homotopy algorithm.

    plavec : constant Standard_Complex_Vectors.Vector := Vector_Rep(plane);
    solvec : Standard_Complex_Vectors.Vector(1..plavec'last+1);
    lochom : Poly_Sys(homsys'range) := Localize(locmap,homsys);
    plaloc : constant Standard_Complex_Matrices.Matrix
           := Localize(locmap,plane);
    solloc : constant Standard_Complex_Vectors.Vector
           := Vector_Rep(locmap,plaloc);
    locsol : Standard_Complex_Vectors.Vector(1..solloc'last+1);
    evahom : Standard_Complex_Vectors.Vector(homsys'range);
    evaloc : Standard_Complex_Vectors.Vector(lochom'range);

  begin
    if outlog
     then put_line(file,"The localization pattern :"); put(file,locmap);
          put_line(file,"The localized homotopy : "); put_line(file,lochom);
    end if;
   -- put_line(file,"The plane in local coordinates : "); put(file,plaloc,2);
    solvec(plavec'range) := plavec;
    solvec(solvec'last) := Create(0.0);
   -- put_line(file,"The solution vector : "); put_line(file,solvec);
    evahom := Eval(homsys,solvec);
   -- put_line(file,"Evaluation of solution vector at homotopy for t=0 :");
   -- put_line(file,evahom);   
    put(file,"Residual of start solution :"); 
    put(file,Max_Norm(evahom),3); put_line(file,".");
 -- put_line(file,"The localized vector representation of the plane at t=0 :");
   -- put_line(file,solloc);
    locsol(solloc'range) := solloc;
    locsol(locsol'last) := Create(0.0);
    evaloc := Eval(lochom,locsol);
    put(file,"Residual of localized plane at t=0 :");
    put(file,Max_Norm(evaloc),3); put_line(file,".");
   -- Linear_Polynomial_Continuation(file,lochom,locsol,report,outlog);
    Determinantal_Polynomial_Continuation(file,lochom,locsol,report,outlog);
 -- put_line(file,"The localized vector representation of the plane at t=1 :");
   -- put_line(file,locsol);
    evaloc := Eval(lochom,locsol);
   -- put_line(file,"Evaluation of localized vector at homotopy for t=1 :");
   -- put_line(file,evaloc);
    put(file,"Residual of localized plane at t=1 :");
    put(file,Max_Norm(evaloc),3); put_line(file,".");
    plane := Matrix_Rep(locmap,locsol(1..locsol'last-1));
   -- put_line(file,"The solution plane at t=1 :");
   -- put(file,plane,2);
    solvec(plavec'range) := Vector_Rep(plane);
    solvec(solvec'last) := Create(1.0);
    evahom := Eval(homsys,solvec);
   -- put_line(file,"Evaluation of solution vector at homotopy for t=1 :");
   -- put_line(file,evahom);
    put(file,"Residual of target solution :");
    put(file,Max_Norm(evahom),3); put_line(file,".");
    Clear(lochom);
  end Trace_Paths;

  procedure Trace_Paths
              ( file : in file_type; homsys : in Poly_Sys;
                locmap : in Standard_Natural_Matrices.Matrix;
                report,outlog : in boolean;
                planes : in Standard_Complex_VecMats.VecMat ) is

    lochom : Poly_Sys(homsys'range) := Localize(locmap,homsys);
    locsols : Standard_Complex_VecVecs.VecVec(planes'range);
    evaloc : Standard_Complex_Vectors.Vector(lochom'range);
    evahom : Standard_Complex_Vectors.Vector(homsys'range);
    lastvar : integer32;

  begin
    if outlog then
      put_line(file,"The localized homotopy :");
      put_line(file,lochom);
    end if;
    for i in planes'range loop               -- solution planes into vectors
      declare
        plaloc : constant Standard_Complex_Matrices.Matrix
               := Localize(locmap,planes(i).all);
        solloc : constant Standard_Complex_Vectors.Vector
               := Vector_Rep(locmap,plaloc);
        locsol : Standard_Complex_Vectors.Vector(1..solloc'last+1);
      begin
        locsol(solloc'range) := solloc;
        locsol(locsol'last) := Create(0.0);
        locsols(i) := new Standard_Complex_Vectors.Vector'(locsol);
        if outlog then
          evaloc := Eval(lochom,locsol);
          put(file,"Residual of localized plane at start :");
          put(file,Max_Norm(evaloc),3); put_line(file,".");
        end if;
      end;
    end loop;
    lastvar := locsols(locsols'first)'last;
    Determinantal_Polynomial_Continuation
      (file,lochom,lastvar,locsols,report,outlog);
    for i in locsols'range loop             -- solution vectors into planes
      planes(i).all := Matrix_Rep(locmap,locsols(i)(1..locsols(i)'last-1));
      if outlog then
        evaloc := Eval(lochom,locsols(i).all);
        put(file,"Residual of localized plane at target :");
        put(file,Max_Norm(evaloc),3); put_line(file,".");
        declare
          plavec : constant Standard_Complex_Vectors.Vector
                 := Vector_Rep(planes(i).all);
          solvec : Standard_Complex_Vectors.Vector
                     (plavec'first..plavec'last+1);
        begin
          solvec(plavec'range) := plavec;
          solvec(solvec'last) := Create(1.0);
          evahom := Eval(homsys,solvec);
         put(file,"Residual of target solution :");
          put(file,Max_Norm(evahom),3); put_line(file,".");
        end;
      end if;
    end loop;
    Clear(lochom); Clear(locsols);
  end Trace_Paths;

  procedure Quantum_Trace_Paths
              ( m,p,q : in natural32; nd : in Node;
              -- (m,p,q) are needed for the symbol table
                homsys : in Poly_Sys; conpar,s_mode : in natural32;
                locmap : in Standard_Natural_Matrices.Matrix;
                planes : in Standard_Complex_VecMats.VecMat ) is

    lochom : Poly_Sys(homsys'range)
           := Column_Localize(nd.top,nd.bottom,locmap,homsys);
    locsols : Standard_Complex_VecVecs.VecVec(planes'range);
    addmix : integer32;

  begin
    if nd.tp = mixed
     then addmix := 1;
     else addmix := 0;
    end if;
    for i in planes'range loop               -- solution planes into vectors
      declare
        plaloc : constant Standard_Complex_Matrices.Matrix
               := Localize(locmap,planes(i).all);
        solloc : constant Standard_Complex_Vectors.Vector
               := Column_Vector_Rep(locmap,plaloc);
        locsol : Standard_Complex_Vectors.Vector(1..solloc'last+2+addmix);
      begin
        locsol(solloc'range) := solloc;
        if s_mode = 0
         then locsol(locsol'last-1) := Create(0.0); -- start value for s is 0
         else locsol(locsol'last-1) := Create(1.0); -- start value for s is 1
        end if;
        locsol(locsol'last)   := Create(0.0);   -- start value for t is 0
        if addmix = 1
         then locsol(locsol'last-2) := Create(1.0);   -- additional s-value
        end if;
        locsols(i) := new Standard_Complex_Vectors.Vector'(locsol);
      end;
    end loop;
    Determinantal_Polynomial_Continuation(lochom,integer32(conpar),locsols);
    for i in locsols'range loop             -- solution vectors into planes
      planes(i).all
        := Column_Matrix_Rep(locmap,locsols(i)(1..locsols(i)'last-2-addmix));
      declare
        plavec : constant Standard_Complex_Vectors.Vector
               := Column_Vector_Rep(nd.top,nd.bottom,planes(i).all);
        solvec : Standard_Complex_Vectors.Vector
                   (plavec'first..plavec'last+2+addmix);
      begin
        solvec(plavec'range) := plavec;
        solvec(solvec'last)   := Create(1.0);   -- target value for t
        solvec(solvec'last-1) := locsols(i)(locsols(i)'last-1);
        if addmix = 1
         then solvec(solvec'last-2)   := locsols(i)(locsols(i)'last-2);
        end if;
      end;
    end loop;
    Clear(lochom); Clear(locsols);
  end Quantum_Trace_Paths;

  procedure Quantum_Trace_Paths
              ( file : in file_type; m,p,q : in natural32; nd : in Node;
              -- (m,p,q) are needed for the symbol table
                xpm : in Standard_Complex_Poly_Matrices.Matrix;
                s : in Standard_Complex_Vectors.Vector; ip : in VecMat;
                homsys : in Poly_Sys; conpar,s_mode : in natural32;
                locmap : in Standard_Natural_Matrices.Matrix;
                report,outlog : in boolean;
                planes : in Standard_Complex_VecMats.VecMat ) is

    lochom : Poly_Sys(homsys'range)
           := Column_Localize(nd.top,nd.bottom,locmap,homsys);
    locsols : Standard_Complex_VecVecs.VecVec(planes'range);
    evaloc : Standard_Complex_Vectors.Vector(lochom'range);
    evahom : Standard_Complex_Vectors.Vector(homsys'range);
    addmix : integer32;

  begin
    if outlog then
      put_line(file,"The localization map : ");   put(file,locmap);
      if nd.tp = mixed
       then Two_Set_up_Symbol_Table(m,p,q,nd.top,nd.bottom);
       else One_Set_up_Symbol_Table(m,p,q,nd.top,nd.bottom);
      end if;
      Reduce_Symbols(nd.top,nd.bottom,locmap);
      put_line(file,"The localized homotopy : "); put_line(file,lochom);
    end if;
    if nd.tp = mixed
     then addmix := 1;
     else addmix := 0;
    end if;
    for i in planes'range loop               -- solution planes into vectors
      declare
        plaloc : constant Standard_Complex_Matrices.Matrix
               := Localize(locmap,planes(i).all);
        solloc : constant Standard_Complex_Vectors.Vector
               := Column_Vector_Rep(locmap,plaloc);
        locsol : Standard_Complex_Vectors.Vector(1..solloc'last+2+addmix);
      begin
        locsol(solloc'range) := solloc;
        if s_mode = 0
         then locsol(locsol'last-1) := Create(0.0); -- start value for s is 0
         else locsol(locsol'last-1) := Create(1.0); -- start value for s is 1
        end if;
        locsol(locsol'last)   := Create(0.0);   -- start value for t is 0
        if addmix = 1
         then locsol(locsol'last-2) := Create(1.0);   -- additional s-value
        end if;
        locsols(i) := new Standard_Complex_Vectors.Vector'(locsol);
        if outlog
         then evaloc := Eval(lochom,locsol);
              put(file,"Residual of localized plane at start :");
              put(file,Max_Norm(evaloc),3); put_line(file,".");
        end if;
      end;
    end loop;
    Determinantal_Polynomial_Continuation
      (file,lochom,integer32(conpar),locsols,report,outlog);
    for i in locsols'range loop             -- solution vectors into planes
      planes(i).all
        := Column_Matrix_Rep(locmap,locsols(i)(1..locsols(i)'last-2-addmix));
      if outlog then
        evaloc := Eval(lochom,locsols(i).all);
        put(file,"Residual of localized plane at target :");
        put(file,Max_Norm(evaloc),3); put_line(file,".");
        put_line(file,"The computed plane : ");
        put(file,planes(i).all,2);
      end if;
      declare
        plavec : constant Standard_Complex_Vectors.Vector
               := Column_Vector_Rep(nd.top,nd.bottom,planes(i).all);
        solvec : Standard_Complex_Vectors.Vector
                   (plavec'first..plavec'last+2+addmix);
      begin
        solvec(plavec'range) := plavec;
        solvec(solvec'last)   := Create(1.0);   -- target value for t
        solvec(solvec'last-1) := locsols(i)(locsols(i)'last-1);
        if addmix = 1
         then solvec(solvec'last-2)   := locsols(i)(locsols(i)'last-2);
        end if;
        if outlog then
          evahom := Eval(homsys,solvec);
          put(file,"Residual of target solution :");
          put(file,Max_Norm(evahom),3); put_line(file,".");
        end if;
      end;
    end loop;
    if outlog then
      declare
        dim : constant integer32 := locsols(1)'last-2-addmix;
      begin
        Verify_Determinants(file,dim,nd,xpm,locmap,locsols,s,ip);
      end;
    end if;
    Clear(lochom); Clear(locsols);
  end Quantum_Trace_Paths;

end Pieri_Continuation;
