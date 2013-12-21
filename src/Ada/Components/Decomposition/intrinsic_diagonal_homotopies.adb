with integer_io;                        use integer_io;
with Timing_Package;                    use Timing_Package;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Random_Matrices;
with Standard_Natural_Vectors;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Standard_Condition_Tables;         use Standard_Condition_Tables;
with Standard_Continuation_Data;        use Standard_Continuation_Data;
with Continuation_Parameters;
with Drivers_for_Poly_Continuation;     use Drivers_for_Poly_Continuation;
with Standard_Plane_Representations;    use Standard_Plane_Representations;
with Standard_Moving_Planes;            use Standard_Moving_Planes;
with Standard_Cascading_Planes;         use Standard_Cascading_Planes;
with Standard_Intrinsic_Solutions;      use Standard_Intrinsic_Solutions;
with Standard_Intrinsic_Newton;         use Standard_Intrinsic_Newton;
with Standard_Intrinsic_Continuation;   use Standard_Intrinsic_Continuation;
with Standard_Singular_Filters;         use Standard_Singular_Filters;
with Extrinsic_Diagonal_Homotopies;     use Extrinsic_Diagonal_Homotopies;

package body Intrinsic_Diagonal_Homotopies is

  sia_size : constant natural := 20000;

-- AUXILIARIES :

  procedure Refine_Roots
              ( file : in file_type; f : in Eval_Poly_Sys;
                jm : in Eval_Jaco_Mat; p : in Matrix;
                sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Applies Newton's method to the solutions in intrinsic format.

    res,res_last : Solution_List;
    epsax : constant double_float := 1.0E-13;
    epsrx : constant double_float := 1.0E-13;
    epsaf : constant double_float := 1.0E-13;
    epsrf : constant double_float := 1.0E-13;
    tol_accept : constant double_float := 1.0E-8;
    incax,incrx,resaf,resrf : double_float;
    cnt,nit : natural;
    fail : boolean;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    table_res : Standard_Natural_Vectors.Vector(0..15) := Create(15);
    table_cor : Standard_Natural_Vectors.Vector(0..15) := Create(15);

  begin
    put(file,"Refining "); put(file,Length_Of(sols),1);
    put(file," solutions ...");
    cnt := 0;
    for i in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
     -- put(file,"Refining solution "); put(file,i,1); put_line(file,"...");
     -- Affine_LU_Newton(file,f,jm,p,ls.v,epsax,epsrx,epsaf,epsrf,
     --                  incax,incrx,resaf,resrf,nit,5,fail);
      Affine_LU_Newton(f,jm,p,ls.v,epsax,epsrx,epsaf,epsrf,
                       incax,incrx,resaf,resrf,nit,5,fail);
      ls.err := incax; ls.res := resaf;
      Update_Corrector(table_cor,ls.all);
      Update_Residuals(table_res,ls.all);
      if ls.err < tol_accept or ls.res < tol_accept
       then Append(res,res_last,ls.all);
            cnt := cnt + 1;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    put(file," kept "); put(file,cnt,1); put_line(file," solutions.");
    Write_Corrector(file,table_cor);
    Write_Residuals(file,table_res);
    Clear(sols);
    sols := res;
  end Refine_Roots;

  generic
    with function f ( x : Vector ) return Vector;
    with function jf ( x : Vector ) return Matrix;
  procedure Generic_Refine_Roots
              ( file : in file_type; n : in natural;
                p : in Matrix; sols : in out Solution_List );

  -- DESCRIPTION :
  --   Generic version of the Refine_Roots from above.
  --   Instead of f and jf, it needs n (#equations in f) on entry.

  procedure Generic_Refine_Roots
              ( file : in file_type; n : in natural;
                p : in Matrix; sols : in out Solution_List ) is

    res,res_last : Solution_List;
    epsax : constant double_float := 1.0E-13;
    epsrx : constant double_float := 1.0E-13;
    epsaf : constant double_float := 1.0E-13;
    epsrf : constant double_float := 1.0E-13;
    tol_accept : constant double_float := 1.0E-8;
    incax,incrx,resaf,resrf : double_float;
    cnt,nit : natural;
    fail : boolean;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    table_cor : Standard_Natural_Vectors.Vector(0..15) := Create(15);
    table_res : Standard_Natural_Vectors.Vector(0..15) := Create(15);

   -- procedure Newton is new Reporting_Affine_LU_Newton(f,jf);
    procedure Newton is new Silent_Affine_LU_Newton(f,jf);

  begin
    put(file,"Refining "); put(file,Length_Of(sols),1);
    put(file," solutions ...");
    cnt := 0;
    for i in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
     -- put(file,"Refining solution "); put(file,i,1); put_line(file,"...");
     -- Newton(file,n,p,ls.v,epsax,epsrx,epsaf,epsrf,
     --        incax,incrx,resaf,resrf,nit,5,fail);
      Newton(n,p,ls.v,epsax,epsrx,epsaf,epsrf,
             incax,incrx,resaf,resrf,nit,5,fail);
      ls.err := incax; ls.res := resaf;
      Update_Corrector(table_cor,ls.all);
      Update_Residuals(table_res,ls.all);
      if ls.err < tol_accept or ls.res < tol_accept
       then Append(res,res_last,ls.all);
            cnt := cnt + 1;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    put(file," kept "); put(file,cnt,1); put_line(file," solutions.");
    Write_Corrector(file,table_cor);
    Write_Residuals(file,table_res);
    Clear(sols);
    sols := res;
  end Generic_Refine_Roots;

  procedure Track_Paths
              ( file : in file_type; f : in Eval_Poly_Sys;
                jm : in Eval_Jaco_Mat;
                start_plane,target_plane : in Matrix;
                sols : in out Solution_List ) is

    timer : Timing_Widget;
    pp : Pred_Pars;
    cp : Corr_Pars;
    gamma : Complex_Number := Random1;
    sials,tmp : Solu_Info_Array_List;

    function Path ( t : Complex_Number ) return Matrix is
    begin
      return Moving_Plane(start_plane,target_plane,gamma,t);
    end Path;
    procedure Cont is new Reporting_LU_Continue(Path);

  begin
    pp := Continuation_Parameters.Create_for_Path;
    cp := Continuation_Parameters.Create_for_Path;
    Set_Continuation_Parameter(sols,Create(0.0));
    sials := Create(sols,sia_size);
    tmp := sials;
    tstart(timer);
    while not Is_Null(tmp) loop
      declare
        sia : Link_to_Solu_Info_Array := Head_Of(tmp);
      begin
        Cont(file,f,jm,sia.all,pp,cp);
        Set_Head(tmp,sia);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Continuation in Intrinsic Coordinates");
    Deep_Clear(sols);
    sols := Concat(sials);
  end Track_Paths;

  generic
    with function f ( x : Vector ) return Vector;
    with function jf ( x : Vector ) return Matrix;
  procedure Generic_Track_Paths 
              ( file : in file_type; ne,nv : in natural;
                start_plane,target_plane : in Matrix;
                sols : in out Solution_List );

  -- DESCRIPTION :
  --   Generic version of the Track_Paths from above.
  --   Instead of f and jf, it needs ne (#equations)
  --   and nv (#variables) as input parameters.

  procedure Generic_Track_Paths 
              ( file : in file_type; ne,nv : in natural;
                start_plane,target_plane : in Matrix;
                sols : in out Solution_List ) is

    timer : Timing_Widget;
    pp : Pred_Pars;
    cp : Corr_Pars;
    gamma : Complex_Number := Random1;
    sials,tmp : Solu_Info_Array_List;

    function Path ( t : Complex_Number ) return Matrix is
    begin
      return Moving_Plane(start_plane,target_plane,gamma,t);
    end Path;
    procedure Cont is new G_Reporting_LU_Continue(f,jf,Path);

  begin
    pp := Continuation_Parameters.Create_for_Path;
    cp := Continuation_Parameters.Create_for_Path;
    Set_Continuation_Parameter(sols,Create(0.0));
    sials := Create(sols,sia_size);
    tmp := sials;
    tstart(timer);
    while not Is_Null(tmp) loop
      declare
        sia : Link_to_Solu_Info_Array := Head_Of(tmp);
      begin
        Cont(file,ne,nv,sia.all,pp,cp);
        Set_Head(tmp,sia);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Continuation in Intrinsic Coordinates");
    Deep_Clear(sols);
    sols := Concat(sials);
  end Generic_Track_Paths;

  procedure Extract_Halves ( n : in natural; x : in Vector;
                             x1,x2 : out Vector ) is

  -- DESCRIPTION :
  --   Puts the first n elements of x in x1 and the last n elements in x2.

  begin
    for i in 1..n loop
      x1(i) := x(i);
      x2(i) := x(i+n);
    end loop;
  end Extract_Halves;

  function Stack_Vectors
             ( m,n1,n2 : natural; y1,y2 : Vector ) return Vector is

  -- DESCRIPTION :
  --    Returns a vector of range 1..m, obtained by stacking y2
  --    (of range 1..n2) after y1 (of range 1..n1).

    res : Vector(1..m);

  begin
    for i in 1..n1 loop
      res(i) := y1(i);
    end loop;
    for i in 1..n2 loop
      res(n1+i) := y2(i);
    end loop;
    return res;
  end Stack_Vectors;

  function Stack_Matrices
              ( n,n2,m,nA,nB : in natural; A,B : Matrix ) return Matrix is

  -- DESCRIPTION :
  --   Returns a matrix with m rows and n2 columns with on its
  --   diagonal the nA-by-n matrix A and the nB-by-n matrix B.

    res : Matrix(1..m,1..n2);

  begin
    for i in 1..nA loop
      for j in 1..n loop
        res(i,j) := A(i,j);
        res(i,j+n) := Create(0.0);
      end loop;
    end loop;
    for i in 1..nB loop
      for j in 1..n loop
        res(nA+i,j) := Create(0.0);
        res(nA+i,j+n) := B(i,j);
      end loop;
    end loop;
    return res;
  end Stack_Matrices;

-- TARGET ROUTINES :

  function Minimal_Intersection_Dimension
             ( n,a,b : natural ) return natural is

    res : natural := a + b - n;

  begin
    if res < 0
     then return 0;
     else return res;
    end if;
  end Minimal_Intersection_Dimension;

  function Collapse ( n : natural; p : Matrix ) return Matrix is

  -- DESCRIPTION :
  --   The matrix on return has the first n rows of p.

     res : Matrix(1..n,p'range(2));

  begin
    for i in 1..n loop
      for j in p'range(2) loop
        res(i,j) := p(i,j);
      end loop;
    end loop;
    return res;
  end Collapse;

  procedure Singular_Filter
              ( file : in file_type;
                sols : in out Solution_List; tol : in double_float ) is

    len : constant natural := Length_Of(sols);
    sing,regu : Solution_List;

  begin
    R_Singular_Filter(file,sols,tol,sing,regu);
   -- S_Singular_Filter(sols,tol,sing,regu);
    put(file,"Tested "); put(file,len,1); put(file," solutions, found ");
    put(file,Length_Of(regu),1); put(file," regular and ");
    put(file,Length_Of(sing),1); put_line(file," singular.");
    Clear(sols);
    sols := regu;
  end Singular_Filter;

  procedure Intrinsic_Diagonal_Homotopy
              ( file : in file_type; a,b : in natural;
                p1,p2 : in Poly_Sys; sols1,sols2 : in Solution_List;
                plane1,plane2 : in Matrix;
                sols : out Solution_List; plane : out Matrix ) is

    n : constant natural := Number_of_Unknowns(p1(p1'first));
    n2 : constant natural := 2*n;
    apb : constant natural := a+b;
    m : constant natural := n2-apb;
    h0 : constant natural := Minimal_Intersection_Dimension(n,a,b);
    AA : Matrix(1..apb,1..n) := Standard_Random_Matrices.Random_Matrix(apb,n);
    dA : Matrix(1..apb,1..n2) := Double_to_Diagonal(AA);
    BB : Matrix(1..apb,1..n) := Standard_Random_Matrices.Random_Matrix(apb,n);
    CC : Matrix(1..n,1..n2) := Standard_Random_Matrices.Random_Matrix(n,n2);
    d : Vector(1..n) := Random_Vector(1,n);
    offset : Vector(1..n2) := Compute_Offset(CC,d);
    target : Matrix(1..apb,0..n2) := Target_Space(n,n2,apb,b,dA,BB,CC,d);
    gentar : Matrix(1..n2,0..m) := Generators(target);
    tarpla : Matrix(1..n2,0..m) := Orthogonalize(gentar);
    stapla : Matrix(1..n2,0..m) := Start_Space(plane1,plane2);
    targt1 : Matrix(1..apb,0..n2);
    gentr1,tarpl1 : Matrix(1..n2,0..m);
    f : Poly_Sys(1..m) := Product(n,n,p1,p2);
    ef : Eval_Poly_Sys(f'range) := Create(f);
    sjm : Jaco_Mat(f'range,1..n2) := Create(f);
    sjf : Eval_Jaco_Mat(sjm'range(1),sjm'range(2)) := Create(sjm);
    psols : Solution_List := Product(sols1,sols2);
   -- esols : Solution_List;
    timer : Timing_Widget;
    tolsing : constant double_float := 1.0E-8;

  begin
    tstart(timer);
    Shift_Offset(tarpla,offset);
    targt1 := Target_Space(n,n2,apb,b-1,dA,BB,CC,d);
    gentr1 := Generators(targt1);
    tarpl1 := Orthogonalize(gentr1);
    Shift_Offset(tarpl1,offset);
    Intersect(target,targt1,tarpla,tarpl1);
   -- new_line(file);
   -- put_line(file,"Path Tracking to Start the Cascade");
   -- new_line(file);
   -- put_line(file,"Refining the roots at the start...");
    Refine_Roots(file,ef,sjf,stapla,psols);
    new_line(file);
    put(file,"Tracking "); put(file,Length_Of(psols),1);
    put_line(file," to start the cascade ...");
    Track_Paths(file,ef,sjf,stapla,tarpla,psols);
    new_line(file);
    put_line(file,"Refining the roots at the target...");
    Refine_Roots(file,ef,sjf,tarpla,psols);
    new_line(file);
    put(file,"Starting the Cascade at Dimension ");
    put(file,b-1,1); put_line(file,".");
   -- new_line(file);
   -- put_line(file,"Refining the roots at the start...");
   -- Refine_Roots(file,ef,sjf,tarpla,psols);
    Set_Continuation_Parameter(psols,Create(0.0));
    put(file,"Tracking "); put(file,Length_Of(psols),1);
    put_line(file," paths to the top dimensional intersection component ...");
    Track_Paths(file,ef,sjf,tarpla,tarpl1,psols);
   -- new_line(file);
   -- put_line(file,"Refining the roots at the target...");
    Refine_Roots(file,ef,sjf,tarpl1,psols);
    Singular_Filter(file,psols,tolsing);
    if not Is_Null(psols) then
      target := targt1;
      tarpla := tarpl1;
     -- running the cascade down to h0
      for h in reverse h0..b-2 loop
        new_line(file);
        put(file,"Starting the Cascade at Dimension ");
        put(file,h,1); put_line(file,".");
        new_line(file);
        targt1 := Target_Space(n,n2,apb,h,dA,BB,CC,d);
        gentr1 := Generators(targt1);
        tarpl1 := Orthogonalize(gentr1);
        Shift_Offset(tarpl1,offset);
        gentar := tarpla;
        Intersect(target,targt1,tarpla,tarpl1);
        Transform(psols,gentar,tarpla);
       -- new_line(file);
       -- put_line(file,"Refining the roots at the start...");
        Refine_Roots(file,ef,sjf,tarpla,psols);
        Set_Continuation_Parameter(psols,Create(0.0));
        put(file,"Tracking "); put(file,Length_Of(psols),1);
        put_line(file," paths to intersection components ...");
        Track_Paths(file,ef,sjf,tarpla,tarpl1,psols);
       -- new_line(file);
       -- put_line(file,"Refining the roots at the target...");
        Refine_Roots(file,ef,sjf,tarpl1,psols);
        Singular_Filter(file,psols,tolsing);
        exit when Is_Null(psols);
        target := targt1; tarpla := tarpl1;
      end loop;
    end if;
    plane := Collapse(n,tarpl1);
    if not Is_Null(psols) 
     then sols := psols;
   --       esols := Expand(psols,tarpl1);
   --       put_line(file,"THE SOLUTIONS at the end : ");
   --       put(file,Length_Of(esols),Head_Of(esols).n,esols);
   --       Clear(esols);
    end if;
    Clear(f); Clear(ef); Clear(sjm); Clear(sjf);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Intrinsic Diagonal Homotopy Cascade");
  end Intrinsic_Diagonal_Homotopy;

  procedure Generic_Diagonal_Homotopy
              ( file : in file_type; nefA,nefB,n,a,b : in natural;
                sols1,sols2 : in Solution_List; plane1,plane2 : in Matrix;
                sols : out Solution_List; plane : out Matrix ) is

    n2 : constant natural := 2*n;
    apb : constant natural := a+b;
    m : constant natural := n2-apb;
    h0 : constant natural := Minimal_Intersection_Dimension(n,a,b);
    AA : Matrix(1..apb,1..n) := Standard_Random_Matrices.Random_Matrix(apb,n);
    dA : Matrix(1..apb,1..n2) := Double_to_Diagonal(AA);
    BB : Matrix(1..apb,1..n) := Standard_Random_Matrices.Random_Matrix(apb,n);
    CC : Matrix(1..n,1..n2) := Standard_Random_Matrices.Random_Matrix(n,n2);
    d : Vector(1..n) := Random_Vector(1,n);
    offset : Vector(1..n2) := Compute_Offset(CC,d);
    target : Matrix(1..apb,0..n2) := Target_Space(n,n2,apb,b,dA,BB,CC,d);
    gentar : Matrix(1..n2,0..m) := Generators(target);
    tarpla : Matrix(1..n2,0..m) := Orthogonalize(gentar);
    stapla : Matrix(1..n2,0..m) := Start_Space(plane1,plane2);
    targt1 : Matrix(1..apb,0..n2);
    gentr1,tarpl1 : Matrix(1..n2,0..m);
    psols : Solution_List := Product(sols1,sols2);
   -- esols : Solution_List;
    timer : Timing_Widget;
    tolsing : constant double_float := 1.0E-8;

    function Product ( x : Vector ) return Vector is

    -- DESCRIPTION :
    --   Returns the result of the evaluation of the functions at x.

      res : Vector(1..m);
      x1,x2 : Vector(1..n);
      y1 : Vector(1..nefA);
      y2 : Vector(1..nefB);

    begin
      Extract_Halves(n,x,x1,x2);
      y1 := fA(x1);
      y2 := fB(x2);
      res := Stack_Vectors(m,nefA,nefB,y1,y2);
      return res;
    end Product;

    function Jacobi ( x : Vector ) return Matrix is

    -- DESCRIPTION :
    --   Returns the value of the Jacobi matrix at the product.

      res : Matrix(1..m,1..n2);
      x1,x2 : Vector(1..n);
      y1 : Matrix(1..nefA,1..n);
      y2 : Matrix(1..nefB,1..n);

    begin
      Extract_Halves(n,x,x1,x2);
      y1 := jfA(x1);
      y2 := jfB(x2);
      res := Stack_Matrices(n,n2,m,nefA,nefB,y1,y2);
      return res;
    end Jacobi;

    procedure Track is new Generic_Track_Paths(Product,Jacobi);
    procedure Refine is new Generic_Refine_Roots(Product,Jacobi);

  begin
    tstart(timer);
    Shift_Offset(tarpla,offset);
    targt1 := Target_Space(n,n2,apb,b-1,dA,BB,CC,d);
    gentr1 := Generators(targt1);
    tarpl1 := Orthogonalize(gentr1);
    Shift_Offset(tarpl1,offset);
    Intersect(target,targt1,tarpla,tarpl1);
   -- new_line(file);
   -- put_line(file,"Path Tracking to Start the Cascade");
   -- new_line(file);
   -- put_line(file,"Refining the roots at the start...");
    Refine(file,m,stapla,psols);
    new_line(file);
    put(file,"Tracking "); put(file,Length_Of(psols),1);
    put_line(file," paths to start the cascade ...");
    Track(file,m,n2,stapla,tarpla,psols);
   -- new_line(file);
   -- put_line(file,"Refining the roots at the target...");
    Refine(file,m,tarpla,psols);
    new_line(file);
    put(file,"Starting the Cascade at Dimension ");
    put(file,b-1,1); put_line(file,".");
   -- new_line(file);
   -- put_line(file,"Refining the roots at the start...");
   -- Refine(file,m,tarpla,psols);
    Set_Continuation_Parameter(psols,Create(0.0));
    put(file,"Tracking "); put(file,Length_Of(psols),1);
    put_line(file," to the top dimensional intersection component ...");
    Track(file,m,n2,tarpla,tarpl1,psols);
   -- new_line(file);
   -- put_line(file,"Refining the roots at the target...");
    Refine(file,m,tarpl1,psols);
    Singular_Filter(file,psols,tolsing);
    if not Is_Null(psols) then
      target := targt1;
      tarpla := tarpl1;
     -- running the cascade down to h0
      for h in reverse h0..b-2 loop
        new_line(file);
        put(file,"Running the Cascade at Dimension ");
        put(file,h,1); put_line(file,".");
        new_line(file);
        targt1 := Target_Space(n,n2,apb,h,dA,BB,CC,d);
        gentr1 := Generators(targt1);
        tarpl1 := Orthogonalize(gentr1);
        Shift_Offset(tarpl1,offset);
        gentar := tarpla;
        Intersect(target,targt1,tarpla,tarpl1);
        Transform(psols,gentar,tarpla);
       -- new_line(file);
       -- put_line(file,"Refining the roots at the start...");
        Refine(file,m,tarpla,psols);
        Set_Continuation_Parameter(psols,Create(0.0));
        put(file,"Tracking "); put(file,Length_Of(psols),1);
        put_line(file," to intersection components ...");
        Track(file,m,n2,tarpla,tarpl1,psols);
       -- new_line(file);
       -- put_line(file,"Refining the roots at the target...");
        Refine(file,m,tarpl1,psols);
        Singular_Filter(file,psols,tolsing);
        exit when Is_Null(psols);
        target := targt1; tarpla := tarpl1;
      end loop;
    end if;
    plane := Collapse(n,tarpl1);
    if not Is_Null(psols) 
     then sols := psols;
   --       esols := Expand(psols,tarpl1);
   --       put_line(file,"THE SOLUTIONS at the end : ");
   --       put(file,Length_Of(esols),Head_Of(esols).n,esols);
   --       Clear(esols);
    end if;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Intrinsic Diagonal Homotopy Cascade");
  end Generic_Diagonal_Homotopy;

end Intrinsic_Diagonal_Homotopies;
