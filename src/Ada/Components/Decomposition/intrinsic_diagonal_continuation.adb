with Timing_Package;                    use Timing_Package;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Random_Matrices;
with Standard_Natural_Vectors;
with Symbol_Table;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Standard_Condition_Tables;         use Standard_Condition_Tables;
with Standard_Homotopy;
with Standard_Continuation_Data;        use Standard_Continuation_Data;
with Standard_IncFix_Continuation;      use Standard_IncFix_Continuation;
with Continuation_Parameters;           use Continuation_Parameters;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets;
with Planes_and_Polynomials;            use Planes_and_Polynomials;
with Standard_Plane_Representations;    use Standard_Plane_Representations;
with Standard_Point_Coordinates;        use Standard_Point_Coordinates;
with Standard_Moving_Planes;            use Standard_Moving_Planes;
with Standard_Plane_Operations;         use Standard_Plane_Operations;
with Standard_Cascading_Planes;         use Standard_Cascading_Planes;
--with Standard_Rescaling_Coordinates;    use Standard_Rescaling_Coordinates;
with Standard_Intrinsic_Solutions;      use Standard_Intrinsic_Solutions;
with Standard_Intrinsic_Newton;         use Standard_Intrinsic_Newton;
with Standard_Intrinsic_Continuation;   use Standard_Intrinsic_Continuation;
with Standard_Solution_Splitters;       use Standard_Solution_Splitters;
with Standard_Diagonal_Polynomials;     use Standard_Diagonal_Polynomials;
with Standard_Diagonal_Solutions;       use Standard_Diagonal_Solutions;

package body Intrinsic_Diagonal_Continuation is

  sia_size : constant natural32 := 20000;

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
    cnt,nit : natural32;
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
    Write_Corrector_Table(file,table_cor);
    Write_Residuals_Table(file,table_res);
    Clear(sols);
    sols := res;
  end Refine_Roots;

  generic

    with function f ( x : Vector ) return Vector;
    with function jf ( x : Vector ) return Matrix;

  procedure Generic_Refine_Roots
              ( file : in file_type; n : in integer32;
                p : in Matrix; sols : in out Solution_List );

  -- DESCRIPTION :
  --   Generic version of the Refine_Roots from above.
  --   Instead of f and jf, it needs n (#equations in f) on entry.

  procedure Generic_Refine_Roots
              ( file : in file_type; n : in integer32;
                p : in Matrix; sols : in out Solution_List ) is

    res,res_last : Solution_List;
    epsax : constant double_float := 1.0E-13;
    epsrx : constant double_float := 1.0E-13;
    epsaf : constant double_float := 1.0E-13;
    epsrf : constant double_float := 1.0E-13;
    tol_accept : constant double_float := 1.0E-8;
    incax,incrx,resaf,resrf : double_float;
    cnt,nit : natural32;
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
      Newton(natural32(n),p,ls.v,epsax,epsrx,epsaf,epsrf,
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
    Write_Corrector_Table(file,table_cor);
    Write_Residuals_Table(file,table_res);
    Clear(sols);
    sols := res;
  end Generic_Refine_Roots;

  procedure Intrinsic_Track_Paths
              ( file : in file_type; report : in boolean;
                f : in Eval_Poly_Sys; jm : in Eval_Jaco_Mat;
                start_plane,target_plane : in Matrix;
                sols : in out Solution_List ) is

    timer : Timing_Widget;
    pp : Pred_Pars;
    cp : Corr_Pars;
    gamma : constant Complex_Number := Random1;
    sials,tmp : Solu_Info_Array_List;

    function Path ( t : Complex_Number ) return Matrix is
    begin
      return Moving_Plane(start_plane,target_plane,gamma,t);
     -- the random gamma does not seem to matter much...
     -- return Moving_Plane(start_plane,target_plane,t);
    end Path;
    procedure S_Cont is new Silent_Affine_LU_Continue(Path);
    procedure R_Cont is new Reporting_Affine_LU_Continue(Path);
   -- procedure S_Cont is new Silent_QR_Continue(Path);
   -- procedure R_Cont is new Reporting_QR_Continue(Path);

  begin
    pp := Continuation_Parameters.Create_for_Path;
    cp := Continuation_Parameters.Create_for_Path;
    Set_Continuation_Parameter(sols,Create(0.0));
    sials := Create(sols,sia_size);
    tmp := sials;
    tstart(timer);
    while not Is_Null(tmp) loop
      declare
        sia : constant Link_to_Solu_Info_Array := Head_Of(tmp);
      begin
        if report
         then R_Cont(file,f,jm,sia.all,pp,cp);
         else S_Cont(f,jm,sia.all,pp,cp);
        end if;
        Set_Head(tmp,sia);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Continuation in Intrinsic Coordinates");
    Deep_Clear(sols);
    sols := Concat(sials);
  end Intrinsic_Track_Paths;

  procedure Recentered_Path_Tracking
              ( file : in file_type; report : in boolean;
                f : in Eval_Poly_Sys; jm : in Eval_Jaco_Mat;
                start_plane,target_plane : in Matrix;
                sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   This procedure is similar to "Intrinsic_Track_Paths"
  --   except that local intrinsic coordinates are used.

  -- ON ENTRY :
  --   file     for writing diagnostics and extra output;
  --   report   if false, then the path trackers will be silent;
  --   f        polynomial system ready for evaluation;
  --   jm       Jacobian matrix of f;
  --   start_plane is the starting plane at the continuation;
  --   target_plane is the target plane;
  --   sols     solutions in intrinsic coordinates at the start plane.

  -- ON RETURN :
  --   sols     solutions in intrinsic coordinates at the target plane.

    timer : Timing_Widget;
    pp : Pred_Pars;
    cp : Corr_Pars;
    sials,tmp : Solu_Info_Array_List;
   -- fail : boolean;
    esols : Solution_List := Expand(sols,start_plane);

  begin
   -- Check_Orthonormality(target_plane,1.0E-12,fail);
   -- if fail
   --  then put_line("Target plane did not pass the orthonormality test!");
   --  else put_line("Target plane passed the orthonormality test.");
   -- end if;
    pp := Continuation_Parameters.Create_for_Path;
    cp := Continuation_Parameters.Create_for_Path;
    Set_Continuation_Parameter(esols,Create(0.0));
    sials := Create(esols,sia_size);
    tmp := sials;
    tstart(timer);
    while not Is_Null(tmp) loop
      declare
        sia : constant Link_to_Solu_Info_Array := Head_Of(tmp);
      begin
        if report
         then Reporting_Local_LU_Continue
                (file,f,jm,start_plane,target_plane,false,sia.all,pp,cp);
         else Silent_Local_LU_Continue
                (f,jm,start_plane,target_plane,false,sia.all,pp,cp);
        end if;
        Set_Head(tmp,sia);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Continuation in Local Intrinsic Coordinates");
    Deep_Clear(esols);
    esols := Concat(sials);
    sols := Project(esols,target_plane);
  end Recentered_Path_Tracking;

  generic

    with function f ( x : Vector ) return Vector;
    with function jf ( x : Vector ) return Matrix;

  procedure Generic_Track_Paths 
              ( file : in file_type; ne,nv : in integer32;
                start_plane,target_plane : in Matrix;
                sols : in out Solution_List );

  -- DESCRIPTION :
  --   Generic version of the Track_Paths from above.
  --   Instead of f and jf, it needs ne (#equations)
  --   and nv (#variables) as input parameters.

  procedure Generic_Track_Paths 
              ( file : in file_type; ne,nv : in integer32;
                start_plane,target_plane : in Matrix;
                sols : in out Solution_List ) is

    timer : Timing_Widget;
    pp : Pred_Pars;
    cp : Corr_Pars;
    gamma : constant Complex_Number := Random1;
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
        sia : constant Link_to_Solu_Info_Array := Head_Of(tmp);
      begin
        Cont(file,natural32(ne),natural32(nv),sia.all,pp,cp);
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

  procedure Extract_Halves ( n : in integer32; x : in Vector;
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
             ( m,n1,n2 : integer32; y1,y2 : Vector ) return Vector is

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
             ( n,n2,m,nA,nB : in integer32; A,B : Matrix ) return Matrix is

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
             ( n,a,b : integer32 ) return integer32 is

    res : constant integer32 := a + b - n;

  begin
    if res < 0
     then return 0;
     else return res;
    end if;
  end Minimal_Intersection_Dimension;

  function Collapse ( n : integer32; p : Matrix ) return Matrix is

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

    len : constant natural32 := Length_Of(sols);
    sing,regu : Solution_List;

  begin
    Reporting_Singular_Filter(file,sols,tol,sing,regu);
   -- Silent_Singular_Filter(sols,tol,sing,regu);
    put(file,"Tested "); put(file,len,1); put(file," solutions, found ");
    put(file,Length_Of(regu),1); put(file," regular and ");
    put(file,Length_Of(sing),1); put_line(file," singular.");
    Clear(sols);
    sols := regu;
  end Singular_Filter;

  procedure Scan_Failed_Paths
              ( file : in file_type; sols : in Solution_List;
                p1,p2 : in Matrix ) is

    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    cnt : natural32 := 0;
    p : Matrix(p1'range(1),p2'range);

  begin
    new_line(file);
    put_line(file,"THE FAILED SOLUTION PATHS : ");
    for i in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      if REAL_PART(ls.t) < 1.0 then
        cnt := cnt + 1;
        put(file,"solution "); put(file,i,1);
        put(file," is failed path ");
        put(file,cnt,1); put_line(file," :");
        put(file,ls.all);
        new_line(file);
        for j in p'range(1) loop
          for k in p'range(2) loop
            p(j,k) := p1(j,k)*(Create(1.0)-ls.t) + p2(j,k)*ls.t;
          end loop;
        end loop;
        declare
          x : constant Vector := Affine_Expand(ls.v,p);
        begin
          put_line(file,"The expanded vector : ");
          put_line(file,x);
        end;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Scan_Failed_Paths;

--  procedure Copy_Offset ( p : in Matrix; q : out Matrix ) is
--
--  -- DESCRIPTION :
--  --   Copies offset vector from p to q.
--
--  begin
--    for i in q'range(1) loop
--      q(i,0) := p(i,0);
--    end loop;
--  end Copy_Offset;

--  procedure Copy_Directions ( p : in Matrix; q : out Matrix ) is
--
--  -- DESCRIPTION :
--  --   Copies directions from p to q.
--
--  begin
--    for i in q'range(1) loop
--      for j in 1..q'last(2) loop
--        q(i,j) := p(i,j);
--      end loop;
--    end loop;
--  end Copy_Directions;

--  procedure Copy_Direction ( p : in Matrix; k : in natural;
--                             q : out Matrix ) is
--
--  -- DESCRIPTION :
--  --   Copies the kth direction vector from p to q.
--
--  -- REQUIRED : k <= p'last(2) = q'last(2).
--
--  begin
--    for i in q'range(1) loop
--      q(i,k) := p(i,k);
--    end loop;
--  end Copy_Direction;

--  procedure Start_Intrinsic_Cascade_the_Long_Way
--              ( file : in file_type; report : in boolean;
--                ef : in Eval_Poly_Sys; sjf : in Eval_Jaco_Mat;
--                stapla,tarpla : in Matrix; psols : in out Solution_List ) is
--
--  -- DESCRIPTION :
--  --   This is an experimental procedure to overcome the problems
--  --   with starting the cascade using intrinsic coordinates.
--  --   The suffix "the_Long_Way" reflects the attempt to move over
--  --   the directions one after the other, which failed badly.
--
--    len1,len2 : natural;
--    backsols : Solution_List;
--    wrkpla,endpla : Matrix(stapla'range(1),stapla'range(2));
--
--  begin
--   -- new_line(file);
--   -- put_line(file,"Path Tracking to Start the Cascade");
--   -- new_line(file);
--   -- put_line(file,"Refining the roots at the start...");
--    Refine_Roots(file,ef,sjf,stapla,psols);
--    wrkpla := stapla;
--    Copy_Directions(stapla,endpla); -- Copy_Offset(stapla,endpla);
--    Copy_Offset(tarpla,endpla);     -- Copy_Directions(tarpla,endpla);
--    for i in 0..tarpla'last(2) loop
--      len1 := Length_Of(psols);
--      new_line(file);
--      put(file,"Tracking "); put(file,len1,1);
--      put(file," paths to start the cascade, stage ");
--      put(file,i,1); put_line(file," ...");
--      Intrinsic_Track_Paths(file,report,ef,sjf,wrkpla,endpla,psols);
--      Copy(psols,backsols);
--      new_line(file);
--      put_line(file,"Refining the roots at the target...");
--      Refine_Roots(file,ef,sjf,endpla,psols);
--      len2 := Length_Of(psols);
--      if len1 = len2
--       then put(file,"No paths lost in Start of the Cascade, stage ");
--            put(file,i,1); put_line(file,".  OK");
--       else put(file,len1-len2,1);
--            put(file," paths lost in Start of the Cascade, stage ");
--            put(file,i,1); put_line(file,".  BUG!!!");
--            Scan_Failed_Paths(file,backsols,wrkpla,endpla);
--      end if;
--      wrkpla := endpla;
--     -- Copy_Directions(tarpla,endpla); --Copy_Offset(tarpla,endpla);
--      if i < tarpla'last(2)
--       then Copy_Direction(tarpla,i+1,endpla);
--      end if;
--    end loop;
--  end Start_Intrinsic_Cascade_the_Long_Way;

  procedure Start_Intrinsic_Cascade
              ( file : in file_type; report : in boolean;
                ef : in Eval_Poly_Sys; sjf : in Eval_Jaco_Mat;
                stapla,tarpla : in Matrix; psols : in out Solution_List ) is

  -- DESCRIPTION :
  --   This is the homotopy used to start up the cascade.
  --   There should be no diverging paths in this homotopy!

    len1,len2 : natural32;
    backsols : Solution_List;

  begin
    new_line(file);
    put_line(file,"Path Tracking to Start the Cascade");
    new_line(file);
    put_line(file,"Refining the roots at the start...");
    Refine_Roots(file,ef,sjf,stapla,psols);
    len1 := Length_Of(psols);
    new_line(file);
    put(file,"Tracking "); put(file,len1,1);
    put_line(file," paths to start the cascade ... ");
   -- Intrinsic_Track_Paths(file,report,ef,sjf,stapla,tarpla,psols);
    Recentered_Path_Tracking(file,report,ef,sjf,stapla,tarpla,psols);
    Copy(psols,backsols);
    new_line(file);
    put_line(file,"Refining the roots at the target...");
    Refine_Roots(file,ef,sjf,tarpla,psols);
    len2 := Length_Of(psols);
    if len1 = len2 then
      put_line(file,"No paths lost in Start of the Cascade.  OK");
    else
      put(file,len1-len2,1);
      put(file," paths lost in Start of the Cascade.  BUG!!!");
      Scan_Failed_Paths(file,backsols,stapla,tarpla);
    end if;
    Clear(backsols);
  end Start_Intrinsic_Cascade;

  procedure Intrinsic_Diagonal_Homotopy
              ( file : in file_type; report : in boolean; a,b : in natural32;
                p1,p2 : in Poly_Sys; sols1,sols2 : in Solution_List;
                plane1,plane2 : in Matrix;
                sols : out Solution_List; plane : out Matrix ) is

    n : constant integer32 := integer32(Number_of_Unknowns(p1(p1'first)));
    n2 : constant integer32 := 2*n;
    apb : constant integer32 := integer32(a+b);
    m : constant integer32 := n2-apb;
    h0 : constant integer32
       := Minimal_Intersection_Dimension(n,integer32(a),integer32(b));
    AA : Matrix(1..apb,1..n);
    dA : Matrix(1..apb,1..n2);
    BB : Matrix(1..apb,1..n);
    CC : Matrix(1..n,1..n2);
    d : constant Vector := Random_Vector(1,n);
    offset : Vector(1..n2);
    target : Matrix(1..apb,0..n2);
    gentar : Matrix(1..n2,0..m);
    tarpla : Matrix(1..n2,0..m);
    stapla : Matrix(1..n2,0..m);
    targt1 : Matrix(1..apb,0..n2);
    gentr1,tarpl1 : Matrix(1..n2,0..m);
    f : Poly_Sys(1..m) := Product(n,n,p1,p2);
    ef : Eval_Poly_Sys(f'range) := Create(f);
    sjm : Jaco_Mat(f'range,1..n2) := Create(f);
    sjf : Eval_Jaco_Mat(sjm'range(1),sjm'range(2)) := Create(sjm);
    psols,esols : Solution_List;
    timer : Timing_Widget;
    tolsing : constant double_float := 1.0E-8;

  begin
   -- put_line("Intrinsic_Diagonal_Homotopy sets up planes:");
   -- put_line("  1. random orthogonal matrix for AA ...");
    AA := Standard_Random_Matrices.Random_Orthogonal_Matrix
            (natural32(apb),natural32(n));
    dA := Double_to_Diagonal(AA);
   -- put_line("  2. random orthogonal matrix for BB ...");
    BB := Standard_Random_Matrices.Random_Orthogonal_Matrix
            (natural32(apb),natural32(n));
   -- put_line("  3. random orthogonal matrix for CC ...");
   -- put("  n = "); put(n,1); put("  n2 = "); put(n2,1); new_line;
    CC := Standard_Random_Matrices.Random_Matrix(natural32(n),natural32(n2));
   -- put_line("computing the offset ...");
    offset := Compute_Offset(CC,d);
   -- put_line("computing target space ...");
    target := Target_Space(n,n2,apb,integer32(b),dA,BB,CC,d);
   -- put_line("computing generators of the target ...");
    gentar := Generators(target);
   -- put_line("orthogonalizing target plane ...");
    tarpla := Orthogonalize(gentar);
   -- put_line("computing start plane ...");
    stapla := Start_Space(plane1,plane2);
    tstart(timer);
    Shift_Offset(tarpla,offset);
    targt1 := Target_Space(n,n2,apb,integer32(b-1),dA,BB,CC,d);
    gentr1 := Generators(targt1);
    tarpl1 := Orthogonalize(gentr1);
    Shift_Offset(tarpl1,offset);
    put_line(file,"Intersecting planes...");
   -- Intersect(file,target,targt1,tarpla,tarpl1);
    Intersect(target,targt1,tarpla,tarpl1);
   -- put_line("making the product of the solution lists ...");
    psols := Product(sols1,sols2);
    Start_Intrinsic_Cascade(file,report,ef,sjf,stapla,tarpla,psols);
    new_line(file);
    put(file,"Starting the Cascade at Dimension ");
    put(file,b-1,1); put_line(file,".");
    new_line(file);
    put_line(file,"Refining the roots at the start...");
    Refine_Roots(file,ef,sjf,tarpla,psols);
    Set_Continuation_Parameter(psols,Create(0.0));
    put(file,"Tracking "); put(file,Length_Of(psols),1);
    put_line(file," paths to the top dimensional intersection component ...");
   -- Intrinsic_Track_Paths(file,report,ef,sjf,tarpla,tarpl1,psols);
    Recentered_Path_Tracking(file,report,ef,sjf,tarpla,tarpl1,psols);
    new_line(file);
    put_line(file,"Refining the roots at the target...");
    Refine_Roots(file,ef,sjf,tarpl1,psols);
    Singular_Filter(file,psols,tolsing);
    if not Is_Null(psols) then
      target := targt1;
      tarpla := tarpl1;
     -- running the cascade down to h0
      for h in reverse h0..integer32(b)-2 loop
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
       -- Intrinsic_Track_Paths(file,report,ef,sjf,tarpla,tarpl1,psols);
        Recentered_Path_Tracking(file,report,ef,sjf,tarpla,tarpl1,psols);
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
              ( file : in file_type; nefA,nefB,n,a,b : in integer32;
                sols1,sols2 : in Solution_List; plane1,plane2 : in Matrix;
                sols : out Solution_List; plane : out Matrix ) is

    n2 : constant integer32 := 2*n;
    apb : constant integer32 := a+b;
    m : constant integer32 := n2-apb;
    h0 : constant integer32 := Minimal_Intersection_Dimension(n,a,b);
    AA : constant Matrix(1..apb,1..n)
       := Standard_Random_Matrices.Random_Orthogonal_Matrix
            (natural32(apb),natural32(n));
    dA : constant Matrix(1..apb,1..n2) := Double_to_Diagonal(AA);
    BB : constant Matrix(1..apb,1..n)
       := Standard_Random_Matrices.Random_Orthogonal_Matrix
            (natural32(apb),natural32(n));
    CC : constant Matrix(1..n,1..n2)
       := Standard_Random_Matrices.Random_Matrix(natural32(n),natural32(n2));
    d : constant Vector := Random_Vector(1,n);
    offset : constant Vector(1..n2) := Compute_Offset(CC,d);
    target : Matrix(1..apb,0..n2) := Target_Space(n,n2,apb,b,dA,BB,CC,d);
    gentar : Matrix(1..n2,0..m) := Generators(target);
    tarpla : Matrix(1..n2,0..m) := Orthogonalize(gentar);
    stapla : constant Matrix(1..n2,0..m) := Start_Space(plane1,plane2);
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

end Intrinsic_Diagonal_Continuation;
