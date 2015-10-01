with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with QuadDobl_Random_Vectors;            use QuadDobl_Random_Vectors;
with QuadDobl_Complex_Vector_Norms;      use QuadDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Matrices;          use QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Linear_Solvers;    use QuadDobl_Complex_Linear_Solvers;
with Standard_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with QuadDobl_Solution_Diagnostics;      use QuadDobl_Solution_Diagnostics;
with QuadDobl_Condition_Tables;
with QuadDobl_Condition_Report;

package body QuadDobl_Root_Refiners is

-- ROOT ACCOUNTING :

  procedure Write_Info
              ( file : in file_type; zero : in Solution;
                initres : in quad_double;
                i,numb,nbdef : in natural32;
                fail,infty : in boolean ) is
  begin
    put(file,"solution "); put(file,i,1); put(file," : ");
    put(file,"   start residual : "); put(file,initres,3);
    if nbdef = 0
     then put(file,"   #iterations : "); put(file,numb,1);
     else put(file,"   #deflations : "); put(file,nbdef,1);
    end if;
    if infty then
      put_line(file,"   at infinity");
    elsif fail then
      put_line(file,"   failure");
    else
      put_line(file,"   success");
    end if;
    put(file,zero);
  end Write_Info;

  procedure Root_Accounting
              ( file : in file_type;
                h1,h2 : in QuadDobl_Complex_Vectors.Vector;
                pl : in out Point_List; ls : in Link_to_Solution;
                nb : in natural32; sa : in out Solution_Array;
                fail,infty,deflate : in boolean;
                tolsing,tolclus : in double_float; nbfail,nbinfty,
                nbreal,nbcomp,nbreg,nbsing,nbclus : in out natural32 ) is
  begin
    if infty then
      put_line(file," at infinity ==");
      nbinfty := nbinfty + 1;
    elsif fail then
      put_line(file," no solution ==");
      nbfail := nbfail + 1;
      ls.m := 0;
    else
      if Is_Real(ls.all,1.0E-13)
       then put(file," real ");    nbreal := nbreal + 1;
       else put(file," complex "); nbcomp := nbcomp + 1;
      end if;
      if deflate or sa(integer32(nb)).rco < tolsing then
        if ls.m <= 1 then  -- do not change input multiplicity field
          declare       -- to determine multiplicity, count clusters
            m : natural32; -- := Multiplicity(ls.all,sa,tolclus);
          begin
            QuadDobl_Condition_Report.Multiplicity
              (ls.all,nb,sa,tolclus,h1,h2,pl,m);
            if ((m = 1) and (not deflate))
             then m := 0;
            end if;
            ls.m := integer32(m);
          end;
        end if;
        if deflate then
          if ls.m = 1
           then put_line(file,"single ==");
           else put_line(file,"multiple ==");
                nbsing := nbsing + 1;
          end if;
          nbreg := nbreg + 1;
        else
          put_line(file,"singular ==");
          nbsing := nbsing + 1;
        end if;
      elsif sa(integer32(nb)).rco > tolsing then
        declare
         -- nb2 : constant natural32 := Is_Clustered(ls.all,nb,sa,tolclus);
          nb2 : natural32;
        begin
          QuadDobl_Condition_Report.Is_Clustered
            (ls.all,nb,sa,tolclus,h1,h2,pl,nb2);
          if nb2 = nb then
            put_line(file,"regular ==");
            nbreg := nbreg + 1;
            -- ls.m := 1;
          else
            put(file,"clustered : ");
            put(file,nb2,1);
            put_line(file," ==");
            nbclus := nbclus + 1;
            ls.m := -integer32(nb2);
            sa(integer32(nb2)).m := -integer32(nb);
          end if;
        end;	   
      end if;  
    end if;
  end Root_Accounting;

  procedure Write_Type
              ( file : in file_type; 
                h1,h2 : in QuadDobl_Complex_Vectors.Vector;
                pl : in out Point_List; ls : in Link_to_Solution;
                nb : in natural32; sa : in out Solution_Array;
                fail,infty,deflate : in boolean;
                tolsing,tolclus : in double_float; nbfail,nbinfty,
                nbreal,nbcomp,nbreg,nbsing,nbclus : in out natural32 ) is

    nb2 : natural32;

  begin
    if infty then
      put_line(file," at infinity ==");
      nbinfty := nbinfty + 1;
    elsif fail then
      put_line(file," no solution ==");
      nbfail := nbfail + 1;
      ls.m := 0;
    else
      if Is_Real(ls.all,1.0E-13)
       then put(file," real ");    nbreal := nbreal + 1;
       else put(file," complex "); nbcomp := nbcomp + 1;
      end if;
      if deflate or sa(integer32(nb)).rco < tolsing then
        if deflate then
          if ls.m = 1 then
            put_line(file,"single ==");
          else
            nb2 := Is_Clustered(ls.all,nb,sa,tolclus);
           -- QuadDobl_Condition_Report.Is_Clustered
           --   (ls.all,nb,sa,tolclus,h1,h2,pl,nb2);
            if nb2 = nb then
              put_line(file,"multiple ==");
            else
              put(file,"multiple: ");
              put(file,nb2,1); put_line(file," ==");
              nbsing := nbsing + 1;
            end if;
          end if;
          nbreg := nbreg + 1;
        else
          put_line(file,"singular ==");
          nbsing := nbsing + 1;
        end if;
      elsif sa(integer32(nb)).rco > tolsing then
        if ls.m > 0 then
          put_line(file,"regular ==");
          nbreg := nbreg + 1;
        else
          put(file,"clustered : ");
          put(file,ls.m,1);
          put_line(file," ==");
          nbclus := nbclus + 1;
        end if;
      end if;  
    end if;
  end Write_Type;

  procedure Write_Type
              ( file : in file_type; ls : in Link_to_Solution;
                fail,infty : in boolean;
                tolsing : in double_float;
                nbfail,nbinfty : in out natural32;
                nbreal,nbcomp,nbreg,nbsing : in out natural32 ) is
  begin
    if infty then
      put_line(file," at infinity =="); nbinfty := nbinfty + 1;
    elsif fail then
      put_line(file," no solution =="); nbfail := nbfail + 1;
      ls.m := 0;
    else
      if Is_Real(ls.all,1.0E-13)
       then put(file," real ");    nbreal := nbreal + 1;
       else put(file," complex "); nbcomp := nbcomp + 1;
      end if;
      if ls.rco < tolsing
       then put_line(file,"singular =="); nbsing := nbsing + 1;
       else put_line(file,"regular ==");  nbreg := nbreg + 1;
      end if;  
    end if;
  end Write_Type;

  procedure Multiplicity
              ( h1,h2 : in QuadDobl_Complex_Vectors.Vector;
                pl : in out Point_List; ls : in Link_to_Solution;
                nb : in natural32; sa : in out Solution_Array;
                fail,infty,deflate : in boolean;
                tolsing,tolclus : in double_float ) is
  begin
    if infty then
      null;
    elsif fail then
      ls.m := 0;
    elsif sa(integer32(nb)).rco < tolsing or deflate then
      if ls.m <= 1 then -- do not change input multiplicity field
        declare   -- to determine multiplicity count clustered solutions
         -- m : constant natural32 := Multiplicity(ls.all,sa,tolclus);
          m : natural32;
        begin
          QuadDobl_Condition_Report.Multiplicity
            (ls.all,nb,sa,tolclus,h1,h2,pl,m);
          if ((m = 1) and (not deflate))
           then ls.m := 0;
           else ls.m := integer32(m);
          end if;
        end;
      end if;
    else  -- sa(nb).rco > tolsing, check for clustering
      declare
        nb2 : constant natural32 := Is_Clustered(ls.all,nb,sa,tolclus);
      begin
        if nb2 /= nb
         then ls.m := -integer32(nb2); sa(integer32(nb2)).m := -integer32(nb);
        end if;
      end;	   
    end if;
  end Multiplicity;

  procedure Multiplicity
              ( h1,h2 : in QuadDobl_Complex_Vectors.Vector;
                pl : in out Point_List; ls : in Link_to_Solution;
                nb : in natural32; sols : in out Solution_List;
                fail,infty,deflate : in boolean;
                tolsing,tolclus : in double_float ) is
  begin
    if infty then
      null;
    elsif fail then
      ls.m := 0;
    elsif ls.rco < tolsing or deflate then
      if ls.m <= 1 then -- do not change input multiplicity field
        declare   -- to determine multiplicity count clustered solutions
         -- m : constant natural32 := Multiplicity(ls.all,sols,tolclus);
          m : natural32;
        begin
          QuadDobl_Condition_Report.Multiplicity
            (ls.all,nb,sols,tolclus,h1,h2,pl,m);
          if ((m = 1) and (not deflate))
           then ls.m := 0;
           else ls.m := integer32(m);
          end if;
        end;
      end if;
    else  -- ls.rco > tolsing, check for clustering
      declare
       -- nb2 : constant natural32 := Is_Clustered(ls.all,nb,sols,tolclus);
        nb2 : natural32;
      begin
        QuadDobl_Condition_Report.Is_Clustered
          (ls.all,nb,sols,tolclus,h1,h2,pl,nb2);
        if nb2 /= nb then
          ls.m := -integer32(nb2);
          Change_Multiplicity(sols,nb2,-integer32(nb));
        end if;
      end;	   
    end if;
  end Multiplicity;

  procedure Write_Global_Info
              ( file : in file_type; tot,nbfail,nbinfty,
                nbreal,nbcomp,nbreg,nbsing,nbclus : in natural32 ) is
  begin
    Standard_Complex_Solutions_io.put_bar(file);
    put(file,"A list of "); put(file,tot,1);
    put_line(file," solutions has been refined :");
    put(file,"Number of regular solutions     : "); put(file,nbreg,1);
    put_line(file,".");
    put(file,"Number of singular solutions    : "); put(file,nbsing,1);
    put_line(file,".");
    put(file,"Number of real solutions        : "); put(file,nbreal,1);
    put_line(file,".");
    put(file,"Number of complex solutions     : "); put(file,nbcomp,1);
    put_line(file,".");
    put(file,"Number of clustered solutions   : "); put(file,nbclus,1);
    put_line(file,".");
    put(file,"Number of solutions at infinity : "); put(file,nbinfty,1);
    put_line(file,".");
    put(file,"Number of failures              : "); put(file,nbfail,1);
    put_line(file,".");
    Standard_Complex_Solutions_io.put_bar(file);
  end Write_Global_Info;

-- ONE NEWTON STEP :

  procedure Write_Diagnostics
              ( file : in file_type; step : natural32;
                err,rco,res : in quad_double ) is
  begin
    put(file,"Step "); put(file,step,4); put(file," : ");
    put(file," |errxa| : "); put(file,err,3);
    put(file," est rco : "); put(file,rco,3);
    put(file," |errfa| : "); put(file,res,3); new_line(file);
  end Write_Diagnostics;

  procedure QuadDobl_Newton_Step
              ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys; 
                jf : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out quad_double ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;

    y : QuadDobl_Complex_Vectors.Vector(f'range) := eval(f,x);
    A : Matrix(f'range,f'range) := eval(jf,x);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2));
    info : integer32;
    Anorm : constant quad_double := Norm1(A);

  begin
    QuadDobl_Complex_Vectors.Min(y);
    lufac(A,A'last(1),ipvt,info);
    estco(A,A'last(1),ipvt,Anorm,rco);
    lusolve(A,A'last(1),ipvt,y);
    QuadDobl_Complex_Vectors.Add(x,y);
    err := Max_Norm(y);
    y := eval(f,x);
    res := Max_Norm(y);
  end QuadDobl_Newton_Step;

  procedure QuadDobl_Newton_Step
              ( f : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys; 
                jf : in QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out quad_double ) is

    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;

    y : QuadDobl_Complex_Vectors.Vector(f'range) := eval(f,x);
    A : Matrix(f'range,f'range) := eval(jf,x);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2));
    info : integer32;
    Anorm : constant quad_double := Norm1(A);

  begin
    QuadDobl_Complex_Vectors.Min(y);
    lufac(A,A'last(1),ipvt,info);
    estco(A,A'last(1),ipvt,Anorm,rco);
    lusolve(A,A'last(1),ipvt,y);
    QuadDobl_Complex_Vectors.Add(x,y);
    err := Max_Norm(y);
    y := eval(f,x);
    res := Max_Norm(y);
  end QuadDobl_Newton_Step;

  procedure QuadDobl_Newton_Step
              ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys; 
                jf : in QuadDobl_Jacobian_Circuits.Circuit;
                x : in out QuadDobl_Complex_Vectors.Vector;
                wrk : in out QuadDobl_Complex_VecVecs.VecVec;
                err,rco,res : out quad_double ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Jacobian_Circuits;

    y : QuadDobl_Complex_Vectors.Vector(f'range);
    A : Matrix(f'range,f'range);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2));
    info : integer32;
    Anorm : quad_double;

  begin
    EvalDiff(jf,x,wrk,y,A);
    QuadDobl_Complex_Vectors.Min(y);
    lufac(A,A'last(1),ipvt,info);
    Anorm := Norm1(A);
    estco(A,A'last(1),ipvt,Anorm,rco);
    lusolve(A,A'last(1),ipvt,y);
    QuadDobl_Complex_Vectors.Add(x,y);
    err := Max_Norm(y);
    y := eval(f,x);
    res := Max_Norm(y);
  end QuadDobl_Newton_Step;

  procedure Silent_Newton
              ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in double_float; numit : in out natural32;
                max : in natural32; fail : out boolean ) is
  begin
    fail := true;
    while numit < max loop
      numit := numit + 1;
      QuadDobl_Newton_Step(f,jf,x.v,x.err,x.rco,x.res);
      if (x.err < epsxa) or (x.res < epsfa)
       then fail := false; exit;
      end if;
    end loop;
  end Silent_Newton;

  procedure Reporting_Newton
              ( file : in file_type;
                f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in double_float; numit : in out natural32;
                max : in natural32; fail : out boolean ) is
  begin
    fail := true;
    while numit < max loop
      numit := numit + 1;
      QuadDobl_Newton_Step(f,jf,x.v,x.err,x.rco,x.res);
      Write_Diagnostics(file,numit,x.err,x.rco,x.res);
      if (x.err < epsxa) or (x.res < epsfa)
       then fail := false; exit;
      end if;
    end loop;
  end Reporting_Newton;

  procedure Silent_Newton
              ( f : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in  QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in double_float; numit : in out natural32;
                max : in natural32; fail : out boolean ) is
  begin
    fail := true;
    while numit < max loop
      numit := numit + 1;
      QuadDobl_Newton_Step(f,jf,x.v,x.err,x.rco,x.res);
      if (x.err < epsxa) or (x.res < epsfa)
       then fail := false; exit;
      end if;
    end loop;
  end Silent_Newton;

  procedure Reporting_Newton
              ( file : in file_type;
                f : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in  QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in double_float; numit : in out natural32;
                max : in natural32; fail : out boolean ) is
  begin
    fail := true;
    while numit < max loop
      numit := numit + 1;
      QuadDobl_Newton_Step(f,jf,x.v,x.err,x.rco,x.res);
      Write_Diagnostics(file,numit,x.err,x.rco,x.res);
      if (x.err < epsxa) or (x.res < epsfa)
       then fail := false; exit;
      end if;
    end loop;
  end Reporting_Newton;

  procedure Silent_Newton
              ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  QuadDobl_Jacobian_Circuits.Circuit;
                x : in out QuadDobl_Complex_Solutions.Solution;
                wrk : in out QuadDobl_Complex_VecVecs.VecVec;
                epsxa,epsfa : in double_float; numit : in out natural32;
                max : in natural32; fail : out boolean ) is
  begin
    fail := true;
    while numit < max loop
      numit := numit + 1;
      QuadDobl_Newton_Step(f,jf,x.v,wrk,x.err,x.rco,x.res);
      if (x.err < epsxa) or (x.res < epsfa)
       then fail := false; exit;
      end if;
    end loop;
  end Silent_Newton;

  procedure QuadDobl_Root_Refiner
              ( f : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                s : in QuadDobl_Complex_Solutions.Link_to_Solution ) is
  begin
    for i in 1..5 loop
      QuadDobl_Newton_Step(f,jf,s.v,s.err,s.rco,s.res);
     -- put("err : "); put(s.err,3);
     -- put(" = rco : "); put(s.rco,3);
     -- put(" = res : "); put(s.res,3); new_line;
    end loop;
  end QuadDobl_Root_Refiner;

  procedure QuadDobl_Root_Refiner
              ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in QuadDobl_Jacobian_Circuits.Circuit;
                s : in QuadDobl_Complex_Solutions.Link_to_Solution;
                wrk : in out QuadDobl_Complex_VecVecs.VecVec ) is
  begin
    for i in 1..5 loop
      QuadDobl_Newton_Step(f,jf,s.v,wrk,s.err,s.rco,s.res);
     -- put("err : "); put(s.err,3);
     -- put(" = rco : "); put(s.rco,3);
     -- put(" = res : "); put(s.res,3); new_line;
    end loop;
  end QuadDobl_Root_Refiner;

  procedure QuadDobl_Root_Refiner
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                s : in out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;
    use QuadDobl_Complex_Solutions;

    f : Eval_Laur_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,p'range) := Create(p);
    jf : Eval_Jaco_Mat(p'range,p'range) := Create(jm);
    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      QuadDobl_Root_Refiner(f,jf,ls);
      tmp := Tail_Of(tmp);
    end loop;
    Clear(f); Clear(jm); Clear(jf);
  end QuadDobl_Root_Refiner;

  procedure QuadDobl_Root_Refiner
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Jacobian_Circuits;
    use QuadDobl_Complex_Solutions;

    f : Eval_Poly_Sys(p'range) := Create(p);
    jf : Circuit := Create(p);
    nm : constant integer32 := integer32(Number_of_Monomials(jf));
    wrk : QuadDobl_Complex_VecVecs.VecVec(1..nm) := WorkSpace(jf);
    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      QuadDobl_Root_Refiner(f,jf,ls,wrk);
      tmp := Tail_Of(tmp);
    end loop;
    QuadDobl_Complex_Poly_SysFun.Clear(f);
    QuadDobl_Jacobian_Circuits.Clear(jf);
    QuadDobl_Complex_VecVecs.Clear(wrk);
  end QuadDobl_Root_Refiner;

-- THE MAIN ROOT REFINERS :

  procedure Silent_Root_Refiner
               ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 s : in out QuadDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32 ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Complex_Solutions;
    
    n : constant integer32 := p'last;
    f : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,p'range) := Create(p);
    jf : Eval_Jaco_Mat(p'range,p'range) := Create(jm);
    sa : Solution_Array(1..integer32(Length_Of(s))) := Create(s);
    nb : natural32;
    fail,infty : boolean;
    h1 : constant QuadDobl_Complex_Vectors.Vector
       := QuadDobl_Random_Vectors.Random_Vector(1,n);
    h2 : constant QuadDobl_Complex_Vectors.Vector
       := QuadDobl_Random_Vectors.Random_Vector(1,n);
    pl : Point_List;

  begin
    for i in sa'range loop
      nb := 0;
      sa(i).res := Sum_Norm(Eval(f,sa(i).v));
      infty := At_Infinity(sa(i).all,false,1.0E+8);
      if not infty and sa(i).res < 0.1 and sa(i).err < 0.1 then
        Silent_Newton(f,jf,sa(i).all,epsxa,epsfa,nb,max,fail);
        Multiplicity(h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,
                     infty,false,tolsing,epsxa);
        numit := numit + nb;
      else
        fail := true;
      end if;
    end loop;
    Deep_Clear(s); s := Create(sa); Clear(sa);
    Clear(pl); Clear(f); Clear(jm); Clear(jf);
  end Silent_Root_Refiner;

  procedure Silent_Root_Refiner
               ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 s,refs : in out QuadDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32 ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Complex_Solutions;
    
    n : constant integer32 := p'last;
    f : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,p'range) := Create(p);
    jf : Eval_Jaco_Mat(p'range,p'range) := Create(jm);
    sa : Solution_Array(1..integer32(Length_Of(s))) := Create(s);
    refs_last : Solution_List;
    nb : natural32;
    fail,infty : boolean;
    h1 : constant QuadDobl_Complex_Vectors.Vector
       := QuadDobl_Random_Vectors.Random_Vector(1,n);
    h2 : constant QuadDobl_Complex_Vectors.Vector
       := QuadDobl_Random_Vectors.Random_Vector(1,n);
    pl : Point_List;

  begin
    for i in sa'range loop
      nb := 0;
      sa(i).res := Sum_Norm(Eval(f,sa(i).v));
      infty := At_Infinity(sa(i).all,false,1.0E+8);
      if not infty and sa(i).res < 0.1 and sa(i).err < 0.1 then
        Silent_Newton(f,jf,sa(i).all,epsxa,epsfa,nb,max,fail);
        Multiplicity(h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,
                     infty,false,tolsing,epsxa);
        numit := numit + nb;
        if not fail
         then Append(refs,refs_last,sa(i).all);
        end if;
      else
        fail := true;
      end if;
    end loop;
    Deep_Clear(s); s := Create(sa); Clear(sa);
    Clear(pl); Clear(f); Clear(jm); Clear(jf);
  end Silent_Root_Refiner;

  procedure Silent_Root_Refiner
               ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 s : in out QuadDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32 ) is

    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;
    use QuadDobl_Complex_Solutions;
    
    n : constant integer32 := p'last;
    f : Eval_Laur_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,p'range) := Create(p);
    jf : Eval_Jaco_Mat(p'range,p'range) := Create(jm);
    sa : Solution_Array(1..integer32(Length_Of(s))) := Create(s);
    nb : natural32;
    fail,infty : boolean;
    h1 : constant QuadDobl_Complex_Vectors.Vector
       := QuadDobl_Random_Vectors.Random_Vector(1,n);
    h2 : constant QuadDobl_Complex_Vectors.Vector
       := QuadDobl_Random_Vectors.Random_Vector(1,n);
    pl : Point_List;

  begin
    for i in sa'range loop
      nb := 0;
      sa(i).res := Sum_Norm(Eval(f,sa(i).v));
      infty := At_Infinity(sa(i).all,false,1.0E+8);
      if not infty and sa(i).res < 0.1 and sa(i).err < 0.1 then
        Silent_Newton(f,jf,sa(i).all,epsxa,epsfa,nb,max,fail);
        Multiplicity(h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,
                     infty,false,tolsing,epsxa);
        numit := numit + nb;
      else
        fail := true;
      end if;
    end loop;
    Deep_Clear(s); s := Create(sa); Clear(sa);
    Clear(pl); Clear(f); Clear(jm); Clear(jf);
  end Silent_Root_Refiner;

  procedure Silent_Root_Refiner
               ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 s,refs : in out QuadDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32 ) is

    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;
    use QuadDobl_Complex_Solutions;
    
    n : constant integer32 := p'last;
    f : Eval_Laur_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,p'range) := Create(p);
    jf : Eval_Jaco_Mat(p'range,p'range) := Create(jm);
    sa : Solution_Array(1..integer32(Length_Of(s))) := Create(s);
    refs_last : Solution_List;
    nb : natural32;
    fail,infty : boolean;
    h1 : constant QuadDobl_Complex_Vectors.Vector
       := QuadDobl_Random_Vectors.Random_Vector(1,n);
    h2 : constant QuadDobl_Complex_Vectors.Vector
       := QuadDobl_Random_Vectors.Random_Vector(1,n);
    pl : Point_List;

  begin
    for i in sa'range loop
      nb := 0;
      sa(i).res := Sum_Norm(Eval(f,sa(i).v));
      infty := At_Infinity(sa(i).all,false,1.0E+8);
      if not infty and sa(i).res < 0.1 and sa(i).err < 0.1 then
        Silent_Newton(f,jf,sa(i).all,epsxa,epsfa,nb,max,fail);
        Multiplicity(h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,
                     infty,false,tolsing,epsxa);
        numit := numit + nb;
        if not fail
         then Append(refs,refs_last,sa(i).all);
        end if;
      else
        fail := true;
      end if;
    end loop;
    Deep_Clear(s); s := Create(sa); Clear(sa);
    Clear(pl); Clear(f); Clear(jm); Clear(jf);
  end Silent_Root_Refiner;

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 s : in out QuadDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 wout : in boolean ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Complex_Solutions;
    
    n : constant integer32 := p'last;
    f : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,p'range) := Create(p);
    jf : Eval_Jaco_Mat(p'range,p'range) := Create(jm);
    sa : Solution_Array(1..integer32(Length_Of(s))) := Create(s);
    nb : natural32;
    nbfail,nbinfty,nbreg,nbsing,nbclus,nbreal,nbcomp : natural32 := 0;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..60)
                      := QuadDobl_Condition_Tables.Create(60); 
    fail,infty : boolean;
    h1 : constant QuadDobl_Complex_Vectors.Vector
       := QuadDobl_Random_Vectors.Random_Vector(1,n);
    h2 : constant QuadDobl_Complex_Vectors.Vector
       := QuadDobl_Random_Vectors.Random_Vector(1,n);
    pl : Point_List;
    initres : quad_double;

  begin
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,sa'last,1); put(file," "); put(file,n,1); new_line(file);
    Standard_Complex_Solutions_io.put_bar(file);
    for i in sa'range loop
      nb := 0;
      sa(i).res := Sum_Norm(Eval(f,sa(i).v));
      initres := sa(i).res;
      infty := At_Infinity(sa(i).all,false,1.0E+8);
      if not infty and sa(i).res < 0.1 and sa(i).err < 0.1 then
        if wout
         then Reporting_Newton(file,f,jf,sa(i).all,epsxa,epsfa,nb,max,fail);
         else Silent_Newton(f,jf,sa(i).all,epsxa,epsfa,nb,max,fail);
        end if;
        numit := numit + nb;
      else
        fail := true;
      end if;
      Multiplicity(h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,
                   infty,false,tolsing,epsxa);
      Write_Info(file,sa(i).all,initres,natural32(i),nb,0,fail,infty);
      Write_Type
        (file,h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,infty,false,
         tolsing,epsxa,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
      QuadDobl_Condition_Tables.Update_Corrector(t_err,sa(i).all);
      QuadDobl_Condition_Tables.Update_Condition(t_rco,sa(i).all);
      QuadDobl_Condition_Tables.Update_Residuals(t_res,sa(i).all);
    end loop;
    Write_Global_Info
      (file,natural32(sa'last),
       nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
    QuadDobl_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    Deep_Clear(s); s := Create(sa); Clear(sa);
    Clear(pl); Clear(f); Clear(jm); Clear(jf);
  end Reporting_Root_Refiner;

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 s,refs : in out QuadDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 wout : in boolean ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Complex_Solutions;
    
    n : constant integer32 := p'last;
    f : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,p'range) := Create(p);
    jf : Eval_Jaco_Mat(p'range,p'range) := Create(jm);
    sa : Solution_Array(1..integer32(Length_Of(s))) := Create(s);
    refs_last : Solution_List;
    nb : natural32;
    nbfail,nbinfty,nbreg,nbsing,nbclus,nbreal,nbcomp : natural32 := 0;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..60)
                      := QuadDobl_Condition_Tables.Create(60); 
    fail,infty : boolean;
    h1 : constant QuadDobl_Complex_Vectors.Vector
       := QuadDobl_Random_Vectors.Random_Vector(1,n);
    h2 : constant QuadDobl_Complex_Vectors.Vector
       := QuadDobl_Random_Vectors.Random_Vector(1,n);
    pl : Point_List;
    initres : quad_double;

  begin
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,sa'last,1); put(file," "); put(file,n,1); new_line(file);
    Standard_Complex_Solutions_io.put_bar(file);
    for i in sa'range loop
      nb := 0;
      sa(i).res := Sum_Norm(Eval(f,sa(i).v));
      initres := sa(i).res;
      infty := At_Infinity(sa(i).all,false,1.0E+8);
      if not infty and sa(i).res < 0.1 and sa(i).err < 0.1 then
        if wout
         then Reporting_Newton(file,f,jf,sa(i).all,epsxa,epsfa,nb,max,fail);
         else Silent_Newton(f,jf,sa(i).all,epsxa,epsfa,nb,max,fail);
        end if;
        numit := numit + nb;
      else
        fail := true;
      end if;
      Multiplicity(h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,
                   infty,false,tolsing,epsxa);
      Write_Info(file,sa(i).all,initres,natural32(i),nb,0,fail,infty);
      Write_Type
        (file,h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,infty,false,
         tolsing,epsxa,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
      QuadDobl_Condition_Tables.Update_Corrector(t_err,sa(i).all);
      QuadDobl_Condition_Tables.Update_Condition(t_rco,sa(i).all);
      QuadDobl_Condition_Tables.Update_Residuals(t_res,sa(i).all);
      if not fail
       then Append(refs,refs_last,sa(i).all);
      end if;
    end loop;
    Write_Global_Info
      (file,natural32(sa'last),
       nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
    QuadDobl_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    Deep_Clear(s); s := Create(sa); Clear(sa);
    Clear(pl); Clear(f); Clear(jm); Clear(jf);
  end Reporting_Root_Refiner;

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 s : in out QuadDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 wout : in boolean ) is

    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;
    use QuadDobl_Complex_Solutions;
    
    n : constant integer32 := p'last;
    f : Eval_Laur_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,p'range) := Create(p);
    jf : Eval_Jaco_Mat(p'range,p'range) := Create(jm);
    sa : Solution_Array(1..integer32(Length_Of(s))) := Create(s);
    nb : natural32;
    nbfail,nbinfty,nbreg,nbsing,nbclus,nbreal,nbcomp : natural32 := 0;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..60)
                      := QuadDobl_Condition_Tables.Create(60); 
    fail,infty : boolean;
    h1 : constant QuadDobl_Complex_Vectors.Vector
       := QuadDobl_Random_Vectors.Random_Vector(1,n);
    h2 : constant QuadDobl_Complex_Vectors.Vector
       := QuadDobl_Random_Vectors.Random_Vector(1,n);
    pl : Point_List;
    initres : quad_double;

  begin
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,sa'last,1); put(file," "); put(file,n,1); new_line(file);
    Standard_Complex_Solutions_io.put_bar(file);
    for i in sa'range loop
      nb := 0;
      sa(i).res := Sum_Norm(Eval(f,sa(i).v));
      initres := sa(i).res;
      infty := At_Infinity(sa(i).all,false,1.0E+8);
      if not infty and sa(i).res < 0.1 and sa(i).err < 0.1 then
        if wout
         then Reporting_Newton(file,f,jf,sa(i).all,epsxa,epsfa,nb,max,fail);
         else Silent_Newton(f,jf,sa(i).all,epsxa,epsfa,nb,max,fail);
        end if;
        numit := numit + nb;
      else
        fail := true;
      end if;
      Multiplicity(h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,
                   infty,false,tolsing,epsxa);
      Write_Info(file,sa(i).all,initres,natural32(i),nb,0,fail,infty);
      Write_Type
        (file,h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,infty,false,
         tolsing,epsxa,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
      QuadDobl_Condition_Tables.Update_Corrector(t_err,sa(i).all);
      QuadDobl_Condition_Tables.Update_Condition(t_rco,sa(i).all);
      QuadDobl_Condition_Tables.Update_Residuals(t_res,sa(i).all);
    end loop;
    Write_Global_Info
      (file,natural32(sa'last),
       nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
    QuadDobl_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    Deep_Clear(s); s := Create(sa); Clear(sa);
    Clear(pl); Clear(f); Clear(jm); Clear(jf);
  end Reporting_Root_Refiner;

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 s,refs : in out QuadDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 wout : in boolean ) is

    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;
    use QuadDobl_Complex_Solutions;
    
    n : constant integer32 := p'last;
    f : Eval_Laur_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,p'range) := Create(p);
    jf : Eval_Jaco_Mat(p'range,p'range) := Create(jm);
    sa : Solution_Array(1..integer32(Length_Of(s))) := Create(s);
    refs_last : Solution_List;
    nb : natural32;
    nbfail,nbinfty,nbreg,nbsing,nbclus,nbreal,nbcomp : natural32 := 0;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..60)
                      := QuadDobl_Condition_Tables.Create(60); 
    fail,infty : boolean;
    h1 : constant QuadDobl_Complex_Vectors.Vector
       := QuadDobl_Random_Vectors.Random_Vector(1,n);
    h2 : constant QuadDobl_Complex_Vectors.Vector
       := QuadDobl_Random_Vectors.Random_Vector(1,n);
    pl : Point_List;
    initres : quad_double;

  begin
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,sa'last,1); put(file," "); put(file,n,1); new_line(file);
    Standard_Complex_Solutions_io.put_bar(file);
    for i in sa'range loop
      nb := 0;
      sa(i).res := Sum_Norm(Eval(f,sa(i).v));
      initres := sa(i).res;
      infty := At_Infinity(sa(i).all,false,1.0E+8);
      if not infty and sa(i).res < 0.1 and sa(i).err < 0.1 then
        if wout
         then Reporting_Newton(file,f,jf,sa(i).all,epsxa,epsfa,nb,max,fail);
         else Silent_Newton(f,jf,sa(i).all,epsxa,epsfa,nb,max,fail);
        end if;
        numit := numit + nb;
      else
        fail := true;
      end if;
      Multiplicity(h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,
                   infty,false,tolsing,epsxa);
      Write_Info(file,sa(i).all,initres,natural32(i),nb,0,fail,infty);
      Write_Type
        (file,h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,infty,false,
         tolsing,epsxa,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
      QuadDobl_Condition_Tables.Update_Corrector(t_err,sa(i).all);
      QuadDobl_Condition_Tables.Update_Condition(t_rco,sa(i).all);
      QuadDobl_Condition_Tables.Update_Residuals(t_res,sa(i).all);
      if not fail
       then Append(refs,refs_last,sa(i).all);
      end if;
    end loop;
    Write_Global_Info
      (file,natural32(sa'last),
       nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
    QuadDobl_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    Deep_Clear(s); s := Create(sa); Clear(sa);
    Clear(pl); Clear(f); Clear(jm); Clear(jf);
  end Reporting_Root_Refiner;

end QuadDobl_Root_Refiners;
