with Timing_Package;                    use Timing_Package;
with Communications_with_User;          use Communications_with_User;
with File_Scanning,Numbers_io;          use File_Scanning;
with Characters_and_Numbers;            use Characters_and_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Quad_Double_Numbers_io;            use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with QuadDobl_Complex_Numbers;          use QuadDobl_Complex_Numbers;
with QuadDobl_Random_Numbers;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Complex_Vectors;
with Quad_Double_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vector_Norms;     use QuadDobl_Complex_Vector_Norms;
with Standard_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Polynomial_Convertors;    use QuadDobl_Polynomial_Convertors;
with QuadDobl_Homotopy;
with QuadDobl_Complex_Solutions_io;     use QuadDobl_Complex_Solutions_io;
with Continuation_Parameters;
with QuadDobl_Continuation_Data_io;     use QuadDobl_Continuation_Data_io;
with QuadDobl_Path_Trackers;            use QuadDobl_Path_Trackers;
with Main_Poly_Continuation;            use Main_Poly_Continuation;
with Lexicographic_Root_Enumeration;
with Total_Degree_Start_Systems;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Lists;
with Standard_Complex_Prod_Systems;
with Standard_Complex_Prod_Systems_io;
with Standard_Complex_Prod_Planes;
with Standard_Linear_Product_System;
with QuadDobl_Linear_Product_System;

package body Drivers_to_Track_QuadDobl_Paths is

  procedure Report_Kind ( kind : in natural32 ) is

  -- DESCRIPTION :
  --   Writes to screen a short summary of the kind of path 
  --   that was tracked, depending on the value of kind:
  --   if 0, then "diverged"; if 1, then "converged to regular";
  --   if 2, then "converged to singular".

  begin
    case kind is
      when 0 => put_line("  diverged");
      when 1 => put_line("  regular");
      when 2 => put_line("  singular");
      when others => put_line("  kind unknown");
    end case;
  end Report_Kind;

  procedure Track ( file : in file_type; report : in boolean;
                    p,q : in Poly_Sys; len,dim : in natural32;
                    f : access procedure ( s : in Solu_Info ) := null ) is

    gamma : constant Complex_Number := QuadDobl_Random_Numbers.Random1;
    timer : Timing_Widget;
    ls : Link_to_Solution;
    s : Solu_Info;
    kind,ind : natural32;
    continue : boolean;
    pp1 : constant Continuation_Parameters.Pred_Pars
        := Continuation_Parameters.Create_for_Path;
    pp2 : constant Continuation_Parameters.Pred_Pars
        := Continuation_Parameters.Create_End_Game;
    cp1 : constant Continuation_Parameters.Corr_Pars
        := Continuation_Parameters.Create_for_Path;
    cp2 : constant Continuation_Parameters.Corr_Pars
        := Continuation_Parameters.Create_End_Game;
    t1 : constant Complex_Number := Create(integer(1));
    tol : constant double_float := 1.0E-14;
    tol_zero : constant double_float := 1.0E-20;
    w : integer32 := 1;
    v : Quad_Double_Vectors.Link_to_Vector;
    e : quad_double := create(0.0);
   -- order : constant natural := Continuation_Parameters.endext_order;
    cnt,nbfail,nbregu,nbsing : natural32 := 0;

    procedure Track_along_Path is
      new Linear_Single_Normal_Reporting_Continue
            (Max_Norm,QuadDobl_Homotopy.Eval,
                      QuadDobl_Homotopy.Diff,QuadDobl_Homotopy.Diff);

    procedure Track_End_Game is
      new Linear_Single_Conditioned_Reporting_Continue
            (Max_Norm,QuadDobl_Homotopy.Eval,
                      QuadDobl_Homotopy.Diff,QuadDobl_Homotopy.Diff);

  begin
    QuadDobl_Homotopy.Create(p,q,2,gamma);
    tstart(timer);
    for i in 1..len loop
      Get_Next(dim,ls,ind,continue);
      exit when (ls = null);
      ls.t := Create(integer(0));
      s := Shallow_Create(ls);
      Track_along_Path(file,s,t1,tol,false,pp1,cp1,f=>f);
      Track_End_Game(file,s,t1,tol,false,0,w,v,e,pp2,cp2,f=>f);
      cnt := ind-1;
      Write_Next_Solution
        (file,cnt,s,tol_zero,tol_zero,nbfail,nbregu,nbsing,kind);
      text_io.flush(file);
      if report then Report_Kind(kind); end if;
      exit when not continue;
      Clear(ls);
    end loop;
    tstop(timer);
    put(file,"Number of path failures      : ");
    put(file,nbfail,1); new_line(file);
    put(file,"Number of regular solutions  : ");
    put(file,nbregu,1); new_line(file);
    put(file,"Number of singular solutions : ");
    put(file,nbsing,1); new_line(file);
    new_line(file);
    print_times(file,timer,"tracking " & Convert(integer32(len)) & " paths");
  end Track;

  procedure Cheater_Track
              ( outfile,solfile : in file_type; report : in boolean;
                p,q : in Poly_Sys; len,dim : in natural32;
                f : access procedure ( s : in Solu_Info ) := null ) is

    ndp : natural32;
    cnt : natural32 := 0;

    procedure Read_Next_Start_Solution
                ( n : in natural32; ls : out Link_to_Solution;
                  ind : out natural32; continue : out boolean ) is
    begin
      cnt := cnt + 1;
      if report
       then put(cnt,ndp);
      end if;
      Read_Next(solfile,n,ls);
      ind := cnt;
      continue := true;
    end Read_Next_Start_Solution;
    procedure Cheater_Homotopy_Track is new Track(Read_Next_Start_Solution);

  begin
    if report
     then ndp := Numbers_io.Number_of_Decimal_Places(len);
    end if;
    Cheater_Homotopy_Track(outfile,report,p,q,len,dim,f);
  end Cheater_Track;

  procedure Total_Degree_Track
              ( file : in file_type; report : in boolean;
                p,q : in Poly_Sys; len,dim : in natural32;
                f : access procedure ( s : in Solu_Info ) := null ) is

    d : constant Standard_Natural_Vectors.Vector(q'range)
      := Total_Degree_Start_Systems.Degrees(q);
    cp : constant Standard_Natural_Vectors.Vector(d'first..d'last-1)
       := Lexicographic_Root_Enumeration.Consecutive_Products(d);
    c : constant QuadDobl_Complex_Vectors.Vector(q'range)
      := Total_Degree_Start_Systems.Coefficients(q);
    cnt : natural32 := 0;
    ndp : natural32;

    procedure Compute_Next_Start_Solution
                ( n : in natural32; ls : out Link_to_Solution;
                  ind : out natural32; continue : out boolean ) is

      s : Standard_Natural_Vectors.Vector(d'range);
      r : QuadDobl_Complex_Vectors.Vector(d'range);

    begin
      cnt := cnt + 1;
      s := Lexicographic_Root_Enumeration.Root_Map(n,cnt,d,cp);
      r := Total_Degree_Start_Systems.Root(d,s,c);
      if report
       then put(cnt,ndp); put("  s = "); put(s); 
      end if;
      ls := new Solution'(Total_Degree_Start_Systems.Create(r));
      ind := cnt;
      continue := true;
    end Compute_Next_Start_Solution;
    procedure Total_Degree_Homotopy_Track is 
      new Track(Compute_Next_Start_Solution);

  begin
    if report
     then ndp := Numbers_io.Number_of_Decimal_Places(len);
    end if;
    Total_Degree_Homotopy_Track(file,report,p,q,len,dim,f);
  end Total_Degree_Track;

  procedure Get_Lex_Linear_Product_Start_Solution
              ( report : in boolean; n : in natural32;
                d,cp : in Standard_Natural_Vectors.Vector;
                ind : in natural32; ndp : in natural32;
                tol : in double_float;
                ls : out Link_to_Solution; fail : out boolean ) is

    s : Standard_Natural_Vectors.Vector(d'range);
    r : QuadDobl_Complex_Vectors.Vector(d'range);
    rcond : quad_double;

  begin
    s := Lexicographic_Root_Enumeration.Root_Map(n,ind,d,cp);
    QuadDobl_Linear_Product_System.Solve(s,tol,rcond,fail,r);
    if report then
      put(ind,ndp); put("  s = "); put(s);
      if fail then put(" fail "); else put(" okay "); end if;
      put(" rcond ="); put(rcond,3);
      if fail then put_line("  no path"); end if;
    end if;
    if not fail
     then ls := new Solution'(Total_Degree_Start_Systems.Create(r));
    end if;
  end Get_Lex_Linear_Product_Start_Solution;

  procedure Next_Lex_Linear_Product_Start_Solution
              ( report : in boolean; n : in natural32;
                d,cp : in Standard_Natural_Vectors.Vector;
                cnt : in out natural32; len,ndp : in natural32;
                tol : in double_float;
                ls : out Link_to_Solution; fail : out boolean ) is

    s : Standard_Natural_Vectors.Vector(d'range);
    r : QuadDobl_Complex_Vectors.Vector(d'range);
    rcond : quad_double;

  begin
    while cnt < len loop
      cnt := cnt + 1;
      s := Lexicographic_Root_Enumeration.Root_Map(n,cnt,d,cp);
      QuadDobl_Linear_Product_System.Solve(s,tol,rcond,fail,r);
      if report then
        put(cnt,ndp); put("  s ="); put(s);
        if fail then put(" fail "); else put(" okay "); end if;
        put(" rcond ="); put(rcond,3);
        if fail then put_line("  no path"); end if;
      end if;
      if not fail
       then ls := new Solution'(Total_Degree_Start_Systems.Create(r,rcond));
      end if;
      exit when not fail;
    end loop;
  end Next_Lex_Linear_Product_Start_Solution;

  procedure Get_Next_Linear_Product_Start_Solution
              ( report : in boolean; n : in natural32;
                s : in out Standard_Natural_Vectors.Vector;
                d : in Standard_Natural_Vectors.Vector;
                cnt : in out natural32; len,ndp : in natural32;
                tol : in double_float;
                ls : out Link_to_Solution; fail : out boolean ) is

    r : QuadDobl_Complex_Vectors.Vector(s'range);
    rcond : quad_double;

  begin
    if cnt = 0
     then QuadDobl_Linear_Product_System.Get_First(tol,s,fail);
     else QuadDobl_Linear_Product_System.Get_Next(tol,d,s,fail);
    end if;
    if not fail then
      cnt := cnt + 1;
      QuadDobl_Linear_Product_System.Solve(s,tol,rcond,fail,r);
      if report then
        put(cnt,ndp); put("  s ="); put(s);
        if fail then put(" fail "); else put(" okay "); end if;
        put(" rcond ="); put(rcond,3);
        if fail then put_line(" no path"); end if;
      end if;
      if not fail 
       then ls := new Solution'(Total_Degree_Start_Systems.Create(r,rcond));
      end if;
    end if;
  end Get_Next_Linear_Product_Start_Solution;

  procedure Get_Next_Linear_Product_Start_Solution
              ( report : in boolean; n : in natural32;
                cnt : in out natural32; len,ndp : in natural32;
                tol : in double_float;
                ls : out Link_to_Solution; fail : out boolean ) is

    r : QuadDobl_Complex_Vectors.Vector(1..integer32(n));
    s : Standard_Natural_Vectors.Vector(1..integer32(n));
    rcond : quad_double;

  begin
    if cnt = 0
     then QuadDobl_Linear_Product_System.Get_First(tol,s,fail);
     else QuadDobl_Linear_Product_System.Get_Next(tol,s,fail);
    end if;
    if not fail then
      cnt := cnt + 1;
      QuadDobl_Linear_Product_System.Solve(s,tol,rcond,fail,r);
      if report then
        put(cnt,ndp); put("  s ="); put(s);
        if fail then put(" fail "); else put(" okay "); end if;
        put(" rcond ="); put(rcond,3);
        if fail then put_line(" no path"); end if;
      end if;
      if not fail 
       then ls := new Solution'(Total_Degree_Start_Systems.Create(r,rcond));
      end if;
    end if;
  end Get_Next_Linear_Product_Start_Solution;

  procedure Linear_Product_Track
              ( file : in file_type; report : in boolean;
                p,q : in Poly_Sys; len,dim : in natural32;
                f : access procedure ( s : in Solu_Info ) := null ) is

    d : constant Standard_Natural_Vectors.Vector(q'range)
      := Total_Degree_Start_Systems.Degrees(q);
    s : Standard_Natural_Vectors.Vector(q'range) := (q'range => 0);
   -- cp : Standard_Natural_Vectors.Vector(d'first..d'last-1)
   --    := Lexicographic_Root_Enumeration.Consecutive_Products(d);
    cnt : natural32 := 0;
    tol : constant double_float := 1.0E-10;
    ndp : natural32;

    procedure Compute_Next_Start_Solution
                ( n : in natural32; ls : out Link_to_Solution;
                  ind : out natural32; continue : out boolean ) is

      fail : boolean;

    begin
     -- Next_Lex_Linear_Product_Start_Solution
     --   (report,n,d,cp,cnt,len,ndp,tol,ls,fail);
      Get_Next_Linear_Product_Start_Solution
        (report,n,s,d,cnt,len,ndp,tol,ls,fail);
      ind := cnt;
      if cnt > len
       then continue := false;
       else continue := true;
      end if;
    end Compute_Next_Start_Solution;
    procedure Linear_Product_Homotopy_Track is 
      new Track(Compute_Next_Start_Solution);

  begin
    if report
     then ndp := Numbers_io.Number_of_Decimal_Places(len);
    end if;
    Linear_Product_Homotopy_Track(file,report,p,q,len,dim,f);
  end Linear_Product_Track;

  procedure Read_Target_System
              ( file : in file_type; p : out Link_to_Poly_Sys ) is
  begin
    get(file,p);
  exception
    when others => put_line("Exception raised while reading target system.");
                   raise;
  end Read_Target_System;

  procedure Read_Start_System
              ( file : in file_type; p : out Link_to_Poly_Sys ) is
  begin
    get(file,p);
  exception
    when others => put_line("Exception raised while reading start system.");
                   raise;
  end Read_Start_System;

  procedure Hyperplane_Coefficients
              ( p : in Standard_Complex_Polynomials.Poly;
                fail : out boolean;
                h : out Standard_Complex_Vectors.Vector;
                qd_h : out QuadDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Returns in h the coefficients of the hyperplane p,
  --   fail becomes false if p is nonlinear.

    use Standard_Complex_Polynomials;

    procedure Scan_Term ( t : in Term; continue : out boolean ) is

      ind : constant integer32
          := Standard_Complex_Prod_Planes.Linear_Index(t.dg.all);

    begin
      if ind < 0 or ind > h'last then
        fail := true;
        continue := false;
        return;
      else
        declare
          re : constant double_float
             := Standard_Complex_Numbers.REAL_PART(t.cf);
          im : constant double_float
             := Standard_Complex_Numbers.IMAG_PART(t.cf);
          dd_re : constant double_double := create(re);
          dd_im : constant double_double := create(im);
          qd_re : constant quad_double := create(dd_re);
          qd_im : constant quad_double := create(dd_im);
        begin
          h(ind) := t.cf;
          qd_h(ind) := QuadDobl_Complex_Numbers.Create(qd_re,qd_im);
        end;
        continue := true;
      end if;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    for i in h'range loop
      h(i) := Standard_Complex_Numbers.Create(integer(0));
      qd_h(i) := Create(integer(0));
    end loop;
    fail := false;
    Scan_Terms(p);
  end Hyperplane_Coefficients;

  procedure Store ( p : in Standard_Complex_Prod_Systems.Prod_Sys;
                    fail : out boolean ) is

  -- DESCRIPTION :
  --   Stores the linear product start system defined by p.

    use Standard_Complex_Polynomials,Standard_Complex_Poly_Lists;
    q : Poly;
    tmp : Prod_Poly;
    h : Standard_Complex_Vectors.Vector(0..p'last);
    qd_h : QuadDobl_Complex_Vectors.Vector(0..p'last);

  begin
    Standard_Linear_Product_System.Init(natural32(p'last));
    QuadDobl_Linear_Product_System.Init(natural32(p'last));
    for i in p'range loop
      tmp := p(i);
      while not Is_Null(tmp) loop
        q := Head_Of(tmp);
        if Degree(q) > 1 then
          fail := true; return;
        else
          Hyperplane_Coefficients(q,fail,h,qd_h);
          if fail then
            return;
          else
            Standard_Linear_Product_System.Add_Hyperplane(natural32(i),h);
            QuadDobl_Linear_Product_System.Add_Hyperplane(natural32(i),qd_h);
            tmp := Tail_Of(tmp);
          end if;
        end if;
      end loop;
    end loop;
  end Store;

  procedure Read_Linear_Product_Start_System
              ( file : in file_type;
                p : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                fail : out boolean ) is

  -- DESCRIPTION :
  --   A random linear product start system is read from file,
  --   using standard complex coefficients.

    lq : Standard_Complex_Prod_Systems.Link_to_Prod_Sys;

  begin
    Standard_Complex_Prod_Systems_io.get(file,lq);
    Store(lq.all,fail);
    if fail then
      put_line("Storing the system as a linear product-system failed!");
      put_line("Please check start system on file and try again later.");
    else
      p := new Standard_Complex_Poly_Systems.Poly_Sys'
                 (Standard_Complex_Prod_Systems.Expand(lq.all));
    end if;
  exception
    when others =>
      put_line("Exception raised while reading linear-product start system.");
      raise;
  end Read_Linear_Product_Start_System;

  procedure Write_Linear_Product_Start_System
              ( file : in file_type; n : in natural32 ) is

    rq : constant Standard_Complex_Prod_Systems.Prod_Sys(1..integer32(n))
       := Standard_Complex_Prod_Planes.Create;

  begin
    Standard_Complex_Prod_Systems_io.put_line(file,n,rq);
  end Write_Linear_Product_Start_System;

  procedure Scan_for_Start_Solutions
              ( file : in file_type; len,dim : out natural32 ) is

    found : boolean;

  begin
    Scan_and_Skip(file,"SOLUTIONS",found);
    len := 0; dim := 0;
    if found
     then get(file,len); get(file,dim);
    end if;
  exception
    when others => put_line("Exception raised while scanning for solutions.");
                   raise;
  end Scan_for_Start_Solutions;

  procedure Read_Systems_and_Track
              ( pfile,qfile,ofile : in file_type;
                f : access procedure ( s : in Solu_Info ) := null ) is

    target,start : Link_to_Poly_Sys;
    ans : character;
    report : boolean;
    len,dim : natural32;

  begin
    Read_Target_System(pfile,target);
    put(ofile,target'last,1);
    new_line(ofile);
    put(ofile,target.all);
    new_line;
    Read_Start_System(qfile,start);
    put("Do you want to monitor progress of solver on screen ? (y/n) "); 
    Ask_Yes_or_No(ans);
    report := (ans = 'y');
    new_line(ofile);
    put_line(ofile,"THE START SYSTEM :");
    put_line(ofile,start.all);
    Scan_for_Start_Solutions(qfile,len,dim);
    Continuation_Parameters.Tune(0,64);
    Driver_for_Continuation_Parameters(ofile);
    new_line;
    Driver_for_Process_io(ofile);
    new_line;
    put_line("See the output file for results...");
    new_line;
    new_line(ofile);
    put_line(ofile,"THE SOLUTIONS :");
    Write_First(ofile,len,dim);
    Cheater_Track(ofile,qfile,report,target.all,start.all,len,dim,f);
  end Read_Systems_and_Track;

  procedure Read_Systems_and_Track
              ( target : in Poly_Sys; qfile,ofile : in file_type;
                kind : in character;
                f : access procedure ( s : in Solu_Info ) := null ) is

    start : Link_to_Poly_Sys;
    start_lps : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    len,dim : natural32;
    ans : character;
    fail,report : boolean;

    use Total_Degree_Start_Systems;

  begin
    if kind = '1' or kind = '3' then
      Read_Start_System(qfile,start);
    else
      Read_Linear_Product_Start_System(qfile,start_lps,fail);
      if fail then return; end if;
      start := new Poly_Sys'
                     (Standard_Poly_Sys_to_QuadDobl_Complex(start_lps.all));
    end if;
    put("Do you want to monitor progress of solver on screen ? (y/n) ");
    Ask_Yes_or_No(ans);
    report := (ans = 'y');
    new_line(ofile);
    put_line(ofile,"THE START SYSTEM :");
    if kind = '1' or kind = '3'
     then put_line(ofile,start.all);
     else Write_Linear_Product_Start_System(ofile,natural32(start_lps'last));
    end if;
    if kind = '1' then
      dim := natural32(start'last);
      len := Product(Degrees(start.all));
    elsif kind = '2' then
      dim := natural32(start_lps'last);
      len := Product(Degrees(start_lps.all));
    else
      Scan_for_Start_Solutions(qfile,len,dim);
    end if;
    new_line;
    Continuation_Parameters.Tune(0,64);
    Driver_for_Continuation_Parameters(ofile);
    new_line;
    Driver_for_Process_io(ofile);
    new_line;
    put_line("See the output file for results...");
    new_line;
    new_line(ofile);
    put_line(ofile,"THE SOLUTIONS :");
    Write_First(ofile,len,dim);
    case kind is
       when '1'
         => Total_Degree_Track(ofile,report,target,start.all,len,dim,f);
       when '2'
         => Linear_Product_Track(ofile,report,target,start.all,len,dim,f);
       when '3'
         => Cheater_Track(ofile,qfile,report,target,start.all,len,dim,f);
       when others => null;
    end case;
  end Read_Systems_and_Track;

end Drivers_to_Track_QuadDobl_Paths;
