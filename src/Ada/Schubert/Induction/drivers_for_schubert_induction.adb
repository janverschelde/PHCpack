with Communications_with_User;           use Communications_with_User;
with File_Scanning,Timing_Package;       use File_Scanning,Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Multprec_Natural_Numbers_io;        use Multprec_Natural_Numbers_io;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Natural_Matrices;          use Standard_Natural_Matrices;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Standard_Random_Matrices;
with Symbol_Table;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
--with Standard_Complex_Poly_SysFun;
with Matrix_Indeterminates;
with Standard_Homotopy;
with Standard_IncFix_Continuation;       use Standard_IncFix_Continuation;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Brackets_io;                        use Brackets_io;
with Bracket_Monomials_io;               use Bracket_Monomials_io;
with Drivers_for_Pieri_Homotopies;       use Drivers_for_Pieri_Homotopies;
with Checker_Boards_io;                  use Checker_Boards_io;
with Checker_Moves;                      use Checker_Moves;
with Checker_Posets,Checker_Posets_io;   use Checker_Posets,Checker_Posets_io;
with Checker_Localization_Patterns;      use Checker_Localization_Patterns;
with Intersection_Posets_io;             use Intersection_Posets_io;
with Moving_Flag_Homotopies;
with Moving_Flag_Continuation;
with Resolve_Schubert_Problems;          use Resolve_Schubert_Problems;

package body Drivers_for_Schubert_Induction is

  function Prompt_for_Bracket_Monomial return Bracket_Monomial is

    res : Bracket_Monomial;

  begin
    put_line("A product of brackets is for example :");
    put_line("  [2 4 6]*[2 4 6]*[2 4 6]; or [2 4 6]^3;");
    put_line("Give a product of brackets, terminate by a semicolon ;");
    get(res);
    return res;
  end Prompt_for_Bracket_Monomial;

  function Is_Isolated 
              ( b : Standard_Natural_Vectors.Vector ) return boolean is
  begin
    for i in b'range loop
      if integer32(b(i)) /= i
       then return false;
      end if;
    end loop;
    return true;
  end Is_Isolated;
 
  function Number_of_Isolated_Solutions
              ( ips : Intersection_Poset ) return Natural_Number is

     res : Natural_Number := Create(integer(0));

  begin
    if not Is_Null(ips.nodes(ips.level)) then
      declare
        psl : constant Poset_List := ips.nodes(ips.level);
        tmp : Poset_List := psl;
        lpnd : Link_to_Poset_Node;
        nd : Link_to_Node;
      begin
        while not Is_Null(tmp) loop
          lpnd := Head_Of(tmp);
          nd := lpnd.ps.white(lpnd.ps.white'last);
          if Is_Isolated(nd.cols)
           then Add(res,nd.coeff);
          end if;
          tmp := Tail_Of(tmp);
        end loop;
      end;
    end if;
    return res;
  end Number_of_Isolated_Solutions;

  procedure Create_Intersection_Poset
              ( n,nb : in integer32; cd : in Array_of_Brackets;
                silent : in boolean; finsum : out Natural_Number ) is

  -- DESCRIPTION :
  --   Creates the intersection poset defined by the nb brackets in cd
  --   for planes in n-space and resolves the intersection condition
  --   imposed by the brackets in cd.

    k : constant integer32 := cd(1)'last;
    p : constant Standard_Natural_Vectors.Vector(1..n)
      := Identity_Permutation(natural32(n));
    r,c,w : Standard_Natural_Vectors.Vector(1..k);
    ps : Poset;
    ips : Intersection_Poset(nb-1);

  begin
    ips.level := 0;
    if not silent then
      put("  the dimension of the planes : "); put(k,1); new_line;
      put("  the number of conditions : "); put(nb,1); new_line;
    end if;
    if nb >= 2 then
      r := Standard_Natural_Vectors.Vector(cd(1).all);
      c := Standard_Natural_Vectors.Vector(cd(2).all);
      if not silent
       then put(cd(1).all); put(" and "); put(cd(2).all);
      end if;
      if not Happy_Checkers(p,r,c) then
        if not silent
         then put_line(" are not happy together.  Please try again.");
        end if;
      else
        if not silent
         then put_line(" form a happy configuration.");
        end if;
        ps := Create(n,r,c);              -- create the first poset
        if not silent
         then Write_Formal_Equation(ps);
        end if;
        ips := Create(nb-1,ps);           -- first poset becomes root
        for i in 3..nb loop
          w := Standard_Natural_Vectors.Vector(cd(i).all);
          Intersect(ips,w,silent);
          if not silent then
            put_line("The new formal equations : ");
            Write_Formal_Equations(ips,ips.level);
          end if;
        end loop;
        if not silent then
          put_line("All formal equations in the intersection poset :");
          Write_Formal_Equations(ips);
          put_line("The intersection condition resolved :");
          Write_Expansion(ips);
        end if;
      end if;
    end if;
    finsum := Final_Sum(ips);
  end Create_Intersection_Poset;

  procedure Create_Intersection_Poset
              ( n : in integer32; bm : in Bracket_Monomial;
                nbsols : out Natural_Number ) is

    nb : constant integer32 := integer32(Number_of_Brackets(bm));
    cd : constant Array_of_Brackets(1..nb) := Create(bm);

  begin
    Create_Intersection_Poset(n,nb,cd,false,nbsols);
  end Create_Intersection_Poset;

  procedure Intersection_Conditions 
               ( bm : in Bracket_Monomial; 
                 rows,cols : out Standard_Natural_Vectors.Link_to_Vector;
                 conds : out Standard_Natural_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Extracts the intersection conditions out of the bracket monomial bm.

  -- ON ENTRY :
  --   bm        bracket monomial stores Schubert problem.

  -- ON RETURN :
  --   rows      first bracket of the bracket monomial bm;
  --   cols      second bracket of the bracket monomial bm;
  --   conds     remainder brackets of bm.

    k,cnt : integer32 := 0;
    extras : Standard_Natural_Matrices.Link_to_Matrix;

    procedure Return_Bracket ( b : in Bracket; continue : out boolean ) is

      bv : constant Standard_Natural_Vectors.Vector(b'range)
         := Standard_Natural_Vectors.Vector(b);

    begin
      cnt := cnt + 1;
     -- put("Bracket "); put(cnt,1); put(" : "); put(b); new_line;
      if cnt = 1 then
        rows := new Standard_Natural_Vectors.Vector'(bv);
        k := bv'last;
      elsif cnt = 2 then
        cols := new Standard_Natural_Vectors.Vector'(bv);
      else
        if cnt = 3 then
          declare
            buffer : Matrix(1..20,1..k);
          begin
            for i in buffer'range(1) loop
              for j in buffer'range(2) loop
                buffer(i,j) := 0;
              end loop;
            end loop;
            extras := new Standard_Natural_Matrices.Matrix'(buffer);
          end;
        end if;
        for i in 1..k loop
          extras(cnt-2,i) := bv(i);
        end loop;
      end if;
      continue := true;
    end Return_Bracket;
    procedure Enum_Brackets is new Enumerate_Brackets(Return_Bracket);

  begin
    Enum_Brackets(bm);
   -- put("k = "); put(k,1);
   -- put(", cnt = "); put(cnt,1); new_line;
    conds := new Standard_Natural_VecVecs.VecVec(1..cnt-2);
    for i in 1..cnt-2 loop
     -- put("copying bracket "); put(i,1); new_line;
      declare
        v : Standard_Natural_Vectors.Vector(1..k);
      begin
        for j in v'range loop
          v(j) := extras(i,j);
        end loop;
        conds(i) := new Standard_Natural_Vectors.Vector'(v);
      end;
    end loop;
  end Intersection_Conditions;

  procedure Resolve_Intersection_Condition ( n : in natural32 ) is

    bm : Bracket_Monomial := Prompt_for_Bracket_Monomial;
    nbsols : Natural_Number;

  begin
    put("Your product : "); put(bm); new_line;
    Create_Intersection_Poset(integer32(n),bm,nbsols);
    Multprec_Natural_Numbers.Clear(nbsols);
    Clear(bm);
  end Resolve_Intersection_Condition;

  function Get_Intersection_Conditions 
             ( k : natural32 )
             return Standard_Natural_VecVecs.Link_to_VecVec is

    res : Standard_Natural_VecVecs.Link_to_VecVec;
    nb : integer32 := 0;

  begin
    put("Give the number of additional intersection conditions : "); get(nb);
    res := new Standard_Natural_VecVecs.VecVec(1..nb);
   -- put("Getting "); put(nb,1); put(" intersection conditions...");
    for i in 1..nb loop
      put("Reading intersection condition "); put(i,1); put_line("...");
      declare
        cond : Standard_Natural_Vectors.Vector(1..integer32(k));
      begin
        Read_Permutation(cond);
        res(i) := new Standard_Natural_Vectors.Vector'(cond);
      end;
    end loop;
    return res;
  end Get_Intersection_Conditions;

  procedure Read_Intersection_Conditions
              ( ip : in Standard_Natural_Vectors.Vector;
                rows,cols : out Standard_Natural_Vectors.Vector ) is
  begin
    new_line;
    loop
      put_line("Reading two intersection conditions...");
      Read_Permutation(rows); Read_Permutation(cols);
      exit when Happy_Checkers(ip,rows,cols);
      put("Your conditions form an unhappy configuration.");
      put_line("  Please try again...");
    end loop;
  end Read_Intersection_Conditions;

  procedure Write_Results
               ( file : in file_type; n,k : in integer32;
                 q,rows,cols : in Standard_Natural_Vectors.Vector;
                 cnds : in Standard_Natural_VecVecs.Link_to_VecVec;
                 vfs : in Standard_Complex_VecMats.VecMat;
                 sols : in Solution_List; fsys : out Link_to_Poly_Sys ) is

    f : Link_to_Poly_Sys;
    mf : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
       := Moving_Flag_Homotopies.Moved_Flag(n);
    locmap : Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,rows,cols);
    dim : constant natural32
        := Checker_Localization_Patterns.Degree_of_Freedom(locmap);
    tmp : Solution_List;
    ls : Link_to_Solution;

  begin
    if Symbol_Table.Number < dim 
     then Matrix_Indeterminates.Initialize_Symbols(dim,locmap);
    end if;
    new_line(file);
    put(file,"Resolved "); put(file,Bracket(rows)); put(file,"*");
    put(file,Bracket(cols));
    put(file," in "); put(file,n,1); put_line(file,"-space,");
    put(file,"with "); put(file,cnds'last,1);
    put(file," fixed conditions : ");
    for i in cnds'range loop
      put(file,Bracket(cnds(i).all));
    end loop;
    put(file," q = "); Standard_Natural_Vectors_io.put(file,q);
    put_line(file," ...");
    new_line(file);
    put_line(file,"THE FIXED FLAGS :");
    for i in vfs'range loop
      put(file,vfs(i).all);
      new_line(file);
    end loop;
    put_line(file,"THE MOVED FLAG :");
    Moving_Flag_Homotopies.Write_Moving_Flag(file,mf);
    if Is_Null(sols) then
      put_line(file,"No solutions found ...");
    else
      new_line(file);
      put(file,"THE "); put(file,k,1);
      put_line(file,"-SOLUTION PLANES :");
      tmp := sols;
      while not Is_Null(tmp) loop
        ls := Head_Of(tmp);
        declare
          x : constant Standard_Complex_Matrices.Matrix(1..n,1..k)
            := Checker_Localization_Patterns.Map(locmap,ls.v);
        begin
          put(file,x);
          new_line(file);
        end;     
        tmp := Tail_Of(tmp);
      end loop;
    end if;
    Moving_Flag_Homotopies.Flag_Conditions(n,k,q,rows,cols,cnds.all,mf,vfs,f);
    put_line(file,"THE POLYNOMIAL SYSTEM :"); put_line(file,f.all);
    if not Is_Null(sols) then
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
    fsys := f;
  end Write_Results;

  function Random_Flags
             ( n,m : integer32 ) return Standard_Complex_VecMats.VecMat is

    res : Standard_Complex_VecMats.VecMat(1..m);

  begin
    for i in res'range loop
      declare
        rf : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
           := Moving_Flag_Homotopies.Random_Flag(n);
      begin
        res(i) := new Standard_Complex_Matrices.Matrix'(rf);
      end;
    end loop;
    return res;
  end Random_Flags;

  procedure Reporting_Moving_Flag_Continuation
              ( n,k : in integer32; tol : in double_float;
                rows,cols : in Standard_Natural_Vectors.Vector;
                cnds : in Standard_Natural_VecVecs.Link_to_VecVec ) is

    file : file_type;

  begin
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    new_line;
    Reporting_Moving_Flag_Continuation(file,n,k,tol,rows,cols,cnds);
  end Reporting_Moving_Flag_Continuation;

  procedure Reporting_Moving_Flag_Continuation
              ( file : in file_type; tune : in boolean;
                n,k : in integer32; tol : in double_float;
                rows,cols : in Standard_Natural_Vectors.Vector;
                cnds : in Standard_Natural_VecVecs.Link_to_VecVec;
                sols : out Solution_List; fsys : out Link_to_Poly_Sys;
                flags : out Standard_Complex_VecMats.VecMat ) is

    ip : constant Standard_Natural_Vectors.Vector(1..n)
       := Identity_Permutation(natural32(n));
    ps : Poset;
    report : boolean;
    timer : Timing_Widget;

  begin
    if tune
     then Moving_Flag_Continuation.Set_Parameters(file,report);
    end if;
    ps := Create(n,rows,cols);
    flags := Random_Flags(n,cnds'last);
    tstart(timer);
    Moving_Flag_Continuation.Track_All_Paths_in_Poset
      (file,n,k,ps,cnds.all,flags,tol,sols);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"tracking all paths");
    Write_Results(file,n,k,ip,rows,cols,cnds,flags,sols,fsys);
  end Reporting_Moving_Flag_Continuation;

  procedure Reporting_Moving_Flag_Continuation
              ( file : in file_type;
                n,k : in integer32; tol : in double_float;
                rows,cols : in Standard_Natural_Vectors.Vector;
                cnds : in Standard_Natural_VecVecs.Link_to_VecVec ) is

    sols : Solution_List;
    fsys : Link_to_Poly_Sys;
    flgs : Standard_Complex_VecMats.VecMat(cnds'range);

  begin
    Reporting_Moving_Flag_Continuation
      (file,true,n,k,tol,rows,cols,cnds,sols,fsys,flgs);
  end Reporting_Moving_Flag_Continuation;

  procedure Run_Moving_Flag_Continuation ( n,k : in integer32 ) is

    p : constant Standard_Natural_Vectors.Vector(1..n)
      := Identity_Permutation(natural32(n));
    rows,cols,cond : Standard_Natural_Vectors.Vector(1..k);
    cnds : Standard_Natural_VecVecs.Link_to_VecVec;
    happy : boolean;
    tol : constant double_float := 1.0E-5;

  begin
    new_line;
    loop
      put_line("Reading rows and columns for white checkers ...");
      Read_Permutation(rows);
      Read_Permutation(cols);
      Check_Happiness(p,rows,cols,happy);
      exit when happy;
      put_line("The white checkers are not happy, please try again!");
    end loop;
    put_line("Reading third intersection condition ...");
    Read_Permutation(cond);
    cnds := new Standard_Natural_VecVecs.VecVec(1..1);
    cnds(1) := new Standard_Natural_Vectors.Vector'(cond);
    skip_line; -- skip enter symbol
    Reporting_Moving_Flag_Continuation(n,k,tol,rows,cols,cnds);
  end Run_Moving_Flag_Continuation;

  procedure Set_Symbol_Table
              ( n,k : in integer32;
                p,rows,cols : in Standard_Natural_Vectors.Vector ) is

    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Column_Pattern(n,k,p,rows,cols);
    dim : constant natural32 := Degree_of_Freedom(locmap);

  begin
    Matrix_Indeterminates.Initialize_Symbols(dim,locmap);
  end Set_Symbol_Table;

  procedure Scan_for_Start_Schubert_Problem
              ( file : in file_type; n : in integer32;
                vf : out Standard_Complex_VecMats.VecMat;
                p : out Link_to_Poly_Sys; sols : out Solution_List;
                fail : out boolean ) is

    found : boolean;

  begin
    new_line;
    put("Scanning the file for "); put(vf'last,1);
    put_line(" fixed flags ...");
    Scan_and_Skip(file,"THE FIXED FLAGS :",found);
    if not found then
      put_line(" failed!?!");
      fail := true;
    else
      for i in vf'range loop
        put("  reading flag "); put(i,1); put(" ...");
        declare
          flag : Standard_Complex_Matrices.Matrix(1..n,1..n);
        begin
          get(file,flag);
          vf(i) := new Standard_Complex_Matrices.Matrix'(flag);
          put_line(" okay.");
        end;
      end loop;
      put("Scanning the file for start system ...");
      Scan_and_Skip(file,"THE POLYNOMIAL SYSTEM :",found);
      if not found then
        put_line(" failed!?!");
        fail := true;
      else
        put_line(" okay.");
        get(file,p);
        put("Scanning the file for solutions ...");
        Scan_and_Skip(file,"THE SOLUTIONS :",found);
        if not found then
          put_line(" failed!?!");
          fail := true;
        else
          put_line(" okay.");
          get(file,sols);
          fail := false;
        end if;
      end if;
    end if;
  end Scan_for_Start_Schubert_Problem;

  function Target_Flags
             ( n,m,t : in integer32 ) 
             return Standard_Complex_VecMats.VecMat is

  -- DESCRIPTION :
  --   Returns a vector of range 1..m of n-by-n matrices, depending on t
  --     t = 0 : the user will be promped to provide the flag coordinates,
  --     t = 1 : complex random numbers are generated for the flags,
  --     t = 2 : real random numbers are generated for the flags.

    res : Standard_Complex_VecMats.VecMat(1..m);

  begin
    for i in 1..m loop
      declare
        flag : Standard_Complex_Matrices.Matrix(1..n,1..n);
      begin
        if t = 0 then
          put("Give "); put(n,1); put("-by-"); put(n,1);
          put(" complex matrix :"); get(flag);
        else
          flag := Standard_Random_Matrices.Random_Matrix
                    (natural32(n),natural32(n));
          if t = 2 then
            for i in flag'range(1) loop
              for j in flag'range(2) loop
                flag(i,j) := Create(REAL_PART(flag(i,j)),0.0);
              end loop;
            end loop;
          end if;
          for i in flag'range(1) loop
            flag(i,i) := Create(1.0);
            for j in i+1..flag'last(2) loop
              flag(i,j) := Create(0.0);
            end loop;
          end loop;
        end if;
        res(i) := new Standard_Complex_Matrices.Matrix'(flag);
      end;
    end loop;
    return res;
  end Target_Flags;

  procedure Run_Cheater_Continuation
              ( file : in file_type; h : in Poly_Sys;
                sols : in out Solution_List ) is

    procedure Rep_Cont is
      new Reporting_Continue(Max_Norm,Standard_Homotopy.Eval,
                             Standard_Homotopy.Diff,Standard_Homotopy.Diff);

  begin
    Standard_Homotopy.Create(h,h'last+1);
    Set_Continuation_Parameter(sols,Create(0.0));
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(sols),1); put(file,"  ");
    put(file,Head_Of(sols).n,1); new_line(file);
    Rep_Cont(file,sols,false,Create(1.0));
  end Run_Cheater_Continuation;

  procedure Run_Cheater
              ( file : in file_type; n,k : in integer32;
                rows,cols : in Standard_Natural_Vectors.Vector;
                cnds : in Standard_Natural_VecVecs.Link_to_VecVec;
                stvf,tgvf : in Standard_Complex_VecMats.VecMat;
                sols : in out Solution_List ) is

    q : constant Standard_Natural_Vectors.Vector(1..n)
      := Identity_Permutation(natural32(n));
    cnd : constant Bracket(1..k) := Bracket(cnds(cnds'first).all);
    stf : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
        := stvf(stvf'first).all;
    tgf : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
        := tgvf(tgvf'first).all;
    nv : constant integer32 := Head_Of(sols).n;
    hom : Link_to_Poly_Sys;
   -- flagpols : Standard_Complex_Poly_Matrices.Matrix
   --          := Moving_Flag_Homotopies.Cheater_Homotopy_Flag(nv,stf,tgf);
   -- a0 : Standard_Complex_Vectors.Vector(1..nv+1) := (1..nv+1 => Create(0.0));
   -- b1 : Standard_Complex_Vectors.Vector(1..nv+1) := (1..nv+1 => Create(0.0));
   -- a,b : Standard_Complex_Matrices.Matrix(1..n,1..n);

  begin
   -- b1(nv+1) := Create(1.0);
   -- for i in flagpols'range(1) loop
   --   for j in flagpols'range(2) loop
   --     a(i,j) := Standard_Complex_Poly_Functions.Eval(flagpols(i,j),a0);
   --     b(i,j) := Standard_Complex_Poly_Functions.Eval(flagpols(i,j),b1);
   --   end loop;
   -- end loop;
   -- put_line(file,"the start flags : ");  put(file,a);
   -- put_line(file,"the target flags : "); put(file,b);
    Moving_Flag_Homotopies.Add_t_Symbol;
    Moving_Flag_Homotopies.Flag_Conditions(n,k,q,rows,cols,cnd,stf,tgf,hom);
   -- put_line(file,"the cheater's homotopy : ");
   -- put_line(file,hom.all);
    declare
     -- sth : constant Poly_Sys(hom'range)
     --     := Standard_Complex_Poly_SysFun.Eval(hom.all,Create(0.0),nv+1);
     -- tgh : constant Poly_Sys(hom'range)
     --     := Standard_Complex_Poly_SysFun.Eval(hom.all,Create(1.0),nv+1);
      squhom : Poly_Sys(1..nv) := Square(natural32(nv),hom.all);
    begin
     -- put_line(file,"the start system :"); put_line(file,sth);
     -- put_line(file,"the target system :"); put_line(file,tgh);
      Run_Cheater_Continuation(file,squhom,sols);
      Standard_Complex_Poly_Systems.Clear(squhom);
    end;
  end Run_Cheater;

  procedure Run_Cheater_Flag_Homotopy
              ( n,k : in integer32;
                rows,cols : in Standard_Natural_Vectors.Vector;
                cnds : in Standard_Natural_VecVecs.Link_to_VecVec;
                inpt : in boolean ) is

    infile,outfile : file_type;
    start,f : Link_to_Poly_Sys;
    startsols : Solution_List;
    stvf,tgvf : Standard_Complex_VecMats.VecMat(cnds'range);
    q : constant Standard_Natural_Vectors.Vector(1..n)
      := Identity_Permutation(natural32(n));
    fail : boolean;
    ans : character;

  begin
    new_line;
    put_line("Reading the name of an input file for start problem ...");
    Read_Name_and_Open_File(infile);
    Set_Symbol_Table(n,k,q,rows,cols);
    Scan_for_Start_Schubert_Problem(infile,n,stvf,start,startsols,fail);
    if not fail then
      new_line;
      if inpt then
        tgvf := Target_Flags(n,stvf'last,0);
      else
        put("Default is complex, generate real flags (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'y'
         then tgvf := Target_Flags(n,stvf'last,2);
         else tgvf := Target_Flags(n,stvf'last,1);
        end if;
      end if;
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(outfile);
      put_line(outfile,"THE START FIXED FLAGS :");
      for i in stvf'range loop
        put(outfile,stvf(i).all);
      end loop;
      new_line(outfile);
      put_line(outfile,"THE START POLYNOMIAL SYSTEM :");
      put_line(outfile,start.all);
      new_line(outfile);
      put_line(outfile,"THE START SOLUTIONS :");
      put(outfile,Length_Of(startsols),
          natural32(Head_Of(startsols).n),startsols);
      new_line(outfile);
      put_line(outfile,"THE FIXED FLAGS :"); 
      for i in tgvf'range loop
        put(outfile,tgvf(i).all);
      end loop;
      Moving_Flag_Homotopies.Flag_Conditions(n,k,q,rows,cols,cnds.all,tgvf,f);
      new_line(outfile);
      put_line(outfile,"THE POLYNOMIAL SYSTEM :");
      put_line(outfile,f.all);
     -- A plain cheater in coefficient space does not work!
     -- Path_Following_with_Cheater(outfile,start.all,f.all,startsols,true);
     -- Instead we run cheater in the parameter space of flags :
      Run_Cheater(outfile,n,k,rows,cols,cnds,stvf,tgvf,startsols);
    end if;
  end Run_Cheater_Flag_Homotopy;

  function Process_Conditions
             ( n,k,m : integer32; conds : Array_of_Brackets )
             return Intersection_Posets.Intersection_Poset is

  -- DESCRIPTION :
  --   Process the m conditions stored in conds on k-planes in n-space.
  --   Returns the intersection poset.

    use Intersection_Posets;
    res : Intersection_Poset(m-1);
    p : constant Standard_Natural_Vectors.Vector(1..n)
      := Identity_Permutation(natural32(n));
    rows,cols : Standard_Natural_Vectors.Vector(1..k);
    ps : Poset;
    ans : character;
    silent : boolean;
    top_roco : Natural_Number;

  begin
    new_line;
    put("Intermediate output during formal root count ? (y/n) "); 
    Ask_Yes_or_No(ans);
    silent := (ans = 'n');
   -- put_line("Reading the first two intersection conditions...");
    rows := Standard_Natural_Vectors.Vector(conds(1).all);
    cols := Standard_Natural_Vectors.Vector(conds(2).all);
   -- Read_Permutation(rows); Read_Permutation(cols);
    if not Happy_Checkers(p,rows,cols) then
      put_line("Your conditions form an unhappy configuration.");
    else
      ps := Create(n,rows,cols);
      res := Create(m-1,ps);
      for k in 3..m loop
       -- put("Reading intersection condition "); put(k,1); put_line("...");
       -- Read_Permutation(cols);
        cols := Standard_Natural_Vectors.Vector(conds(k).all);
        Intersect(res,cols,silent);
      end loop;
    end if;
    top_roco := Final_Sum(res);
    put("The formal root count : "); put(top_roco); new_line;
    Clear(top_roco);
    return res;
  end Process_Conditions;

  function Bracket_to_Vector
             ( b : Bracket ) return Standard_Natural_Vectors.Vector is

  -- DESCRIPTION :
  --   Converts the bracket to a vector of standard natural numbers.

    res : Standard_Natural_Vectors.Vector(b'range);

  begin
    for i in b'range loop
      res(i) := b(i);
    end loop;
    return res;
  end Bracket_to_Vector;

  function Remaining_Intersection_Conditions
             ( b : Array_of_Brackets )
             return Standard_Natural_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns the remaining b'last-2 conditions stored in b.

  -- REQUIRED : b'last > 2.

    res : Standard_Natural_VecVecs.VecVec(b'first..b'last-2);

  begin
    for i in b'first+2..b'last loop
      declare
        bck : constant Link_to_Bracket := b(i);
        bvc : Standard_Natural_Vectors.Vector(bck'range);
      begin
        for j in bvc'range loop
          bvc(j) := bck(j);
        end loop;
        res(i-2) := new Standard_Natural_Vectors.Vector'(bvc);
      end;
    end loop;
    return res;
  end Remaining_Intersection_Conditions;            

  procedure Resolve_Schubert_Problem
              ( n,k : in integer32; bm : in Bracket_Monomial ) is

  -- DESCRIPTION :
  --   Prompts the user for m intersection conditions on k-planes in n-space,
  --   and writes the evolution of the root count from the leaves to the root.

  -- ON ENTRY :
  --   n        ambient space
  --   k        dimension of the solution planes;
  --   bm       product of k-brackets, with conditions on the k-planes.

    use Intersection_Posets;

    file : file_type;
    cnd : constant Array_of_Brackets := Create(bm);
    nbc : constant integer32 := cnd'last;
    ips : Intersection_Poset(nbc-1) := Process_Conditions(n,k,nbc,cnd);
    sps : Solution_Poset(ips.m) := Create(ips);
    top_roco,bottom_roco : Natural_Number;
    q : constant Standard_Natural_Vectors.Vector
      := Identity_Permutation(natural32(n));
    rows : constant Standard_Natural_Vectors.Vector
         := Bracket_to_Vector(cnd(cnd'first).all);
    cols : constant Standard_Natural_Vectors.Vector
         := Bracket_to_Vector(cnd(cnd'first+1).all);
    conds : Standard_Natural_VecVecs.VecVec(1..nbc-2)
          := Remaining_Intersection_Conditions(cnd);
    link2conds : constant Standard_Natural_VecVecs.Link_to_VecVec
               := new Standard_Natural_VecVecs.VecVec'(conds);
    flags : Standard_Complex_VecMats.VecMat(1..nbc-2)
          := Random_Flags(n,nbc-2);
    fsys : Link_to_Poly_Sys;
    sols : Solution_List;
    tol : constant double_float := 1.0E-6;
    ans : character;
    monitor_games : boolean;
    timer : Timing_Widget;

  begin
    new_line;
    put_line("Reading a name for the output file ...");
    Read_Name_and_Create_File(file);
    new_line;
    put("Monitor Littlewood-Richardson homotopies"
      & " in each checker game ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    monitor_games := (ans = 'y');
    top_roco := Final_Sum(ips);
    put("The formal root count : "); put(top_roco); new_line;
    put_line("... running the root counting from the bottom up ...");
    Count_Roots(file,ips,bottom_roco);
    put(" Top down root count : "); put(top_roco); new_line;
    put("Bottom up root count : "); put(bottom_roco); new_line;
    put_line("... resolving the Schubert problem ...");
    new_line;
    put_line("See the output file for results ...");
    new_line;
    tstart(timer);
    Resolve(file,monitor_games,n,k,tol,ips,sps,conds,flags,sols);
    tstop(timer);
    Write_Results(file,n,k,q,rows,cols,link2conds,flags,sols,fsys);
    new_line(file);
    print_times(file,timer,"resolving a Schubert problem");
  end Resolve_Schubert_Problem;

  procedure Solve_Schubert_Problems ( n : in integer32 ) is

    k : integer32;
    bm : Bracket_Monomial := Prompt_for_Bracket_Monomial;
    rows,cols : Standard_Natural_Vectors.Link_to_Vector;
    cnds : Standard_Natural_VecVecs.Link_to_VecVec;
    nbsols : Natural_Number;
    ans : character;
    inpt : boolean;
    tol : double_float := 1.0E-5;

  begin
    Intersection_Conditions(bm,rows,cols,cnds);
    k := rows'last;
    new_line;
    put_line("MENU for Littlewood-Richardson homotopies :");
    put_line("  0. solve a generic instance for random flags;");
    put_line("  1. run a cheater's homotopy to other random flags;");
    put_line("  2. solve a specific instance via cheater to given flags.");
    put("Type 0, 1, or 2 to select from menu : ");
    Ask_Alternative(ans,"012");
    if ans = '0' then
      Resolve_Schubert_Problem(n,k,bm);
    else
      new_line;
      put_line("resolving the intersection conditions ...");
      Create_Intersection_Poset(n,bm,nbsols);
      put("Number of isolated solutions : "); put(nbsols); new_line;
      if nbsols > 0 then
        inpt := (ans = '2');
        Run_Cheater_Flag_Homotopy(n,k,rows.all,cols.all,cnds,inpt);
      end if;
    end if;
    Clear(bm);
  end Solve_Schubert_Problems;

end Drivers_for_Schubert_Induction;
