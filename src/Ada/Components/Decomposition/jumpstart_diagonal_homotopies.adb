with Communications_with_User;           use Communications_with_User;
with File_Scanning;                      use File_Scanning;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Main_Poly_Continuation;             use Main_Poly_Continuation;
with Drivers_to_Track_Standard_Paths;
with Permutations,Permute_Operations;    use Permutations,Permute_Operations;
with Witness_Sets,Witness_Sets_io;       use Witness_Sets,Witness_Sets_io;
with Extrinsic_Diagonal_Homotopies;      use Extrinsic_Diagonal_Homotopies;
with Extrinsic_Diagonal_Homotopies_io;   use Extrinsic_Diagonal_Homotopies_io;

package body Jumpstart_Diagonal_Homotopies is

  function Ask_for_Monitor return boolean is

    ans : character;

  begin
    new_line;
    put("Do you want to see the progress of the tracking on screen ? (y/n) ");
    Ask_Yes_or_No(ans);
    return (ans = 'y');
  end Ask_for_Monitor;

  procedure Read_Degree_of_Witness_Set
              ( file : in file_type; d,n : out natural32 ) is

    found : boolean;

  begin
   -- put_line("  looking for solutions ...");
    Scan_and_Skip(file,"SOLUTIONS",found);
    if not found
     then put_line("found no solutions in the input file");
     else Read_First(file,d,n);
         -- put("degree of the witness set : "); put(d,1); new_line;
    end if;
  end Read_Degree_of_Witness_Set;

  procedure Read_Witness_Set 
              ( file : out file_type; lp : out Link_to_Poly_Sys;
	        n,k,d : out natural32 ) is
  begin
    new_line;
    put_line("Reading the name of a file for a witness set.");
    Read_Name_and_Open_File(file);
    n := 0;
    get(file,n); skip_line(file); -- skip optional #variables on the line
    new_line;
    put("total number of equations and variables : "); put(n,1); new_line;
   -- put_line("  reading the polynomial system ...");
    lp := new Poly_Sys(1..integer32(n));
    Symbol_Table.Init(n);
    get(file,lp.all);
    k := Count_Embed_Symbols(n,"zz");
    Swap_Symbols_to_End(n,k,"zz",lp.all);
    put("dimension of the witness set : "); put(k,1); new_line;
    Read_Degree_of_Witness_Set(file,d,n);
  end Read_Witness_Set;

  procedure Test_Read_Next
               ( file : in file_type; p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Interactive test on the reading of the solutions in a witness set.
  --   The system on input is an embedded system.

  -- REQUIRED : 
  --   The file is properly positioned after the banner,
  --   ready for reading the first solution.

    n : constant natural32 := natural32(p'last);
    ls : Link_to_Solution;
    cnt : natural32 := 0;
    y : Standard_Complex_Vectors.Vector(p'range);
    ans : character;

  begin
    loop
      cnt := cnt + 1;
      new_line;
      put("Solution "); put(cnt,1); put_line(" of the 1st witness set :");
      Read_Next(file,n,ls); put(ls.all); new_line;
      y := Eval(p,ls.v);
      put_line("Residual vector :"); put_line(y);
      put("Do you want to see more solutions ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_Read_Next;

  procedure Test_Read_Next
               ( wf1,wf2 : in file_type; n1,n2 : in natural32;
                 p1,p2 : in Poly_Sys ) is

    ls : Link_to_Solution;
    cnt : natural32 := 0;
    y1 : Standard_Complex_Vectors.Vector(p1'range);
    y2 : Standard_Complex_Vectors.Vector(p2'range);
    ans : character;

  begin
    loop
      cnt := cnt + 1;
      new_line;
      put("Solution "); put(cnt,1); put_line(" of the 1st witness set :");
      Read_Next(wf1,n1,ls); put(ls.all); new_line;
      y1 := Eval(p1,ls.v);
      put_line("Residual vector :"); put_line(y1);
      put("Do you want to see more solutions ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    cnt := 0;
    loop
      cnt := cnt + 1;
      new_line;
      put("Solution "); put(cnt,1); put_line(" of the 2nd witness set :");
      Read_Next(wf2,n2,ls); put(ls.all); new_line;
      y2 := Eval(p2,ls.v);
      put_line("Residual vector :"); put_line(y2);
      put("Do you want to see more solutions ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_Read_Next;

  procedure Read_Two_Witness_Sets 
               ( wf1,wf2 : out file_type; lp1,lp2 : out Link_to_Poly_Sys;
                 n1,n2,k1,k2,d1,d2 : out natural32;
                 sym1,sym2 : out Symbol_Table.Link_to_Array_of_Symbols ) is

  begin
    Read_Witness_Set(wf1,lp1,n1,k1,d1);
   -- Test_Read_Next(wf1,lp1.all);
    sym1 := Extrinsic_Diagonal_Homotopies_io.Get_Link_to_Symbols;
    Read_Witness_Set(wf2,lp2,n2,k2,d2);
    sym2 := Extrinsic_Diagonal_Homotopies_io.Get_Link_to_Symbols;
   -- Test_Read_Next(wf1,wf2,n1,n2,lp1.all,lp2.all);
   -- put("Symbols for 1st witness set : "); Write(sym1.all);
   -- put("Symbols for 2nd witness set : "); Write(sym2.all);
  end Read_Two_Witness_Sets;

  procedure Match_Symbols
              ( n1,n2,k1,k2 : in natural32; p2 : in out Poly_Sys;
                ls1,ls2 : in Symbol_Table.Link_to_Array_of_Symbols;
                solsym : out Symbol_Table.Link_to_Array_of_Symbols ) is

    n : constant natural32 := n1 - k1;
    prm : Permutation(1..integer32(n));
    sa : Symbol_Table.Array_of_Symbols(1..integer32(n2));
    full_prm,inv_prm : Permutation(1..integer32(n2));
    pp2 : Poly_Sys(p2'range);

  begin
    prm := Match_Symbols(ls1(1..integer32(n)),ls2(1..integer32(n)));
   -- put("The matching permutation : ");
   -- put(Standard_Integer_Vectors.Vector(prm)); new_line;
    full_prm(1..integer32(n)) := prm;
    for i in integer32(n)+1..integer32(n2) loop
      full_prm(i) := i;
    end loop;
    inv_prm := inv(full_prm);
    pp2 := p2*inv_prm;
    Copy(pp2,p2); Clear(pp2);
    Symbol_Table.Clear;
    if ls1'last >= ls2'last then
      Assign_Symbol_Table(ls1.all);
      solsym := ls1;
    else
      sa := ls2.all;
      Permute(full_prm,sa);
      Assign_Symbol_Table(sa);
      solsym := new Symbol_Table.Array_of_Symbols'(sa);
    end if;
  exception
    when others =>
       put_line("Exception raised in Jumpstart_Diagonal_Homotopies");
       put_line("in Match_Symbols"); raise;
  end Match_Symbols;

  procedure Reset_to_Reread_Solutions
              ( file : in out file_type; d,n : out natural32 ) is

    found : boolean;

  begin
    Reset(file);
    Scan_and_Skip(file,"SOLUTIONS",found);
    Read_First(file,d,n);
   -- put("Ready to read "); put(d,1);
   -- put(" solutions of dimension "); put(n,1); put_line(" ...");
  end Reset_to_Reread_Solutions;

  procedure Start_Cascade
              ( file,a_f : in file_type; b_f : in out file_type;
                report : in boolean; start,target : in Poly_Sys;
                a_n,b_n,a_dim,b_dim,a_deg,b_deg : in natural32;
                sbt : in Symbol_Table.Link_to_Array_of_Symbols ) is

    ls1,ls2 : Link_to_Solution;
    s1 : Solution(integer32(a_n-a_dim));
    s2 : Solution(integer32(b_n-b_dim));
    dim : constant natural32 := natural32(start'last);
    len : constant natural32 := a_deg*b_deg;
    cnt : natural32 := 0;
    ind1 : natural32 := 0;
    ind2 : natural32 := 0;
    ps : Solution(integer32(dim));
    dd,nn : natural32;

    procedure Get_Next ( n : in natural32; ls : out Link_to_Solution;
                         ind : out natural32; continue : out boolean ) is

    -- DESCRIPTION :
    --   Reads the next solution from file to form the product.

    --  y : Standard_Complex_Vectors.Vector(1..n);

    begin
      if cnt = 0 then
        Read_Next(a_f,a_n,ls1,sbt.all); ind1 := 1;
        s1 := Remove_Embedding(ls1.all,a_dim);
      end if;
      cnt := cnt + 1;
      ind2 := ind2 + 1;
      Read_Next(b_f,b_n,ls2,sbt.all);
      s2 := Remove_Embedding(ls2.all,b_dim);
      ps := Extrinsic_Product(a_dim,b_dim,s1,s2);
      ls := new Solution'(ps);
     -- y := Eval(start,ls.v);
     -- put("The residual vector at solution "); put(cnt,1); put_line(":");
     -- put_line(y);
      ind := cnt;
      if (cnt mod b_deg = 0) and (cnt < len) then
        if report
         then put("Resetting second solution file for next solution ");
              put(cnt+1,1); put_line(".");
        end if;
        Reset_to_Reread_Solutions(b_f,dd,nn);
        Read_Next(a_f,a_n,ls1,sbt.all);
        s1 := Remove_Embedding(ls1.all,a_dim);
        ind1 := ind1 + 1; ind2 := 0;
      end if;
      if report
       then put("Solution path "); put(ind1,1); put("x"); put(ind2,1);
            put(" = "); put(cnt,1);
      end if;
      continue := true;
    end Get_Next;
    procedure Track_Paths is
      new Drivers_to_Track_Standard_Paths.track(Get_Next);

    procedure Enumerate_Products is -- intermediate try out
    begin
      for i in 1..a_deg loop
        Read_Next(a_f,a_n,ls1);
        s1 := Remove_Embedding(ls1.all,a_dim);
        for j in 1..b_deg loop
          Read_Next(b_f,b_n,ls2);
          s2 := Remove_Embedding(ls2.all,b_dim);
          ps := Extrinsic_Product(a_dim,b_dim,s1,s2);
          cnt := cnt + 1;
          put("Solution path "); put(i,1); put("x"); put(j,1);
          put(" = "); put(cnt,1);
          put_line(" ... ");
        end loop;
        if i < a_deg
         then put("Resetting second solution file at solution ");
              put(cnt,1); put_line(".");
              Reset_to_Reread_Solutions(b_f,dd,nn);
        end if;
      end loop;     
    end Enumerate_Products;

  begin
   -- Enumerate_Products;
    Track_Paths(file,report,target,start,len,dim);
  end Start_Cascade;

  procedure Intersect_Witness_Sets
              ( file,a_f : in file_type;
                b_f : in out file_type; a_p,b_p : in Poly_Sys;
                a_n,b_n,a_dim,b_dim,a_deg,b_deg : in natural32;
                s1e,s2e,sbt : in Symbol_Table.Link_to_Array_of_Symbols ) is

    use Symbol_Table;

    cd : constant natural32 := Cascade_Dimension(a_p,b_p,a_dim,b_dim);
    s1 : constant Array_of_Symbols := Remove_Embed_Symbols(sbt.all);
    s11 : constant Array_of_Symbols(s1'range) := Add_Suffix(s1,'1');
    s22 : constant Array_of_Symbols(s1'range) := Add_Suffix(s1,'2');
    start,target : Poly_Sys(1..integer32(cd));
    oc : natural32;
    report : boolean;

  begin
    new_line;
    put("Intersecting sets of dimension ");
    put(a_dim,1); put(" and "); put(b_dim,1); 
    put(" of degrees "); put(a_deg,1); put(" and "); put(b_deg,1);
    put_line(".");
    put("The dimension of the cascade : "); put(cd,1); new_line;
    Extrinsic_Cascade_Homotopy(a_p,b_p,a_dim,b_dim,start,target);
    Symbol_Table.Clear;
    Assign_Symbol_Table(s11,s22);
    Add_Embed_Symbols(cd-Symbol_Table.Number);
    put_line(file,target);
    new_line(file);
    put_line(file,"TITLE : target system to start diagonal homotopy cascade");
    new_line(file);
    put_line(file,start);
    new_line(file);
    put_line(file,"TITLE : start system to start diagonal homotopy cascade");
    new_line(file);
    new_line;
    Driver_for_Continuation_Parameters(file);
    new_line;
    Driver_for_Process_io(file,oc);
    report := Ask_for_Monitor;
    new_line;
    put_line("See the output file for results ...");
    new_line;
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    Standard_Complex_Solutions_io.Write_First(file,a_deg*b_deg,cd);
    Start_Cascade
      (file,a_f,b_f,report,start,target,
       a_n,b_n,a_dim,b_dim,a_deg,b_deg,sbt);
  end Intersect_Witness_Sets;

  procedure Jumpstart_Diagonal_Homotopy is

    outfile,infile1,infile2 : file_type;
    lp1,lp2 : Link_to_Poly_Sys;
    n1,n2,k1,k2,d1,d2 : natural32;
    s1,s2,ss : Symbol_Table.Link_to_Array_of_Symbols;

  begin
    new_line;
    put_line("Jumpstarting a diagonal homotopy to intersect algebraic sets.");
    Read_Two_Witness_Sets(infile1,infile2,lp1,lp2,n1,n2,k1,k2,d1,d2,s1,s2);
    Match_Symbols(n1,n2,k1,k2,lp2.all,s1,s2,ss);
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(outfile);
    if k1 < k2 then
      Intersect_Witness_Sets
        (outfile,infile2,infile1,lp2.all,lp1.all,n2,n1,k2,k1,d2,d1,s2,s1,ss);
     else
      Intersect_Witness_Sets
        (outfile,infile1,infile2,lp1.all,lp2.all,n1,n2,k1,k2,d1,d2,s1,s2,ss);
    end if;
  end Jumpstart_Diagonal_Homotopy;

  procedure Jumpstart_Diagonal_Homotopy
                ( infile : in out file_type;
                  outfile : in file_type; p : in out Poly_Sys ) is

    infile2 : file_type;
    lp2 : Link_to_Poly_Sys;
    n1,n2,k1,k2,d1,d2 : natural32;
    s1,s2,ss : Symbol_Table.Link_to_Array_of_Symbols;

  begin
    n1 := natural32(p'last);
    k1 := Count_Embed_Symbols(n1,"zz");
    Swap_Symbols_to_End(n1,k1,"zz",p);
    new_line;
    put("dimension of the first witness set : "); put(k1,1); new_line;
    Read_Degree_of_Witness_Set(infile,d1,n1);
    s1 := Extrinsic_Diagonal_Homotopies_io.Get_Link_to_Symbols;
    Read_Witness_Set(infile2,lp2,n2,k2,d2);
    s2 := Extrinsic_Diagonal_Homotopies_io.Get_Link_to_Symbols;
    Match_Symbols(n1,n2,k1,k2,lp2.all,s1,s2,ss);
    if k1 < k2 then
      Intersect_Witness_Sets
        (outfile,infile2,infile,lp2.all,p,n2,n1,k2,k1,d2,d1,s2,s1,ss);
     else
      Intersect_Witness_Sets
        (outfile,infile,infile2,p,lp2.all,n1,n2,k1,k2,d1,d2,s1,s2,ss);
    end if;
  end Jumpstart_Diagonal_Homotopy;

  procedure Track_Paths_Down
               ( infile,outfile : in file_type;
                 start : in Poly_Sys; n,k,d : in natural32 ) is

    target : constant Poly_Sys(start'range) := Remove_Slice(start);
    oc : natural32;
    cnt : natural32 := 0;
    tol : constant double_float := 1.0E-8;
    report : boolean;

    procedure Get_Next ( n : in natural32; ls : out Link_to_Solution;
                         ind : out natural32; continue : out boolean ) is

    -- DESCRIPTION :
    --   Reads the next solution from file to form the product.

      abz : double_float;

    begin
      cnt := cnt + 1;
      ind := cnt;
      Read_Next(infile,n,ls);
      abz := AbsVal(ls.v(ls.v'last));
      if report
       then put("Start solution with |z| ="); put(abz,3);
            if abz > tol
             then put(" > ");  put(tol,1); put(" path ");
             else put(" <= "); put(tol,1); put(" no path ");
            end if;
      end if;
      if abz < tol
       then ls.t := Create(1.0);
      end if;
      if report
       then put(cnt,1);
      end if;
      continue := true;
    end Get_Next;
    procedure Track_Paths is
      new Drivers_to_Track_Standard_Paths.track(Get_Next);

  begin
    put_line(outfile,target);
    new_line(outfile);
    put_line(outfile,"TITLE : target system in the cascade one level down");
    new_line(outfile);
    new_line;
    Driver_for_Continuation_Parameters(outfile);
    new_line;
    Driver_for_Process_io(outfile,oc);
    report := Ask_for_Monitor;
    new_line;
    put_line("See the output file for results ...");
    new_line;
    new_line(outfile);
    put_line(outfile,"THE SOLUTIONS :");
    Standard_Complex_Solutions_io.Write_First(outfile,d,n);
    Track_Paths(outfile,report,target,start,d,n);
  end Track_Paths_Down;

  procedure Jumpstart_Cascade_Homotopy is

    infile,outfile : file_type;
    lp : Link_to_Poly_Sys;
    n,k,d : natural32;

  begin
    new_line;
    put_line("Jumpstarting a cascade homotopy to go down one level.");
    Read_Witness_Set(infile,lp,n,k,d);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(outfile);
   -- Test_Read_Next(infile,lp.all);
    Track_Paths_Down(infile,outfile,lp.all,n,k,d);
  end Jumpstart_Cascade_Homotopy;

  procedure Jumpstart_Cascade_Homotopy
               ( infile,outfile : in file_type; p : in out Poly_Sys ) is

    k,d,n : natural32;

  begin
    new_line;
    put_line("Jumpstarting a cascade homotopy to go down one level.");
    n := natural32(p'last);
    k := Count_Embed_Symbols(n,"zz");
    Swap_Symbols_to_End(n,k,"zz",p);
    put("dimension of the witness set : "); put(k,1); new_line;
    Read_Degree_of_Witness_Set(infile,d,n);
    Track_Paths_Down(infile,outfile,p,n,k,d);
  end Jumpstart_Cascade_Homotopy;

  procedure Remove_Last_Slack_Variable
               ( infile,outfile : in file_type; p : in out Poly_Sys ) is

    k,d,n : natural32;
    found : boolean;
    sols : Solution_List;
    np : Poly_Sys(p'first..p'last-1);

  begin
    new_line;
    put_line("Removing last slack variable after going down one level.");
    n := natural32(p'last);
    k := Count_Embed_Symbols(n,"zz");
    Swap_Symbols_to_End(n,k,"zz",p);
    put("dimension of the witness set : "); put(k,1); new_line;
    if k = 0 then
      put_line("no slack variables to remove ...");
    else
      np := Remove_Embedding1(p,1);
      put(outfile,natural32(np'last),np);
      Scan_and_Skip(infile,"SOLUTIONS",found);
      if not found then
        put_line("found no solutions in the input file");
      else
        get(infile,sols);
        d := Length_Of(sols);
        put("degree of the witness set : "); put(d,1); new_line;
        Remove_Component(sols);
        new_line(outfile);
        put_line(outfile,"THE SOLUTIONS :");
        put(outfile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end Remove_Last_Slack_Variable;

end Jumpstart_Diagonal_Homotopies;
