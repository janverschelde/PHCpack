with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Random_Vectors;            use DoblDobl_Random_Vectors;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with DoblDobl_Solution_Diagnostics;      use DoblDobl_Solution_Diagnostics;
with DoblDobl_Select_Solutions;
with DoblDobl_Condition_Tables;          use DoblDobl_Condition_Tables;

package body DoblDobl_Condition_Report is

  procedure Write_Diagnostics ( sols : in Solution_List ) is

    tmp : Solution_List := sols;
    ls : Link_to_Solution;
 
  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      put("  err : "); put(ls.err,3);
      put("  rco : "); put(ls.rco,3);
      put("  res : "); put(ls.res,3); new_line;
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Diagnostics;

  procedure Compute_Condition_Tables ( sols : in Solution_List ) is

  -- NOTE : this code goes three times through the same list...

    e : Standard_Natural_Vectors.Vector(0..30) := Create(30);
    c : Standard_Natural_Vectors.Vector(0..30) := Create(30);
    r : Standard_Natural_Vectors.Vector(0..30) := Create(30);

  begin
    new_line;
    Corrector_Table(e,sols); Write_Corrector_Table(standard_output,e);
    Condition_Table(c,sols); Write_Condition_Table(standard_output,c);
    Residuals_Table(r,sols); Write_Residuals_Table(standard_output,r);
  end Compute_Condition_Tables;

  procedure Write_Condition_Results
               ( file : file_type;
                 e,c,r : in Standard_Natural_Vectors.Vector;
                 cnt_real : in natural32; tol_real : in double_float ) is
  begin
    Write_Corrector_Table(file,e);
    Write_Condition_Table(file,c);
    Write_Residuals_Table(file,r);
    put(file,"Number of real solutions (with tolerance"); 
    put(file,tol_real,3); put(file,") : "); put(file,cnt_real,1);
    put_line(file,".");
  end Write_Condition_Results;

  procedure Write_Condition_Results
               ( file : file_type; i,f : in natural32;
                 e,c,r : in Standard_Natural_Vectors.Vector;
                 cnt_real : in natural32; tol_real : in double_float ) is
  begin
    if i >= f then   -- clear out the frequency updating
      new_line;
    end if;
    new_line;
    Write_Condition_Results(file,e,c,r,cnt_real,tol_real);
  end Write_Condition_Results;

  procedure Count_Clusters
               ( infile : in out file_type; outfile : in file_type;
                 bannered,to_file : in boolean;
                 root : in Link_to_Quad_Node; tol : in double_float ) is

    cnt : natural32 := 0;
    ind : integer32 := 0;
    rv : Standard_Natural_Vectors.Link_to_Vector;
    resfile : file_type;
    ans : character;
    file_created : boolean := false;
    to_vector : boolean := false;

    procedure Count ( lp1,lp2 : in Link_to_Point ) is
    begin
      cnt := cnt + 1;
    end Count;
    procedure Search_Clusters is new DoblDobl_Quad_Trees.Clusters(Count);

    procedure Write ( lp1,lp2 : in Link_to_Point ) is
    begin
      if to_file then
        put(outfile,"  "); put(outfile,lp1.label,1);
        put(outfile,"  "); put(outfile,lp2.label,1); new_line(outfile);
      elsif file_created then
        put(resfile,"  "); put(resfile,lp1.label,1);
        put(resfile,"  "); put(resfile,lp2.label,1); new_line(resfile);
      elsif to_vector then
        ind := ind + 1; rv(ind) := natural32(lp1.label);
        ind := ind + 1; rv(ind) := natural32(lp2.label);
      else
        put("  "); put(lp1.label,1);
        put("  "); put(lp2.label,1); new_line;
      end if;
    end Write;
    procedure Write_Clusters is new DoblDobl_Quad_Trees.Clusters(Write);

  begin
    Search_Clusters(root,tol);
    if cnt = 0 then
      put_line(outfile,"  found no pairs of clustered solutions.");
    elsif cnt = 1 then
      put_line(outfile,"  found one pair of candidate clustered solutions.");
    else
      put(outfile,"  found "); put(outfile,cnt,1);
      put_line(outfile," pairs of candidate clustered solutions.");
    end if;
    if cnt > 0 then
      if to_file then
        put_line(outfile,"The pairs of candidate clustered solutions :");
        Write_Clusters(root,tol);
      else
        new_line;
        put_line("MENU to report candidate clustered pairs of solutions.");
        put_line("  0. do not report which pairs are candidate clustered;");
        put_line
          ("  1. write pairs of candidate clustered solutions to screen;");
        put_line
          ("  2. write pairs of candidate clustered solutions to file;");
        put_line
          ("  3. write the entire candidate clustered solutions to file.");
        put("Type 0, 1, 2, or 3 to select an action : ");
        Ask_Alternative(ans,"0123");
        if ans = '1' then
          new_line;
          Write_Clusters(root,tol);
        elsif ans = '2' or ans = '3' then
          new_line;
          put_line("Reading the name of the output file...");
          Read_Name_and_Create_File(resfile);
          if ans = '2' then
            file_created := true;
            Write_Clusters(root,tol);
          else
            to_vector := true;
            rv := new Standard_Natural_Vectors.Vector(1..integer32(2*cnt));
            Write_Clusters(root,tol);
            new_line;
            put("the indices : "); put(rv); new_line;
            Reset(infile);
            declare
              sa : Solution_Array(rv'range);
              sv : Standard_Natural_Vectors.Vector(rv'range);
              fail : boolean;
            begin
              DoblDobl_Select_Solutions.Select_from_File
                (infile,bannered,rv.all,sv,sa,fail);
              if not fail then
                DoblDobl_Select_Solutions.Write_Selection
                  (resfile,natural32(sa(sa'first).n),rv.all,sv,sa);
              end if;
            end;
          end if;
        end if;
      end if;
    end if;
  end Count_Clusters;

  procedure Count_Clusters
               ( outfile : in file_type; to_file : in boolean;
                 sols : in Solution_List;
                 root : in Link_to_Quad_Node; tol : in double_float ) is

    cnt : natural32 := 0;
    ind : integer32 := 0;
    rv : Standard_Natural_Vectors.Link_to_Vector;
    resfile : file_type;
    ans : character;
    file_created : boolean := false;
    to_vector : boolean := false;

    procedure Count ( lp1,lp2 : in Link_to_Point ) is
    begin
      cnt := cnt + 1;
    end Count;
    procedure Search_Clusters is new DoblDobl_Quad_Trees.Clusters(Count);

    procedure Write ( lp1,lp2 : in Link_to_Point ) is
    begin
      if to_file then
        put(outfile,"  "); put(outfile,lp1.label,1);
        put(outfile,"  "); put(outfile,lp2.label,1); new_line(outfile);
      elsif file_created then
        put(resfile,"  "); put(resfile,lp1.label,1);
        put(resfile,"  "); put(resfile,lp2.label,1); new_line(resfile);
      elsif to_vector then
        ind := ind + 1; rv(ind) := natural32(lp1.label);
        ind := ind + 1; rv(ind) := natural32(lp2.label);
      else
        put("  "); put(lp1.label,1);
        put("  "); put(lp2.label,1); new_line;
      end if;
    end Write;
    procedure Write_Clusters is new DoblDobl_Quad_Trees.Clusters(Write);

  begin
    Search_Clusters(root,tol);
    if cnt = 0 then
      put_line(outfile,"  found no pairs of clustered solutions.");
    elsif cnt = 1 then
      put_line(outfile,"  found one pair of candidate clustered solutions.");
    else
      put(outfile,"  found "); put(outfile,cnt,1);
      put_line(outfile," pairs of candidate clustered solutions.");
    end if;
    if cnt > 0 then
      if to_file then
        put_line(outfile,"The pairs of candidate clustered solutions :");
        Write_Clusters(root,tol);
      else
        new_line;
        put_line("MENU to report candidate clustered pairs of solutions.");
        put_line("  0. do not report which pairs are candidate clustered;");
        put_line
          ("  1. write pairs of candidate clustered solutions to screen;");
        put_line
          ("  2. write pairs of candidate clustered solutions to file;");
        put_line
          ("  3. write the entire candidate clustered solutions to file.");
        put("Type 0, 1, 2, or 3 to select an action : ");
        Ask_Alternative(ans,"0123");
        if ans = '1' then
          new_line;
          Write_Clusters(root,tol);
        elsif ans = '2' or ans = '3' then
          new_line;
          put_line("Reading the name of the output file...");
          Read_Name_and_Create_File(resfile);
          if ans = '2' then
            file_created := true;
            Write_Clusters(root,tol);
          else
            to_vector := true;
            rv := new Standard_Natural_Vectors.Vector(1..integer32(2*cnt));
            Write_Clusters(root,tol);
            new_line;
            put("the indices : "); put(rv); new_line;
            declare
              sa : Solution_Array(rv'range);
              sv : Standard_Natural_Vectors.Vector(rv'range);
            begin
              DoblDobl_Select_Solutions.Select_Solutions(sols,rv.all,sv,sa);
              DoblDobl_Select_Solutions.Write_Selection
                (resfile,natural32(sa(sa'first).n),rv.all,sv,sa);
            end;
          end if;
        end if;
      end if;
    end if;
  end Count_Clusters;

  procedure Write_Cluster_Report
               ( infile : in out file_type; outfile : in file_type;
                 bannered,to_file : in boolean;
                 pl : in Point_List; size : in natural32;
                 tol : in double_float ) is

    root : Link_to_Quad_Node := Create_Root_Leaf(pl,size);
    max_depth : constant natural32 := 10;
    min_size : constant natural32 := 1024;

  begin
    put_line(outfile,"Creating a quadtree to search for clusters ...");
    Create(root,max_depth,min_size);
    Sort_Leaves(root);
    Count_Clusters(infile,outfile,bannered,to_file,root,tol);
  end Write_Cluster_Report;

  procedure Write_Cluster_Report
               ( outfile : in file_type; to_file : in boolean;
                 sols : in Solution_List;
                 pl : in Point_List; tol : in double_float ) is

    size : constant natural32 := Length_Of(sols);
    root : Link_to_Quad_Node := Create_Root_Leaf(pl,size);
    max_depth : constant natural32 := 10;
    min_size : constant natural32 := 1024;

  begin
    put_line(outfile,"Creating a quadtree to search for clusters ...");
    Create(root,max_depth,min_size);
    Sort_Leaves(root);
    Count_Clusters(outfile,to_file,sols,root,tol);
  end Write_Cluster_Report;

  procedure Is_Clustered
               ( s : in Solution; nb : in natural32;
                 sols : in Solution_List; tol : in double_float; 
                 h1,h2 : in DoblDobl_Complex_Vectors.Vector;
                 pl : in out Point_List; val : out natural32 ) is

    pt : constant Point := Create(s.v,h1,h2,integer32(nb));
    lpt : constant Link_to_Point := new Point'(pt);
    lbl : integer32;
    ls : Link_to_Solution;

  begin
    Insert_no_Duplicates(pl,lpt,tol,lbl);
    if lbl = lpt.label then
      val := nb;
    else
      ls := DoblDobl_Complex_Solutions.Retrieve(sols,natural32(lbl));
      if DoblDobl_Solution_Diagnostics.Equal(ls.all,s,tol)
       then val := natural32(lbl);
       else val := nb;
      end if;
    end if;
  end Is_Clustered;

  procedure Is_Clustered
               ( s : in Solution; nb : in natural32;
                 sols : in Solution_Array; tol : in double_float; 
                 h1,h2 : in DoblDobl_Complex_Vectors.Vector;
                 pl : in out Point_List; val : out natural32 ) is

    pt : constant Point := Create(s.v,h1,h2,integer32(nb));
    lpt : constant Link_to_Point := new Point'(pt);
    lbl : integer32;

  begin
    Insert_no_Duplicates(pl,lpt,tol,lbl);
    if lbl = lpt.label then
      val := nb;
    else
      if DoblDobl_Solution_Diagnostics.Equal(sols(lbl).all,s,tol)
       then val := natural32(lbl);
       else val := nb;
      end if;
    end if;
  end Is_Clustered;

  procedure Multiplicity
               ( s : in out Solution; nb : in natural32;
                 sols : in Solution_List; tol : in double_float; 
                 h1,h2 : in DoblDobl_Complex_Vectors.Vector;
                 pl : in out Point_List; val : out natural32 ) is

    pt : constant Point := Create(s.v,h1,h2,integer32(nb));
    lpt : constant Link_to_Point := new Point'(pt);
    cnt,lbl : integer32;
    ptpos : Point_List;
    ls : Link_to_Solution;

  begin
    Insert_with_Duplicates(pl,lpt,tol,cnt,ptpos);
    val := 1;
    if cnt > 1 then
      for i in 1..cnt loop
        ptpos := Tail_Of(ptpos);
        exit when Is_Null(ptpos);
        lbl := Head_Of(ptpos).label;
        ls := DoblDobl_Complex_Solutions.Retrieve(sols,natural32(lbl));
        if DoblDobl_Solution_Diagnostics.Equal(ls.all,s,tol) then
          ls.m := ls.m + 1;
          val := val + 1;
        end if;
      end loop;
      s.m := integer32(val);
    end if;
  end Multiplicity;

  procedure Multiplicity
               ( s : in out Solution; nb : in natural32;
                 sols : in Solution_Array; tol : in double_float; 
                 h1,h2 : in DoblDobl_Complex_Vectors.Vector;
                 pl : in out Point_List; val : out natural32 ) is

    pt : constant Point := Create(s.v,h1,h2,integer32(nb));
    lpt : constant Link_to_Point := new Point'(pt);
    cnt,lbl : integer32;
    ptpos : Point_List;

  begin
    Insert_with_Duplicates(pl,lpt,tol,cnt,ptpos);
    val := 1;
    if cnt > 1 then
      for i in 1..cnt loop
        ptpos := Tail_Of(ptpos);
        exit when Is_Null(ptpos);
        lbl := Head_Of(ptpos).label;
        if DoblDobl_Solution_Diagnostics.Equal(sols(lbl).all,s,tol) then
          sols(lbl).m := sols(lbl).m + 1;
          val := val + 1;
        end if;
      end loop;
      s.m := integer32(val);
    end if;
  end Multiplicity;

  procedure Scan_for_Condition_Tables 
               ( infile : in out file_type; outfile : in file_type; 
                 bannered,to_file : in boolean;
                 len,dim : in natural32;
                 tol_real,tol_clus : in double_float;
                 i_end,f_val : out natural32;
                 e,c,r : out Standard_Natural_Vectors.Vector;
                 nb_real : out natural32; pl : out Point_List ) is

    s : Solution(integer32(dim));
    i : natural32 := 1;
    f : natural32 := 1024;  -- frequency updater
    cnt_real : natural32 := 0;
    h1 : constant DoblDobl_Complex_Vectors.Vector(1..integer32(dim))
       := Random_Vector(1,integer32(dim));
    h2 : constant DoblDobl_Complex_Vectors.Vector(1..integer32(dim))
       := Random_Vector(1,integer32(dim));
    first,last : Point_List;

  begin
    e := DoblDobl_Condition_Tables.Create(30);
    c := DoblDobl_Condition_Tables.Create(30);
    r := DoblDobl_Condition_Tables.Create(30);
    new_line;
    put("Scanning the input file for "); put(len,1); 
    put(" solutions of dimension "); put(dim,1); put_line("...");
    while i <= len loop
      Read_Next(infile,s);
      Update_Corrector(e,s.err);
      Update_Condition(c,s.rco);
      Update_Residuals(r,s.res);
      if Is_Real(s,tol_real)
       then cnt_real := cnt_real + 1;
      end if;
      Append(first,last,h1,h2,integer32(i),s.v);
      i := i+1;
      if i mod f = 0 
       then put(i,1); put(" ... "); f := 2*f;
      end if;
    end loop;
    nb_real := cnt_real; pl := first; i_end := i; f_val := f;
  exception
    when others => 
      put("An error occurred while reading solution ");
      put(i,1); put_line(".");
      nb_real := cnt_real; pl := first; i_end := i; f_val := f;
  end Scan_for_Condition_Tables;

  procedure Scan_for_Condition_Tables 
               ( file : in file_type; sols : in Solution_List;
                 tol_real,tol_clus : in double_float;
                 e,c,r : out Standard_Natural_Vectors.Vector;
                 nb_real : out natural32; pl : out Point_List ) is

    len : constant natural32 := Length_Of(sols);
    tmp : Solution_List := sols;
    ls : Link_to_Solution := Head_Of(sols);
    dim : constant integer32 := ls.n;
    i : natural32 := 1;
    cnt_real : natural32 := 0;
    h1 : constant DoblDobl_Complex_Vectors.Vector(1..dim)
       := Random_Vector(1,dim);
    h2 : constant DoblDobl_Complex_Vectors.Vector(1..dim)
       := Random_Vector(1,dim);
    first,last : Point_List;

  begin
    e := DoblDobl_Condition_Tables.Create(30);
    c := DoblDobl_Condition_Tables.Create(30);
    r := DoblDobl_Condition_Tables.Create(30);
    while i <= len loop
      ls := Head_Of(tmp);
      Update_Corrector(e,ls.err);
      Update_Condition(c,ls.rco);
      Update_Residuals(r,ls.res);
      if Is_Real(ls.all,tol_real)
       then cnt_real := cnt_real + 1;
      end if;
      Append(first,last,h1,h2,integer32(i),ls.v);
      tmp := Tail_Of(tmp); i := i + 1;
    end loop;
    nb_real := cnt_real; pl := first;
  exception
    when others => 
      put("An error occurred while processing solution ");
      put(i,1); put_line("."); pl := first;
  end Scan_for_Condition_Tables;

end DoblDobl_Condition_Report;
