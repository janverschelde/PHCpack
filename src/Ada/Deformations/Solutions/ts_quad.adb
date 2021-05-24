with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with File_Scanning;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Natural_Matrices;
with Standard_Natural_Matrices_io;       use Standard_Natural_Matrices_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_Point_Lists;               use Standard_Point_Lists;
with Standard_Quad_Trees;                use Standard_Quad_Trees;
with Standard_Quad_Trees_io;             use Standard_Quad_Trees_io;
with Standard_Condition_Report;          use Standard_Condition_Report;
with Standard_Select_Solutions;

procedure ts_quad is
 
-- DESCRIPTION :
--   Interactive testing facility on the package Standard_Quad_Trees.

-- FREQUENCY MATRIX :
--   The first hypothesis was that, thanks to random projections,
--   the points would also be distributed at random in the plane
--   and we could just use a fixed evenly spaced grid to distribute
--   the points in buckets.  Well, this hypothesis is false...

  function Classify ( a,dx,x : double_float ) return integer32 is

  -- DESCRIPTION :
  --   Returns 0 if a < x, otherwise returns the index k such that x lies 
  --   between a+(k-1)*dx and a+k*dx.

    fk : double_float;
    k : integer32;

  begin
    if x < a then
      return 0;
    else
      fk := (x-a)/dx;
       k := integer32(fk);
       if double_float(k) > fk
        then return k;
        else return k+1;
       end if;
    end if;
  end Classify;

  procedure Fit ( k : in out integer32; n : in integer32 ) is

  -- DESCRIPTION :
  --   If k = 0 then k is set to 1, if k > n, then k is set to n.

  begin
    if k = 0 then
      k := 1;
    elsif k > n then
      k := n;
    end if;
  end Fit;      

  procedure Frequency_Matrix
              ( pl : in Point_List;
                a,b,c,d : in double_float; n : in natural32 ) is

  -- DESCRIPTION :
  --   Creates a frequency matrix for the points in the list,
  --   using a [a,b]x[c,d] bounding box and n intervals for both x and y.

  -- ON ENTRY :
  --   pl       list of points in the plane;
  --   a        lower bound for the x-coordinate of the points;
  --   b        upper bound for the x-coordinate of the points;
  --   c        lower bound for the y-coordinate of the points;
  --   d        upper bound for the y-coordinate of the points;
  --   n        number of intervals to be used for both x and y.

    timer : Timing_Widget;
    fm : Standard_Natural_Matrices.Matrix(1..integer32(n),1..integer32(n));
    dx : constant double_float := (b-a)/double_float(n);
    dy : constant double_float := (d-c)/double_float(n);
    tmp : Point_List := pl;
    lp : Link_to_Point;
    ind_x,ind_y : integer32;
    sum : natural32;

  begin
    for i in fm'range(1) loop
      for j in fm'range(2) loop
        fm(i,j) := 0;
      end loop;
    end loop;
    tstart(timer);
    while not Is_Null(tmp) loop
      lp := Head_Of(tmp);
      ind_x := Classify(a,dx,lp.x); Fit(ind_x,integer32(n));
      ind_y := Classify(c,dy,lp.y); Fit(ind_y,integer32(n));
      fm(ind_x,ind_y) := fm(ind_x,ind_y) + 1;
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
    new_line;
    put_line("The frequency matrix :"); put(fm);
    sum := 0; -- sanity check ...
    for i in fm'range(1) loop
      for j in fm'range(2) loop
        sum := sum + fm(i,j);
      end loop;
    end loop;
    put("#elements counted in the matrix : "); put(sum,1);
    new_line; new_line;
    print_times(standard_output,timer,"frequency matrix");
  end Frequency_Matrix;

  procedure Bounding_Box ( pl : in Point_List ) is

  -- DESCRIPTION :
  --   Computes a bounding box around the points in the list,
  --   defined by minimal and maximal x and y coordinates.

    timer : Timing_Widget;
    tmp : Point_List := Tail_Of(pl);
    lp : Link_to_Point := Head_Of(pl);
    min_x : double_float := lp.x;
    max_x : double_float := lp.x;
    min_y : double_float := lp.y;
    max_y : double_float := lp.y;
    cnt : natural32 := 1;

  begin
    tstart(timer);
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      lp := Head_Of(tmp);
      if lp.x < min_x then min_x := lp.x; end if;
      if lp.x > max_x then max_x := lp.x; end if;
      if lp.y < min_y then min_y := lp.y; end if;
      if lp.y > max_y then max_y := lp.y; end if;
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
    put("Counted "); put(cnt,1); put_line(" points in the plane.");
    put("Bounding box : ");
    put("["); put(min_x,2); put(","); put(max_x,2); put("]");
    put("x");
    put("["); put(min_y,2); put(","); put(max_y,2); put("]");
    new_line;
    new_line;
    print_times(standard_output,timer,"computing a bounding box");
    Frequency_Matrix(pl,min_x,max_x,min_y,max_y,32);
  end Bounding_Box;

-- end of FREQUENCY MATRIX code

  procedure Quadrant_Partition ( pl : in Point_List ) is

    timer : Timing_Widget;
    cx,cy : double_float;
    ne_cnt,nw_cnt,sw_cnt,se_cnt : natural32;
    ne_pts,nw_pts,sw_pts,se_pts : Point_List;

  begin
    tstart(timer);
    Center(pl,cx,cy);
    tstop(timer);
    new_line;
    put("Center : ("); put(cx); put(","); put(cy); put_line(").");
    new_line;
    print_times(standard_output,timer,"computing the center");
    tstart(timer);
    Partition(pl,cx,cy,ne_cnt,nw_cnt,sw_cnt,se_cnt,
                       ne_pts,nw_pts,sw_pts,se_pts);
    tstop(timer);
    new_line;
    put_line("Distribution of the points into quadrants : ");
    put("  #points in north east quadrant : "); put(ne_cnt,1); new_line;
    put("  #points in north west quadrant : "); put(nw_cnt,1); new_line;
    put("  #points in south west quadrant : "); put(sw_cnt,1); new_line;
    put("  #points in south east quadrant : "); put(se_cnt,1); new_line;
    new_line;
    print_times(standard_output,timer,"partitioning into quadrants");
  end Quadrant_Partition;

  procedure Hash_Solutions
              ( infile : in file_type; len,dim : in natural32;
                pl : out Point_List ) is

    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    s : Solution(integer32(dim));
    h1 : constant Vector(1..integer32(dim)) := Random_Vector(1,integer32(dim));
    h2 : constant Vector(1..integer32(dim)) := Random_Vector(1,integer32(dim));
    freq : natural32 := 1024; -- frequency updater
    first,last : Point_List;

  begin
    put("Hashing solutions ... ");
    tstart(timer);
    for i in 1..len loop
      Read_Next(infile,s);
      if i mod freq = 0
       then put(i,1); put(" ... "); freq := 2*freq;
      end if;
      Append(first,last,h1,h2,integer32(i),s.v);
    end loop;
    tstop(timer);
    new_line;
    new_line;
    print_times(standard_output,timer,"creation of point list");
    pl := first;
  exception when others =>
    put_line("An exception was raised before reached end...");
    raise;
  end Hash_Solutions;

  procedure Read_File_and_Hash_Solutions 
              ( file : in out file_type; bannered : out boolean;
                pl : out Point_List; length : out natural32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a file name for a solution list
  --   and returns a list of points in the plane.

  -- ON RETURN :
  --   file     file type of the input file;
  --   bannered is true if the solutions were preceeded by a banner;
  --   pl       list of points in the plane, obtained by application
  --            of hash functions on the list of solutions;
  --   length   size of the list pl, equals number of solutions.

    use Standard_Select_Solutions;

    dim : natural32;
    fail : boolean;

  begin
    put_line("Reading the name of the input file for the solutions...");
    Read_Name_and_Open_File(file);
    Scan_Banner_Dimensions(file,length,dim,bannered,fail);
    if fail then
      put("Format of the solution list on file is incorrect.");
      put_line("  Please try again.");
    else
      put("Ready to hash "); put(length,1);
      put(" solutions of dimension "); put(dim,1); put_line(".");
      new_line;
      Hash_Solutions(file,length,dim,pl);
    end if;
  end Read_File_and_Hash_Solutions;

  procedure Cluster_Report
              ( infile : in out file_type; bannered : in boolean;
                root : in Link_to_Quad_Node; tol : in double_float ) is

  -- DESCRIPTION :
  --   Scans all sorted leaves of the quad node and reports all pairs
  --   of clustered points.

    use Standard_Complex_Solutions;

    cnt : natural32 := 0;
    ind : integer32 := 0;
    rv : Standard_Natural_Vectors.Link_to_Vector;
    sols : Solution_List;
    ans : character;
    to_file : boolean := false;
    to_vector : boolean := false;
    outfile : file_type;

    procedure Report_Cluster ( lp1,lp2 : in Link_to_Point ) is
    begin
      if to_file then
        if cnt = 0 then
          put_line(outfile,"Labels of pairs of clustered solutions : ");
        end if;
        put(outfile,"  "); put(outfile,lp1.label,1);
        put(outfile,"  "); put(outfile,lp2.label,1);
        new_line(outfile);
      elsif to_vector then
        ind := ind + 1; rv(ind) := natural32(lp1.label);
        ind := ind + 1; rv(ind) := natural32(lp2.label);
      else
        if cnt = 0 then
          put_line("Labels of pairs of clustered solutions : ");
        end if;
        put("  "); put(lp1.label,1);
        put("  "); put(lp2.label,1);
        new_line;
      end if;
      cnt := cnt + 1;
    end Report_Cluster;
    procedure Report_Clusters is
      new Standard_Quad_Trees.Clusters(Report_Cluster);

  begin
    Report_Clusters(root,tol);
    if cnt = 0 then
      put_line("Found no pairs of clustered solutions.");
    elsif cnt = 1 then
      put_line("Found one pair of clustered solutions.");
    else
      put("Found "); put(cnt,1);
      put_line(" pairs of clustered solutions.");
    end if;
    if cnt > 0 then
      put("Do you wish to write the pairs to file ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("MENU to write clustered pairs :");
        put_line("  1. write indices of the clustered pairs to file; or");
        put_line("  2. write solutions to file.");
        put("Type 1 or 2 to choose : "); Ask_Alternative(ans,"12");
        new_line;
        put_line("Reading the name of the output file...");
        Read_Name_and_Create_File(outfile);
        if ans = '1' then
          to_file := true; cnt := 0;
          Report_Clusters(root,tol);
        else
          to_vector := true;
          rv := new Standard_Natural_Vectors.Vector(1..integer32(2*cnt));
          Report_Clusters(root,tol);
          new_line;
          put("The indices : "); put(rv); new_line;
          Reset(infile);
          Standard_Select_Solutions.Select_from_File
            (infile,bannered,rv.all,sols);
          put(outfile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
        end if;
      end if;
    end if;
  end Cluster_Report;

  procedure Create_Quad_Tree
              ( file : in out file_type; bannered : in boolean;
                pl : in Point_List; size : in natural32 ) is

    timer : Timing_Widget;
    root : Link_to_Quad_Node := Create_Root_Leaf(pl,size);
    max_depth,min_size : natural32 := 0;
    ans : character;

  begin
    put("Give the maximal depth of the tree : "); get(max_depth);
    put("Give the minimal splitting size : "); get(min_size);
    tstart(timer);
    Create(root,max_depth,min_size);
    tstop(timer);
    new_line;
    put_line("The tree : "); put(root);
    new_line;
    put("The number of leaves in the tree : ");
    put(Number_of_Leaves(root),1); new_line;
    put("The number of nodes in the tree : ");
    put(Number_of_Nodes(root),1); new_line;
    new_line;
    print_times(standard_output,timer,"creation of a quad tree");
    new_line;
    put("Do you wish to sort the leaves ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      tstart(timer);
      Sort_Leaves(root);
      tstop(timer);
      new_line;
      print_times(standard_output,timer,"sorting the leaves");
      new_line;
      tstart(timer);
      Cluster_Report(file,bannered,root,1.0E-8);
      tstop(timer);
      new_line;
      print_times(standard_output,timer,"reporting the clusters");
    end if;
  end Create_Quad_Tree;

  procedure In_Memory_Solutions is

  -- DESCRIPTION :
  --   This procedure test the quad tree construction for a solution list
  --   entirely loaded into main memory.

    use Standard_Complex_Solutions;

    infile : file_type;
    ans : character;
    sols : Solution_List;
    fail,found : boolean;
    e,c,r : Standard_Natural_Vectors.Vector(0..15);
    nbreal : natural32;
    pl : Point_List;

  begin
    new_line;
    put_line("Reading solution list first into memory...");
    Read_Name_and_Open_File(infile);
    new_line;
    put("Is the solution list preceded by a system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'n' then
      get(infile,sols); fail := false;
    else
      put("Scanning for THE SOLUTIONS banner ...");
      File_Scanning.Scan_and_Skip(infile,"THE SOLUTIONS",found);
      if found then
        put_line(" found.");
        get(infile,sols);
        fail := false;
      else 
        put_line(" not found!");
        fail := true;
      end if;
    end if;
    if not fail then
      put("Read "); put(Length_Of(sols),1); put(" solutions");
      put(" of dimension "); put(Head_Of(sols).n,1); put_line(".");
      Scan_for_Condition_Tables
        (standard_output,sols,1.0E-8,1.0E-8,e,c,r,nbreal,pl);
      Write_Condition_Results(standard_output,e,c,r,nbreal,1.0E-8);
      Write_Cluster_Report(standard_output,false,sols,pl,1.0E-8);
    end if;
  end In_Memory_Solutions;

  procedure Main is

    file : file_type;
    bannered : boolean;
    pl : Point_List;
    sz : natural32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Creating quad trees for lists of solutions.");
    loop
      new_line;
      put_line("MENU for operations on the list of points : ");
      put_line("  0. leave this menu;");
      put_line("  1. compute a bounding box and frequency matrix;");
      put_line("  2. apply once a quadrant partition to the list;");
      put_line("  3. create a quad tree for the solution list on file;");
      put_line("  4. make quad tree after reading ALL solutions first.");
      put("Type 0, 1, 2, 3, or 4 to select an operation : ");
      Ask_Alternative(ans,"01234");
      exit when ans = '0';
      if ans = '4' then
        In_Memory_Solutions;
      else
        new_line;
        Read_File_and_Hash_Solutions(file,bannered,pl,sz);
        new_line;
        if ans = '1' then
          Bounding_Box(pl);
        elsif ans = '2' then
          Quadrant_Partition(pl);
        else
          Create_Quad_Tree(file,bannered,pl,sz);
        end if;
      end if;
    end loop;
  end Main;

begin
  Main;
end ts_quad;
