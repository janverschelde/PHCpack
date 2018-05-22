with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Lists_of_Integer_Vectors;
with Lists_of_Floating_Vectors;
with Arrays_of_Floating_Vector_Lists;   use Arrays_of_Floating_Vector_Lists;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Standard_Integer_VecVecs;
with Standard_Integer64_Matrices;
with Standard_Integer64_Linear_Solvers; use Standard_Integer64_Linear_Solvers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Standard_Poly_Laur_Convertors;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_Systems_io;  use DoblDobl_Complex_Laur_Systems_io;
with DoblDobl_Polynomial_Convertors;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_Systems_io;  use QuadDobl_Complex_Laur_Systems_io;
with QuadDobl_Polynomial_Convertors;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;     use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;     use QuadDobl_Complex_Solutions_io;
with Floating_Lifting_Functions;
with Floating_Mixed_Subdivisions;       use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;    use Floating_Mixed_Subdivisions_io;
with Cell_Stack;                        use Cell_Stack;
with Mixed_Volume;
with MixedVol_Algorithm;                use MixedVol_Algorithm;
with Polyhedral_Start_Systems;
with Mixed_Labels_Queue;
with Semaphore;
with Multitasking;
with Pipelined_Labeled_Cells;           use Pipelined_Labeled_Cells;
with Pipelined_Polyhedral_Trackers;     use Pipelined_Polyhedral_Trackers;
with Pipelined_Polyhedral_Drivers;

procedure ts_mtmva is

-- DESCRIPTION :
--   Running the MixedVol algorithm with a callback function
--   from within a single task.

  procedure Write_Mixed_Cells 
              ( file : in file_type; nbequ,r : in integer32;
                mtype : in Standard_Integer_Vectors.Link_to_Vector;
                mcc : in Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Writes the mixed cells configuration to file.

  -- ON ENTRY :
  --   file     to write the mixed cell configuration on;
  --   nbequ    number of equations and number of variables;
  --   r        number of different supports
  --   mtype    type of mixture of the supports;
  --   mcc      mixed cells.

    mix : constant Standard_Integer_Vectors.Vector := Mixture(r,mtype);

  begin
    put(file,natural32(nbequ),mix,mcc);
  end Write_Mixed_Cells;

  procedure Mixed_Cell_Configuration
              ( nbequ,r,size,nb : in integer32;
                mtype,perm : in Standard_Integer_Vectors.Link_to_Vector;
                Vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                lft : in Standard_Floating_Vectors.Link_to_Vector;
                labels : in Lists_of_Integer_Vectors.List;
                sub : out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Takes the cells in the cell stack and converts the cells into
  --   a mixed cell configuraiton.

  -- ON ENTRY :
  --   nbequ    number of equations in the input Laurent system,
  --            or alternatively, the number of variables;
  --   r        number of different supports
  --   size     the number of labels in a mixed cell;
  --   nb       total number of cells;
  --   mtype    type of mixture of the supports;
  --   perm     permutation on the supports;
  --   Vtx      vertex set of all supports;
  --   lft      lifting values for all points in Vtx;
  --   labels   labels to the coordinates of the mixed cells.

  -- ON RETURN :
  --   sub      a mixed cell configuration.

    use Lists_of_Integer_Vectors;

    last : Mixed_Subdivision;
    tmp : List := labels;
    lbl : Standard_Integer_Vectors.Link_to_Vector;

  begin
    put("Permutation of the supports : "); put(perm); new_line;
    put_line("Creating a regular mixed-cell configuration ...");
    new_line;
    put_line("See the output file for the mixed-cell configuration ...");
    new_line;
    if r < nbequ then
      while not Is_Null(tmp) loop
        lbl := Head_Of(tmp);
        declare
          mic : constant Mixed_Cell
              := Labels_to_Mixed_Cell(nbequ,r,mtype,lbl,Vtx,lft);
        begin
          Append(sub,last,mic);
        end;
        tmp := Tail_Of(tmp);
      end loop;
    else
      while not Is_Null(tmp) loop
        lbl := Head_Of(tmp);
        declare
          mic : constant Mixed_Cell
              := Labels_to_Mixed_Cell(nbequ,r,mtype,perm,lbl,Vtx,lft);
        begin
          Append(sub,last,mic);
        end;
        tmp := Tail_Of(tmp);
      end loop;
    end if;
  end Mixed_Cell_Configuration;

  procedure Mixed_Cell_Configuration
              ( nbequ,r,size,nb : in integer32; stlb : in double_float;
                mtype,perm : in Standard_Integer_Vectors.Link_to_Vector;
                Vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                lft : in Standard_Floating_Vectors.Link_to_Vector;
                labels : in Lists_of_Integer_Vectors.List;
                sub,mcc,stbmcc : out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Takes the cells in the cell stack and converts the cells into
  --   a mixed cell configuraiton, classifying the cells into the
  --   original ones (without artificial origins) and stable ones.

  -- ON ENTRY :
  --   nbequ    number of equations in the input Laurent system,
  --            or alternatively, the number of variables;
  --   r        number of different supports
  --   size     the number of labels in a mixed cell;
  --   nb       total number of cells;
  --   stlb     lifting bound for stable mixed cells;
  --   mtype    type of mixture of the supports;
  --   perm     permutation on the supports;
  --   Vtx      vertex set of all supports;
  --   lft      lifting values for all points in Vtx;
  --   labels   labels to the coordinates of the mixed cells.

  -- ON RETURN :
  --   sub      a mixed cell configuration, all mixed cells;
  --   mcc      mixed cells without artificial origin;
  --   stbmcc   stable mixed cells with artificial origin.

    use Lists_of_Integer_Vectors;

    last,mcclast,stbmcclast : Mixed_Subdivision;
    tmp : List := labels;
    lbl : Standard_Integer_Vectors.Link_to_Vector;

  begin
    put("Permutation of the supports : "); put(perm); new_line;
    put_line("Creating a regular mixed-cell configuration ...");
    new_line;
    put_line("See the output file for the mixed-cell configuration ...");
    new_line;
    if r < nbequ then
      while not Is_Null(tmp) loop
        lbl := Head_Of(tmp);
        declare
          mic : constant Mixed_Cell
              := Labels_to_Mixed_Cell(nbequ,r,mtype,lbl,Vtx,lft);
        begin
          Append(sub,last,mic);
          if Is_Original(mic,stlb) then
            Append(mcc,mcclast,mic);
          elsif Is_Stable(mic.nor.all,stlb,mic.pts.all) then
            Append(stbmcc,stbmcclast,mic);
          end if;
        end;
        tmp := Tail_Of(tmp);
      end loop;
    else
      while not Is_Null(tmp) loop
        lbl := Head_Of(tmp);
        declare
          mic : constant Mixed_Cell
              := Labels_to_Mixed_Cell(nbequ,r,mtype,perm,lbl,Vtx,lft);
        begin
          Append(sub,last,mic);
          if Is_Original(mic,stlb) then
            Append(mcc,mcclast,mic);
          elsif Is_Stable(mic.nor.all,stlb,mic.pts.all) then
            Append(stbmcc,stbmcclast,mic);
          end if;
        end;
        tmp := Tail_Of(tmp);
      end loop;
    end if;
  end Mixed_Cell_Configuration;

  procedure Sequential_Mixed_Volume_Computation
              ( nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                stlb : in double_float; r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                sub,mcc,stbmcc : out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Extracts the supports of the polynomial system and computes
  --   its mixed volume.  The two stages in the sequential version are
  --   (1) producing the labels to the indices in the mixed cells;
  --   (2) processing the labels into a mixed cell configuration.
  --   These two stages are executed one after the other.

  -- ON ENTRY :
  --   nbequ    the number of equations in the input Laurent system;
  --   nbpts    the total number of points in the supports;
  --   ind      ind(k) marks the beginning of the k-th support;
  --   cnt      cnt(k) counts the number of points in the k-th support;
  --   support  vector range 1..nbequ*nbpts with the coordinates of
  --            all points in the supports;
  --   stlb     lifting bound to use for stable mixed volumes,
  --            equals 0.0 if no stable mixed volumes are requested.

  -- ON RETURN :
  --   r        number of distinct supports;
  --   mtype    the type of mixture of the supports;
  --   perm     permutation of the supports;
  --   sub      a mixed cell configuration for a random lifting.

    use Lists_of_Integer_Vectors;

    nb,size : integer32;
    mixvol : natural32;
    idx,sdx,ndx : Standard_Integer_Vectors.Link_to_Vector;
    vtx,spt : Standard_Integer_VecVecs.Link_to_VecVec;
    lft : Standard_Floating_Vectors.Link_to_Vector;
    cells : CellStack;
    cellcnt : natural32 := 0;
    labels,labels_last : List;

    procedure write_labels ( pts : Standard_Integer_Vectors.Link_to_Vector ) is
    begin
      cellcnt := cellcnt + 1;
      put("Cell "); put(cellcnt,1); put(" is spanned by ");
      put(pts.all); new_line;
      Append(labels,labels_last,pts.all); -- the .all is needed for sharing
    end write_labels;

  begin
    mv_upto_pre4mv
      (nbequ,nbpts,ind,cnt,support.all,r,mtype,perm,idx,vtx,sdx,spt,ndx);
    mv_lift(nbequ,stlb,r,idx,vtx,lft);
    size := Mixed_Volume.cell_size(r,mtype);
    Cs_Init(cells,size);
    Mixed_Volume.MixedVol_with_Callback
      (nbequ,r,size,mtype,idx,vtx,lft,nb,cells,mixvol,false,
       write_labels'access);
    put("The mixed volume is "); put(mixvol,1); put_line(".");
    put("There are "); put(nb,1); put_line(" mixed cells.");
    if stlb = 0.0 then
      Mixed_Cell_Configuration(nbequ,r,size,nb,mtype,perm,Vtx,lft,labels,sub);
    else
      Mixed_Cell_Configuration
       (nbequ,r,size,nb,stlb,mtype,perm,Vtx,lft,labels,sub,mcc,stbmcc);
      put("Total number of mixed cells : "); put(Length_Of(sub),1); new_line;
      put("Mixed cells of original supports : ");
      put(Length_Of(mcc),1); new_line;
      put("Stable mixed cells with artificial origins : ");
      put(Length_Of(stbmcc),1); new_line;
    end if;
  end Sequential_Mixed_Volume_Computation;

  function Extract_Exponent_Vectors
              ( dim : integer32; mic : Mixed_Cell )
              return Standard_Integer64_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the matrix with the exponent vectors defined
  --   by the points in the mixed cell mic.
  --   The dimension of the square matrix is determined by dim,
  --   the dimension before the lifting is applied to the points.

    use Lists_of_Floating_Vectors;

    res : Standard_Integer64_Matrices.Matrix(1..dim,1..dim);
    tmp : List;
    first,ptr : Standard_Floating_Vectors.Link_to_Vector;
    row : integer32 := 0;

  begin
    for i in mic.pts'range loop
      first := Head_Of(mic.pts(i));
      tmp := Tail_Of(mic.pts(i));
      while not Is_Null(tmp) loop
        ptr := Head_Of(tmp);
        row := row + 1;
        for col in 1..dim loop
          res(row,col) := integer64(ptr(col)) - integer64(first(col));
        end loop;
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
    return res;
  end Extract_Exponent_Vectors;

  procedure Mixed_Volume_Calculation
              ( file : in file_type; nt : in integer32; stable : in boolean;
                p : in Standard_Complex_Laur_Systems.Laur_Sys ) is

  -- DESCRIPTION :
  --   Extracts the supports of the polynomial system and computes
  --   its mixed volume.  A mixed subdivision is written to file.
 --    The number of tasks equals nt.

    nbequ : constant integer32 := p'last;
    nbpts,r : integer32 := 0;
    stlb : double_float := 0.0;
    cnt,ind : Standard_Integer_Vectors.Vector(1..nbequ);
    sup,mtype,perm : Standard_Integer_Vectors.Link_to_Vector;
    sub,mcc,orgmcc,stbmcc : Mixed_Subdivision;
    mv : natural32;
    ans : character;
    otp : boolean;
    sem : Semaphore.Lock;
    orgcnt,stbcnt : natural32;
    celcnt : natural32 := 0;
    sumvol : natural64 := 0;
    stabmv : natural64 := 0;
    mat : Standard_Integer64_Matrices.Matrix(1..nbequ,1..nbequ);

    procedure Write_Mixed_Volume
                ( idtask,r : in integer32;
                  mtype : in Standard_Integer_Vectors.Link_to_Vector;
                  mic : in out Mixed_Cell ) is

    -- DESCRIPTION :
    --   Writes the mixed volume of the cell mic with type of mixture
    --   defined in r and mtype to screen, to test the callback procedure
    --   in the pipelined production of the mixed cells.
    --   If the stable mixed volume is wanted, then the classification
    --   of the cell is written to screen as well.

      vol,stabvol : natural64;
      stablecell : boolean := false;
      originalcell : boolean := true;

    begin
      mat := Extract_Exponent_Vectors(nbequ,mic);
      Upper_Triangulate(mat);
      vol := Polyhedral_Start_Systems.Volume_of_Diagonal(mat);
      if stable then
        stablecell := Is_Stable(mic.nor.all,stlb,mic.pts.all);
        originalcell := Is_Original(mic,stlb);
        if stablecell
         then stabvol := vol;
         else stabvol := 0;
        end if;
      end if;
      Semaphore.Request(sem);
      celcnt := celcnt + 1;
      if originalcell
       then sumvol := sumvol + vol;
      end if;
      if stable
       then stabmv := stabmv + stabvol;
      end if;
      Semaphore.Release(sem);
      if not stable then
        put_line("the mixed volume of cell " 
                 & Multitasking.to_string(celcnt) & " is "
                 & Multitasking.to_string(integer32(vol)));
      else
        if Is_Original(mic,stlb) then
          put_line("the mixed volume of cell " 
                   & Multitasking.to_string(celcnt) & " is "
                   & Multitasking.to_string(integer32(vol)) & " original");
        elsif stablecell then
          put_line("the mixed volume of cell " 
                   & Multitasking.to_string(celcnt) & " is "
                   & Multitasking.to_string(integer32(vol)) & " stable");
        else
          put_line("the mixed volume of cell " 
                   & Multitasking.to_string(celcnt) & " is "
                   & Multitasking.to_string(integer32(vol)) & " other");
        end if;
      end if;
    end Write_Mixed_Volume;

  begin
    if stable
     then stlb := Floating_Lifting_Functions.Lifting_Bound(p);
    end if;
    Extract_Supports(nbequ,p,nbpts,ind,cnt,sup);
    if nt < 2 then
      Sequential_Mixed_Volume_Computation
        (nbequ,nbpts,ind,cnt,sup,stlb,r,mtype,perm,sub,mcc,stbmcc);
    else
      new_line;
      put("Monitor the progress of the computations ? (y/n) ");
      Ask_Yes_or_No(ans);
      otp := (ans = 'y');
      if otp then
        Pipelined_Mixed_Cells(nt,nbequ,nbpts,otp,ind,cnt,sup,
          r,mtype,perm,mcc,mv,Write_Mixed_Volume'access);
        if stable
         then put("The sum of the volumes of original cells : ");
         else put("The sum of the volumes of all cells : ");
        end if;
        put(sumvol,1); new_line;
        mv := natural32(sumvol);
      else
        Pipelined_Mixed_Cells(nt,nbequ,nbpts,otp,ind,cnt,sup,
          r,mtype,perm,mcc,mv);
      end if;
      if stable then
        put("The total mixed volume : "); put(mv,1); new_line;
        put("The stable mixed volume : "); put(stabmv,1); new_line;
        Split_Original_Cells(mcc,stlb,orgmcc,stbmcc,orgcnt,stbcnt);
        put("Number of cells without artificial origin : ");
        put(Length_Of(orgmcc),1); new_line;
        put("Number of stable mixed cells : ");
        put(Length_Of(stbmcc),1); new_line;
      else
        put("The mixed volume : "); put(mv,1); new_line;
      end if;
    end if;
    Write_Mixed_Cells(file,nbequ,r,mtype,mcc);
  end Mixed_Volume_Calculation;

  procedure Pipelined_Polyhedral_Driver
              ( file : in file_type; nt : in integer32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys ) is

  -- DESCRIPTION :
  --   Extracts the supports of the polynomial system p and computes
  --   its mixed volume and solves a random coefficient system,
  --   in standard double precision. The number of tasks equals nt.
  --   Calls the driver in Pipelined_Polyhedral_Drivers.

    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;
    use Pipelined_Polyhedral_Drivers;

    cfile,qfile : file_type;
    q : Laur_Sys(p'range);
    mv : natural32;
    sols : Solution_List;
    ans : character;
    mfi,rep : boolean;

  begin
    new_line;
    put_line("Reading the file name to write the start system ...");
    Read_Name_and_Create_File(qfile);
    new_line;
    put("Do you want the mixed cell configuration on file ? (y/n) ");
    Ask_Yes_or_No(ans);
    mfi := (ans = 'y');
    if mfi then
      put_line("Reading the file name to write the cell configuration ...");
      Read_Name_and_Create_File(cfile);
    end if;
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    rep := (ans = 'y');
    Pipelined_Polyhedral_Homotopies(file,cfile,qfile,nt,mfi,rep,p,mv,q,sols);
  end Pipelined_Polyhedral_Driver;

  procedure Pipelined_Polyhedral_Driver
              ( file : in file_type; nt : in integer32;
                stable : in boolean; stlb : in double_float;
                p : in Standard_Complex_Laur_Systems.Laur_Sys ) is

  -- DESCRIPTION :
  --   Extracts the supports of the polynomial system p and computes
  --   its mixed volume and solves a random coefficient system,
  --   in standard double precision. The number of tasks equals nt.
  --   Calls the driver in Pipelined_Polyhedral_Drivers.

    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;
    use Pipelined_Polyhedral_Drivers;

    cfile,qfile : file_type;
    q : Laur_Sys(p'range);
    sub,orgmcc,stbmcc : Mixed_Subdivision;
    mv,orgcnt,stbcnt : natural32;
    sols : Solution_List;
    ans : character;
    mfi,rep : boolean;

  begin
    new_line;
    put_line("Reading the file name to write the start system ...");
    Read_Name_and_Create_File(qfile);
    new_line;
    put("Do you want the mixed cell configuration on file ? (y/n) ");
    Ask_Yes_or_No(ans);
    mfi := (ans = 'y');
    if mfi then
      put_line("Reading the file name to write the cell configuration ...");
      Read_Name_and_Create_File(cfile);
    end if;
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    rep := (ans = 'y');
   -- if stable then
   --   Pipelined_Polyhedral_Homotopies
   --     (file,cfile,qfile,nt,mfi,rep,stable,stlb,p,sub,mv,q,sols);
   --   Split_Original_Cells(sub,stlb,orgmcc,stbmcc,orgcnt,stbcnt);
   --   put("Number of cells without artificial origin : ");
   --   put(Length_Of(orgmcc),1); new_line;
   --   put("Number of stable mixed cells : ");
   --   put(Length_Of(stbmcc),1); new_line;
   -- else
      Pipelined_Polyhedral_Homotopies
        (file,cfile,qfile,nt,mfi,rep,p,mv,q,sols);
   -- end if;
  end Pipelined_Polyhedral_Driver;

  procedure Pipelined_Polyhedral_Driver
              ( file : in file_type; nt : in integer32;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys ) is

  -- DESCRIPTION :
  --   Extracts the supports of the polynomial system p and computes
  --   its mixed volume and solves a random coefficient system,
  --   in double double precision. The number of tasks equals nt.
  --   Calls the driver in Pipelined_Polyhedral_Drivers.

    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Solutions;
    use Pipelined_Polyhedral_Drivers;

    cfile,qfile : file_type;
    q : Laur_Sys(p'range);
    mv : natural32;
    sols : Solution_List;
    ans : character;
    mfi,rep : boolean;

  begin
    new_line;
    put_line("Reading the file name to write the start system ...");
    Read_Name_and_Create_File(qfile);
    new_line;
    put("Do you want the mixed cell configuration on file ? (y/n) ");
    Ask_Yes_or_No(ans);
    mfi := (ans = 'y');
    if mfi then
      put_line("Reading the file name to write the cell configuration ...");
      Read_Name_and_Create_File(cfile);
    end if;
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    rep := (ans = 'y');
    Pipelined_Polyhedral_Homotopies(file,cfile,qfile,nt,mfi,rep,p,mv,q,sols);
  end Pipelined_Polyhedral_Driver;

  procedure Pipelined_Polyhedral_Driver
              ( file : in file_type; nt : in integer32;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys ) is

  -- DESCRIPTION :
  --   Extracts the supports of the polynomial system p and computes
  --   its mixed volume and solves a random coefficient system,
  --   in quad double precision. The number of tasks equals nt.
  --   Calls the driver in Pipelined_Polyhedral_Drivers.

    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Solutions;
    use Pipelined_Polyhedral_Drivers;

    cfile,qfile : file_type;
    q : Laur_Sys(p'range);
    mv : natural32;
    sols : Solution_List;
    ans : character;
    mfi,rep : boolean;

  begin
    new_line;
    put_line("Reading the file name to write the start system ...");
    Read_Name_and_Create_File(qfile);
    new_line;
    put("Do you want the mixed cell configuration on file ? (y/n) ");
    Ask_Yes_or_No(ans);
    mfi := (ans = 'y');
    if mfi then
      put_line("Reading the file name to write the cell configuration ...");
      Read_Name_and_Create_File(cfile);
    end if;
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    rep := (ans = 'y');
    Pipelined_Polyhedral_Homotopies(file,cfile,qfile,nt,mfi,rep,p,mv,q,sols);
  end Pipelined_Polyhedral_Driver;

  procedure Random_Coefficient_System
              ( file : in file_type; nt : in integer32; stable : in boolean;
                p : in Standard_Complex_Laur_Systems.Laur_Sys ) is

  -- DESCRIPTION :
  --   Extracts the supports of the polynomial system p and computes
  --   its mixed volume and solves a random coefficient system,
  --   in standard double precision.  The number of tasks equals nt.

    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

    q : Laur_Sys(p'range);
    nbequ : constant integer32 := p'last;
    nbpts,r : integer32 := 0;
    stlb : double_float := 0.0;
    cnt,ind : Standard_Integer_Vectors.Vector(1..nbequ);
    sup,mtype,perm : Standard_Integer_Vectors.Link_to_Vector;
    lif : Link_to_Array_of_Lists;
    mcc : Mixed_Subdivision;
    mv : natural32;
    sols : Solution_List;
    ans : character;

  begin
    if stable
     then stlb := Floating_Lifting_Functions.Lifting_Bound(p);
    end if;
    put("Test the pipelined driver ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      if stable
       then Pipelined_Polyhedral_Driver(file,nt,stable,stlb,p);
       else Pipelined_Polyhedral_Driver(file,nt,p);
      end if;
      return;
    end if;
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    Extract_Supports(nbequ,p,nbpts,ind,cnt,sup);
    if ans = 'y' then
      Reporting_Multitasking_Tracker
        (file,nt,nbequ,nbpts,ind,cnt,sup,r,mtype,perm,lif,mcc,mv,q,sols);
    else
      Silent_Multitasking_Tracker
        (nt,nbequ,nbpts,ind,cnt,sup,r,mtype,perm,lif,mcc,mv,q,sols);
      put(file,q'last,1); new_line(file);
      put(file,q);
    end if;
    new_line;
    put("The mixed volume : "); put(mv,1); new_line;
    new_line;
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
  end Random_Coefficient_System;

  procedure Random_Coefficient_System
              ( file : in file_type; nt : in integer32; stable : in boolean;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys ) is

  -- DESCRIPTION :
  --   Extracts the supports of the polynomial system p and computes
  --   its mixed volume and solves a random coefficient system,
  --   in double double precision.  The number of tasks equals nt.

    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Polynomial_Convertors;
    use DoblDobl_Complex_Solutions;

    stp : Standard_Complex_Laur_Systems.Laur_Sys(p'range)
        := DoblDobl_Complex_to_Standard_Laur_Sys(p);
    q : Laur_Sys(p'range);
    nbequ : constant integer32 := p'last;
    nbpts,r : integer32 := 0;
    stlb : double_float := 0.0;
    cnt,ind : Standard_Integer_Vectors.Vector(1..nbequ);
    sup,mtype,perm : Standard_Integer_Vectors.Link_to_Vector;
    lif : Link_to_Array_of_Lists;
    mcc : Mixed_Subdivision;
    mv : natural32;
    sols : Solution_List;
    ans : character;

  begin
    if stable then
      stlb := Floating_Lifting_Functions.Lifting_Bound(stp);
    else
      put("Test the pipelined driver ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        Pipelined_Polyhedral_Driver(file,nt,p);
        return;
      end if;
    end if;
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    Extract_Supports(nbequ,stp,nbpts,ind,cnt,sup);
    if ans = 'y' then
      Reporting_Multitasking_Tracker
        (file,nt,nbequ,nbpts,ind,cnt,sup,r,mtype,perm,lif,mcc,mv,q,sols);
    else
      Silent_Multitasking_Tracker
        (nt,nbequ,nbpts,ind,cnt,sup,r,mtype,perm,lif,mcc,mv,q,sols);
      put(file,q'last,1); new_line(file);
      put(file,q);
    end if;
    new_line;
    put("The mixed volume : "); put(mv,1); new_line;
    new_line;
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
  end Random_Coefficient_System;

  procedure Random_Coefficient_System
              ( file : in file_type; nt : in integer32; stable : in boolean;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys ) is

  -- DESCRIPTION :
  --   Extracts the supports of the polynomial system p and computes
  --   its mixed volume and solves a random coefficient system,
  --   in quad double precision.  The number of tasks equals nt.

    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Polynomial_Convertors;
    use QuadDobl_Complex_Solutions;

    stp : Standard_Complex_Laur_Systems.Laur_Sys(p'range)
        := QuadDobl_Complex_to_Standard_Laur_Sys(p);
    q : Laur_Sys(p'range);
    nbequ : constant integer32 := p'last;
    nbpts,r : integer32 := 0;
    stlb : double_float := 0.0;
    cnt,ind : Standard_Integer_Vectors.Vector(1..nbequ);
    sup,mtype,perm : Standard_Integer_Vectors.Link_to_Vector;
    lif : Link_to_Array_of_Lists;
    mcc : Mixed_Subdivision;
    mv : natural32;
    sols : Solution_List;
    ans : character;

  begin
    if stable then
      stlb := Floating_Lifting_Functions.Lifting_Bound(stp);
    else
      put("Test the pipelined driver ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        Pipelined_Polyhedral_Driver(file,nt,p);
        return;
      end if;
    end if;
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    Extract_Supports(nbequ,stp,nbpts,ind,cnt,sup);
    if ans = 'y' then
      Reporting_Multitasking_Tracker
        (file,nt,nbequ,nbpts,ind,cnt,sup,r,mtype,perm,lif,mcc,mv,q,sols);
    else
      Silent_Multitasking_Tracker
        (nt,nbequ,nbpts,ind,cnt,sup,r,mtype,perm,lif,mcc,mv,q,sols);
      put(file,q'last,1); new_line(file);
      put(file,q);
    end if;
    new_line;
    put("The mixed volume : "); put(mv,1); new_line;
    new_line;
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
  end Random_Coefficient_System;

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system, computes its supports,
  --   and then calls the routine to compute the mixed volume,
  --   or to solve a random coefficient system with polyhedral homotopies
  --   in double double precision.

    lp : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    file : file_type;
    nt : integer32 := 0;
    ans : character;
    stable : boolean;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(lp);
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    new_line;
    put("Give the number of tasks : "); get(nt);
    new_line;
    put("Do you want stable mixed volumes ? (y/n) ");
    Ask_Yes_or_No(ans);
    stable := (ans = 'y');
    new_line;
    put("Do you want a random coefficient system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' 
     then Random_Coefficient_System(file,nt,stable,lp.all);
     else Mixed_Volume_Calculation(file,nt,stable,lp.all);
    end if;
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system, computes its supports,
  --   and then calls the routine to compute the mixed volume,
  --   or to solve a random coefficient system with polyhedral homotopies
  --   in double double precision.

    lp : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    file : file_type;
    nt : integer32 := 0;
    ans : character;
    stable : boolean;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(lp);
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    new_line;
    put("Give the number of tasks : "); get(nt);
    new_line;
    put("Do you want stable mixed volumes ? (y/n) ");
    Ask_Yes_or_No(ans);
    stable := (ans = 'y');
    new_line;
    put("Do you want a random coefficient system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Random_Coefficient_System(file,nt,stable,lp.all);
    else
      declare
        use DoblDobl_Polynomial_Convertors;
        stp : Standard_Complex_Laur_Systems.Laur_Sys(lp'range)
            := DoblDobl_Complex_to_Standard_Laur_Sys(lp.all);
      begin
        Mixed_Volume_Calculation(file,nt,stable,stp);
      end;
    end if;
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system, computes its supports,
  --   and then calls the routine to compute the mixed volume,
  --   or to solve a random coefficient system with polyhedral homotopies
  --   in quad double precision.

    lp : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    file : file_type;
    nt : integer32 := 0;
    ans : character;
    stable : boolean;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(lp);
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    new_line;
    put("Give the number of tasks : "); get(nt);
    new_line;
    put("Do you want stable mixed volumes ? (y/n) ");
    Ask_Yes_or_No(ans);
    stable := (ans = 'y');
    new_line;
    put("Do you want a random coefficient system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Random_Coefficient_System(file,nt,stable,lp.all);
    else
      declare
        use QuadDobl_Polynomial_Convertors;
        stp : Standard_Complex_Laur_Systems.Laur_Sys(lp'range)
            := QuadDobl_Complex_to_Standard_Laur_Sys(lp.all);
      begin
        Mixed_Volume_Calculation(file,nt,stable,stp);
      end;
    end if;
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the level of the precision
  --   and then calls the proper driver.

    precision : constant character := Prompt_for_Precision;

  begin
    case precision is
      when '0' => Standard_Main;
      when '1' => DoblDobl_Main;
      when '2' => QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_mtmva;
