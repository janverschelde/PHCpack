with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Lists_of_Integer_Vectors;          use Lists_of_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Standard_Integer_VecVecs;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Floating_Mixed_Subdivisions;       use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;    use Floating_Mixed_Subdivisions_io;
with Cell_Stack;                        use Cell_Stack;
with Mixed_Volume;
with MixedVol_Algorithm;                use MixedVol_Algorithm;
with Mixed_Labels_Queue;
with Multitasking;

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

    mix : Standard_Integer_Vectors.Vector(1..r);

  begin
    for i in 1..r loop
      mix(i) := mtype(i-1);
    end loop;
    put(file,natural32(nbequ),mix,mcc);
  end Write_Mixed_Cells;

  procedure Mixed_Cell_Configuration
              ( nbequ,r,size,nb : in integer32;
                mtype,perm : in Standard_Integer_Vectors.Link_to_Vector;
                Vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                lft : in Standard_Floating_Vectors.Link_to_Vector;
                labels : in List; sub : out Mixed_Subdivision ) is

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

    last : Mixed_Subdivision;
    tmp : List := labels;
    lbl : Standard_Integer_Vectors.Link_to_Vector;

  begin
    put("Permutation of the supports : "); put(perm); new_line;
    put_line("Creating a regular mixed-cell configuration ...");
    new_line;
    put_line("See the output file for the mixed subdivision...");
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

  procedure Sequential_Mixed_Volume_Computation
              ( nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                sub : out Mixed_Subdivision ) is

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
  --            all points in the supports.

  -- ON RETURN :
  --   r        number of distinct supports;
  --   mtype    the type of mixture of the supports;
  --   perm     permutation of the supports;
  --   sub      a mixed cell configuration for a random lifting.

    stlb : constant double_float := 0.0;
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
    mv_lift(nbequ,nbpts,ind,cnt,support.all,stlb,r,idx,vtx,lft);
    size := Mixed_Volume.cell_size(r,mtype);
    Cs_Init(cells,size);
    Mixed_Volume.MixedVol_with_Callback
      (nbequ,r,size,mtype,idx,vtx,lft,nb,cells,mixvol,false,
       write_labels'access);
    put("The mixed volume is "); put(mixvol,1); put_line(".");
    put("There are "); put(nb,1); put_line(" mixed cells.");
    Mixed_Cell_Configuration(nbequ,r,size,nb,mtype,perm,Vtx,lft,labels,sub);
  end Sequential_Mixed_Volume_Computation;

  procedure Produce_Cells
              ( nbequ,nbpts,r : in integer32;
                mtype,idx : in Standard_Integer_Vectors.Link_to_Vector;
                vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                lft : in Standard_Floating_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   This code calls the mixedvol algorithm with a callback function.
  --   The callback function places the labels to the points in the
  --   mixed cells into a queue.

  -- ON ENTRY :
  --   nbequ    the number of equations in the input Laurent system;
  --   nbpts    the total number of points in the supports;
  --   r        number of different supports;
  --   mtype    type of mixture;
  --   idx      index to the vertex set;
  --   vtx      vertex points;
  --   lft      lifting values for the vertex points.

    cells : CellStack;
    cellcnt : natural32 := 0;
    CellSize : constant integer32 := Mixed_Volume.cell_size(r,mtype);
    nbcells : integer32 := 0;
    mixvol : natural32;

    procedure write_labels ( pts : Standard_Integer_Vectors.Link_to_Vector ) is
    begin
      cellcnt := cellcnt + 1;
      put_line("Appending cell " & Multitasking.to_string(cellcnt)
                                 & " to the queue");
      Mixed_Labels_Queue.Append(pts);
    end write_labels;

  begin
    put_line("starting the cell production ...");
    Cs_Init(cells,CellSize);
    Mixed_Volume.MixedVol_with_Callback
      (nbequ,r,CellSize,mtype,idx,vtx,lft,nbcells,cells,mixvol,false,
       write_labels'access);
    Mixed_Labels_Queue.Stop;
    put_line("The mixed volume is " & Multitasking.to_string(mixvol) & ".");
    put("There are " & Multitasking.to_string(nbcells) & " mixed cells.");
    Cs_Del(cells);
  end Produce_Cells;

  procedure Process_Cells 
              ( idtask,nbequ,nbpts,r : in integer32;
                mtype,perm : in Standard_Integer_Vectors.Link_to_Vector;
                vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                lft : in Standard_Floating_Vectors.Link_to_Vector;
                mcc : in out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   This code is executed by task with identification number
  --   equal to idtask.  The labels to the mixed cells are converted
  --   into mixed cells and stored in a mixed cell configuration.

  -- ON ENTRY :
  --   idtask   identification number of the task;
  --   nbequ    the number of equations in the input Laurent system;
  --   nbpts    the total number of points in the supports;
  --   r        the number of distinct supports;
  --   mtype    type of mixture;
  --   perm     permutation used to permute the supports;
  --   vtx      coordinates of the vertex points;
  --   lft      lifting values for the vertex points.

  -- ON RETURN :
  --   mcc      the mixed cells processed by task with id idtask.

    labels : Standard_Integer_Vectors.Link_to_Vector;
    mcc_last : Mixed_Subdivision;
    stop : boolean := false;
    cnt : natural32 := 0;

    use Standard_Integer_Vectors;

  begin
    while not stop loop
      labels := Mixed_Labels_Queue.Next;
      if labels = null then -- check if all labels are produced
        if Mixed_Labels_Queue.Stopped     -- production stopped
         then stop := true;   -- all labels have been processed
        end if;
      else
        cnt := cnt + 1;
        put_line("Task " & Multitasking.to_string(idtask)
                         & " processes cell "
                         & Multitasking.to_string(cnt));
        if r < nbequ then
          declare
            mic : constant Mixed_Cell
                := Labels_to_Mixed_Cell(nbequ,r,mtype,labels,Vtx,lft);
          begin
            Append(mcc,mcc_last,mic);
          end;
        else
          declare
            mic : constant Mixed_Cell
                := Labels_to_Mixed_Cell(nbequ,r,mtype,perm,labels,Vtx,lft);
          begin
            Append(mcc,mcc_last,mic);
          end;
        end if;
      end if;
    end loop;
    put_line("Task " & Multitasking.to_string(idtask) & " processed "
                     & Multitasking.to_string(Length_Of(mcc))
                     & " mixed cells.");
  end Process_Cells;

  procedure Multitasked_Mixed_Volume_Computation
              ( ntasks,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                sub : out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Extracts the supports of the polynomial system and computes
  --   its mixed volume.  In this multitasked version, the producing
  --   of the cells is interlacing with the processing of the cells.

  -- ON ENTRY :
  --   ntasks   the number of tasks;
  --   nbequ    the number of equations in the input Laurent system;
  --   nbpts    the total number of points in the supports;
  --   ind      ind(k) marks the beginning of the k-th support;
  --   cnt      cnt(k) counts the number of points in the k-th support;
  --   support  vector range 1..nbequ*nbpts with the coordinates of
  --            all points in the supports.

  -- ON RETURN :
  --   r        number of distinct supports;
  --   mtype    the type of mixture of the supports;
  --   perm     permutation of the supports;
  --   sub      a mixed cell configuration for a random lifting.

    stlb : constant double_float := 0.0; -- no stable mv for now...
    mcc : array(2..ntasks) of Mixed_Subdivision;
    idx,sdx,ndx : Standard_Integer_Vectors.Link_to_Vector;
    vtx,spt : Standard_Integer_VecVecs.Link_to_VecVec;
    lft : Standard_Floating_Vectors.Link_to_Vector;
    sub_last : Mixed_Subdivision;

    procedure do_job ( i,n : in integer32 ) is
    begin
      put_line("In do_job with task " & Multitasking.to_string(i));
      if i = 1
       then Produce_Cells(nbequ,nbpts,r,mtype,idx,vtx,lft);
       else Process_Cells(i,nbequ,nbpts,r,mtype,perm,vtx,lft,mcc(i));
      end if;
    end do_job;
    procedure do_jobs is new Multitasking.Reporting_Workers(do_job);

  begin
    Mixed_Labels_Queue.Start;
    mv_upto_pre4mv
      (nbequ,nbpts,ind,cnt,support.all,r,mtype,perm,idx,vtx,sdx,spt,ndx);
    mv_lift(nbequ,nbpts,ind,cnt,support.all,stlb,r,idx,vtx,lft);
    do_jobs(ntasks);
    for i in mcc'range loop
      Concat(sub,sub_last,mcc(i));
    end loop;
    Standard_Integer_Vectors.Clear(idx);
    Standard_Integer_Vectors.Clear(sdx);
    Standard_Integer_Vectors.Clear(ndx);
    Standard_Floating_Vectors.Clear(lft);
    Standard_Integer_VecVecs.Deep_Clear(vtx);
    Standard_Integer_VecVecs.Deep_Clear(spt);
  end Multitasked_Mixed_Volume_Computation;

  procedure Mixed_Volume_Calculation
              ( file : in file_type; p : in Laur_Sys ) is

  -- DESCRIPTION :
  --   Extracts the supports of the polynomial system and computes
  --   its mixed volume.  A mixed subdivision is written to file.

    nbequ : constant integer32 := p'last;
    nbpts,nt,r : integer32 := 0;
    cnt,ind : Standard_Integer_Vectors.Vector(1..nbequ);
    sup,mtype,perm : Standard_Integer_Vectors.Link_to_Vector;
    mcc : Mixed_Subdivision;

  begin
    put("Give the number of tasks : "); get(nt);
    Extract_Supports(nbequ,p,nbpts,ind,cnt,sup);
    if nt < 2 then
      Sequential_Mixed_Volume_Computation
        (nbequ,nbpts,ind,cnt,sup,r,mtype,perm,mcc);
    else
      Multitasked_Mixed_Volume_Computation
        (nt,nbequ,nbpts,ind,cnt,sup,r,mtype,perm,mcc);
    end if;
    Write_Mixed_Cells(file,nbequ,r,mtype,mcc);
  end Mixed_Volume_Calculation;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system, computes its supports,
  --   and then calls the routine to compute the mixed volume.

    lp : Link_to_Laur_Sys;
    file : file_type;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(lp);
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    new_line;
    Mixed_Volume_Calculation(file,lp.all);
  end Main;

begin
  Main;
end ts_mtmva;
