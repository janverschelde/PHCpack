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
with Mixed_Volume_Computation;
with Cell_Stack;                        use Cell_Stack;
with Mixed_Volume;
with MixedVol_Algorithm;                use MixedVol_Algorithm;
with Mixed_Labels_Queue;
with Semaphore;
with Multitasking;
with Pipelined_Labeled_Cells;           use Pipelined_Labeled_Cells;

procedure ts_mtmva is

-- DESCRIPTION :
--   Running the MixedVol algorithm with a callback function
--   from within a single task.

  function Mixture ( r : integer32;
                     mtype : in Standard_Integer_Vectors.Link_to_Vector )
                   return Standard_Integer_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the type of mixture as a vector where the entries
  --   start the count at one, instead of at zero as in mtype.

    res : Standard_Integer_Vectors.Vector(1..r);

  begin
    for i in 1..r loop
      res(i) := mtype(i-1);
    end loop;
    return res;
  end Mixture;

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
    ans : character;
    otp : boolean;
    sem : Semaphore.Lock;
    celcnt : natural32 := 0;
    sumvol : natural32 := 0;

    procedure Write_Mixed_Volume
                ( r : in integer32;
                  mtype : in Standard_Integer_Vectors.Link_to_Vector;
                  mic : in out Mixed_Cell ) is

    -- DESCRIPTION :
    --   Writes the mixed volume of the cell mic with type of mixture
    --   defined in r and mtype to screen, to test the callback procedure
    --   in the pipelined production of the mixed cells.

      mix : constant Standard_Integer_Vectors.Vector := Mixture(r,mtype);
      vol : natural32;

    begin
      Mixed_Volume_Computation.Mixed_Volume(r,mix,mic,vol);
      Semaphore.Request(sem);
      celcnt := celcnt + 1;
      sumvol := sumvol + vol;
      Semaphore.Release(sem);
      put_line("the mixed volume of cell " 
               & Multitasking.to_string(celcnt) & " is "
               & Multitasking.to_string(vol));
     -- put_line("the mixed volume of a cell : "
     --          & Multitasking.to_string(vol));
    end Write_Mixed_Volume;

  begin
    put("Give the number of tasks : "); get(nt);
    Extract_Supports(nbequ,p,nbpts,ind,cnt,sup);
    if nt < 2 then
      Sequential_Mixed_Volume_Computation
        (nbequ,nbpts,ind,cnt,sup,r,mtype,perm,mcc);
    else
      new_line;
      put("Monitor the progress of the computations ? (y/n) ");
      Ask_Yes_or_No(ans);
      otp := (ans = 'y');
      if otp then
        Pipelined_Mixed_Cells(nt,nbequ,nbpts,otp,ind,cnt,sup,r,mtype,perm,mcc,
          Write_Mixed_Volume'access);
      else
        Pipelined_Mixed_Cells(nt,nbequ,nbpts,otp,ind,cnt,sup,r,mtype,perm,mcc);
      end if;
      put("The sum of the volumes of all cells : ");
      put(sumvol,1); new_line;
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
