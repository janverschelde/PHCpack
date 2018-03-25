with Timing_Package;                     use Timing_Package;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Laur_Poly_Convertors;
with Standard_Poly_Laur_Convertors;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with DoblDobl_Polynomial_Convertors;
with DoblDobl_Laur_Poly_Convertors;
with DoblDobl_Poly_Laur_Convertors;
with QuadDobl_Polynomial_Convertors;
with QuadDobl_Laur_Poly_Convertors;
with QuadDobl_Poly_Laur_Convertors;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Floating_Mixed_Subdivisions_io;     use Floating_Mixed_Subdivisions_io;
with MixedVol_Algorithm;                 use MixedVol_Algorithm;
with Drivers_for_Static_Lifting;
with Pipelined_Labeled_Cells;
with Pipelined_Polyhedral_Trackers;      use Pipelined_Polyhedral_Trackers;

package body Pipelined_Polyhedral_Drivers is

-- SILENT VERSIONS, ON ORDINARY POLYNOMIAL SYSTEMS :

  procedure Pipelined_Polyhedral_Homotopies
              ( nt : in integer32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                mv : out natural32;
                q : out Standard_Complex_Poly_Systems.Poly_Sys;
                qsols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Poly_Laur_Convertors;
    use Standard_Laur_Poly_Convertors;

    lp,lq : Standard_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    lp := Polynomial_to_Laurent_System(p);
    Pipelined_Polyhedral_Homotopies(nt,lp,mv,lq,qsols);
    q := Laurent_to_Polynomial_System(lq);
    Standard_Complex_Laur_Systems.Clear(lp);
    Standard_Complex_Laur_Systems.Clear(lq);
  end Pipelined_Polyhedral_Homotopies;

  procedure Pipelined_Polyhedral_Homotopies
              ( nt : in integer32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                mv : out natural32;
                q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                qsols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Poly_Laur_Convertors;
    use DoblDobl_Laur_Poly_Convertors;

    lp,lq : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    lp := Polynomial_to_Laurent_System(p);
    Pipelined_Polyhedral_Homotopies(nt,lp,mv,lq,qsols);
    q := Laurent_to_Polynomial_System(lq);
    DoblDobl_Complex_Laur_Systems.Clear(lp);
    DoblDobl_Complex_Laur_Systems.Clear(lq);
  end Pipelined_Polyhedral_Homotopies;

  procedure Pipelined_Polyhedral_Homotopies
              ( nt : in integer32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                mv : out natural32;
                q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                qsols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Poly_Laur_Convertors;
    use QuadDobl_Laur_Poly_Convertors;

    lp,lq : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    lp := Polynomial_to_Laurent_System(p);
    Pipelined_Polyhedral_Homotopies(nt,lp,mv,lq,qsols);
    q := Laurent_to_Polynomial_System(lq);
    QuadDobl_Complex_Laur_Systems.Clear(lp);
    QuadDobl_Complex_Laur_Systems.Clear(lq);
  end Pipelined_Polyhedral_Homotopies;

-- VERSIONS WITH OUTPUT TO FILES, ON ORDINARY POLYNOMIAL SYSTEMS :

  procedure Pipelined_Polyhedral_Homotopies
              ( file,cfile,qfile : in file_type; nt : in integer32;
                misufile,contrep : in boolean;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                mv : out natural32;
                q : out Standard_Complex_Poly_Systems.Poly_Sys;
                qsols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Poly_Laur_Convertors;
    use Standard_Laur_Poly_Convertors;

    lp,lq : Standard_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    lp := Polynomial_to_Laurent_System(p);
    Pipelined_Polyhedral_Homotopies
      (file,cfile,qfile,nt,misufile,contrep,lp,mv,lq,qsols);
    q := Laurent_to_Polynomial_System(lq);
    Standard_Complex_Laur_Systems.Clear(lp);
    Standard_Complex_Laur_Systems.Clear(lq);
  end Pipelined_Polyhedral_Homotopies;

  procedure Pipelined_Polyhedral_Homotopies
              ( file,cfile,qfile : in file_type; nt : in integer32;
                misufile,contrep : in boolean;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                mv : out natural32;
                q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                qsols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Poly_Laur_Convertors;
    use DoblDobl_Laur_Poly_Convertors;

    lp,lq : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    lp := Polynomial_to_Laurent_System(p);
    Pipelined_Polyhedral_Homotopies
      (file,cfile,qfile,nt,misufile,contrep,lp,mv,lq,qsols);
    q := Laurent_to_Polynomial_System(lq);
    DoblDobl_Complex_Laur_Systems.Clear(lp);
    DoblDobl_Complex_Laur_Systems.Clear(lq);
  end Pipelined_Polyhedral_Homotopies;

  procedure Pipelined_Polyhedral_Homotopies
              ( file,cfile,qfile : in file_type; nt : in integer32;
                misufile,contrep : in boolean;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                mv : out natural32;
                q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                qsols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Poly_Laur_Convertors;
    use QuadDobl_Laur_Poly_Convertors;

    lp,lq : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    lp := Polynomial_to_Laurent_System(p);
    Pipelined_Polyhedral_Homotopies
      (file,cfile,qfile,nt,misufile,contrep,lp,mv,lq,qsols);
    q := Laurent_to_Polynomial_System(lq);
    QuadDobl_Complex_Laur_Systems.Clear(lp);
    QuadDobl_Complex_Laur_Systems.Clear(lq);
  end Pipelined_Polyhedral_Homotopies;

-- SILENT VERSIONS, ON LAURENT SYSTEMS :

  procedure Pipelined_Polyhedral_Homotopies
              ( nt : in integer32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                mv : out natural32;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                qsols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;
    use Drivers_for_Static_Lifting;

    nbequ : constant integer32 := p'last;
    nbpts : integer32 := 0;
    cnt,ind : Standard_Integer_Vectors.Vector(1..nbequ);
    sup : Standard_Integer_Vectors.Link_to_Vector;
    r : integer32;
    mtype,perm : Standard_Integer_Vectors.Link_to_Vector;
    lif : Link_to_Array_of_Lists;
    mcc : Mixed_Subdivision;

  begin
    Extract_Supports(nbequ,p,nbpts,ind,cnt,sup);
    Silent_Multitasking_Tracker
      (nt,nbequ,nbpts,ind,cnt,sup,r,mtype,perm,lif,mcc,mv,q,qsols);
    Standard_Integer_Vectors.Clear(sup);
    Clear(mcc);
  end Pipelined_Polyhedral_Homotopies;

  procedure Pipelined_Polyhedral_Homotopies
              ( nt : in integer32;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mv : out natural32;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                qsols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Polynomial_Convertors;
    use DoblDobl_Complex_Solutions;
    use Drivers_for_Static_Lifting;

    stp : Standard_Complex_Laur_Systems.Laur_Sys(p'range)
        := DoblDobl_Complex_to_Standard_Laur_Sys(p);
    nbequ : constant integer32 := p'last;
    nbpts : integer32 := 0;
    cnt,ind : Standard_Integer_Vectors.Vector(1..nbequ);
    sup : Standard_Integer_Vectors.Link_to_Vector;
    r : integer32;
    mtype,perm : Standard_Integer_Vectors.Link_to_Vector;
    lif : Link_to_Array_of_Lists;
    mcc : Mixed_Subdivision;

  begin
    Extract_Supports(nbequ,stp,nbpts,ind,cnt,sup);
    Silent_Multitasking_Tracker
      (nt,nbequ,nbpts,ind,cnt,sup,r,mtype,perm,lif,mcc,mv,q,qsols);
    Standard_Integer_Vectors.Clear(sup);
    Clear(mcc);
  end Pipelined_Polyhedral_Homotopies;

  procedure Pipelined_Polyhedral_Homotopies
              ( nt : in integer32;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mv : out natural32;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                qsols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Polynomial_Convertors;
    use QuadDobl_Complex_Solutions;
    use Drivers_for_Static_Lifting;

    stp : Standard_Complex_Laur_Systems.Laur_Sys(p'range)
        := QuadDobl_Complex_to_Standard_Laur_Sys(p);
    nbequ : constant integer32 := p'last;
    nbpts : integer32 := 0;
    cnt,ind : Standard_Integer_Vectors.Vector(1..nbequ);
    sup : Standard_Integer_Vectors.Link_to_Vector;
    r : integer32;
    mtype,perm : Standard_Integer_Vectors.Link_to_Vector;
    lif : Link_to_Array_of_Lists;
    mcc : Mixed_Subdivision;

  begin
    Extract_Supports(nbequ,stp,nbpts,ind,cnt,sup);
    Silent_Multitasking_Tracker
      (nt,nbequ,nbpts,ind,cnt,sup,r,mtype,perm,lif,mcc,mv,q,qsols);
    Standard_Integer_Vectors.Clear(sup);
    Clear(mcc);
  end Pipelined_Polyhedral_Homotopies;

-- VERSIONS WITH OUTPUT TO FILES, ON LAURENT SYSTEMS :

  procedure Pipelined_Polyhedral_Homotopies
              ( file,cfile,qfile : in file_type; nt : in integer32;
                misufile,contrep : in boolean;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                mv : out natural32;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                qsols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;
    use Drivers_for_Static_Lifting;

    timer : Timing_Widget;
    nbequ : constant integer32 := p'last;
    nbpts : integer32 := 0;
    cnt,ind : Standard_Integer_Vectors.Vector(1..nbequ);
    sup : Standard_Integer_Vectors.Link_to_Vector;
    r : integer32;
    mtype,perm : Standard_Integer_Vectors.Link_to_Vector;
    lif : Link_to_Array_of_Lists;
    mcc : Mixed_Subdivision;

  begin
    tstart(timer);
    Extract_Supports(nbequ,p,nbpts,ind,cnt,sup);
    if contrep then
      Reporting_Multitasking_Tracker
        (file,nt,nbequ,nbpts,ind,cnt,sup,r,mtype,perm,lif,mcc,mv,q,qsols);
    else
      Silent_Multitasking_Tracker
        (nt,nbequ,nbpts,ind,cnt,sup,r,mtype,perm,lif,mcc,mv,q,qsols);
    end if;
    tstop(timer);
    put(file,q'last,1); new_line(file);
    put(file,q);
    new_line(file);
    print_times(file,timer,"pipelined polyhedral path trackers");
    declare
      mix : constant Standard_Integer_Vectors.Vector
          := Pipelined_Labeled_Cells.Mixture(r,mtype);
    begin
      Floating_Volume_Computation(file,nbequ,mix,mcc,mv);
      if misufile
       then put(cfile,natural32(nbequ),mix,mcc);
      end if;
    end;
    put(qfile,q'last,1); new_line(qfile);
    put(qfile,q);
    if Length_Of(qsols) > 0 then
      new_line(qfile);
      put_line(qfile,"THE SOLUTIONS :");
      put(qfile,Length_Of(qsols),natural32(Head_Of(qsols).n),qsols);
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,Length_Of(qsols),natural32(Head_Of(qsols).n),qsols);
    end if;
    Standard_Integer_Vectors.Clear(sup);
    Clear(mcc);
  end Pipelined_Polyhedral_Homotopies;

  procedure Pipelined_Polyhedral_Homotopies
              ( file,cfile,qfile : in file_type; nt : in integer32;
                misufile,contrep : in boolean;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mv : out natural32;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                qsols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Polynomial_Convertors;
    use DoblDobl_Complex_Solutions;
    use Drivers_for_Static_Lifting;

    stp : Standard_Complex_Laur_Systems.Laur_Sys(p'range)
        := DoblDobl_Complex_to_Standard_Laur_Sys(p);
    timer : Timing_Widget;
    nbequ : constant integer32 := p'last;
    nbpts : integer32 := 0;
    cnt,ind : Standard_Integer_Vectors.Vector(1..nbequ);
    sup : Standard_Integer_Vectors.Link_to_Vector;
    r : integer32;
    mtype,perm : Standard_Integer_Vectors.Link_to_Vector;
    lif : Link_to_Array_of_Lists;
    mcc : Mixed_Subdivision;

  begin
    tstart(timer);
    Extract_Supports(nbequ,stp,nbpts,ind,cnt,sup);
    if contrep then
      Reporting_Multitasking_Tracker
        (file,nt,nbequ,nbpts,ind,cnt,sup,r,mtype,perm,lif,mcc,mv,q,qsols);
    else
      Silent_Multitasking_Tracker
        (nt,nbequ,nbpts,ind,cnt,sup,r,mtype,perm,lif,mcc,mv,q,qsols);
    end if;
    tstop(timer);
    put(file,q'last,1); new_line(file);
    put(file,q);
    new_line(file);
    print_times(file,timer,"pipelined polyhedral path trackers");
    declare
      mix : constant Standard_Integer_Vectors.Vector
          := Pipelined_Labeled_Cells.Mixture(r,mtype);
    begin
      Floating_Volume_Computation(file,nbequ,mix,mcc,mv);
      if misufile
       then put(cfile,natural32(nbequ),mix,mcc);
      end if;
    end;
    put(qfile,q'last,1); new_line(qfile);
    put(qfile,q);
    if Length_Of(qsols) > 0 then
      new_line(qfile);
      put_line(qfile,"THE SOLUTIONS :");
      put(qfile,Length_Of(qsols),natural32(Head_Of(qsols).n),qsols);
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,Length_Of(qsols),natural32(Head_Of(qsols).n),qsols);
    end if;
    Standard_Integer_Vectors.Clear(sup);
    Clear(mcc);
  end Pipelined_Polyhedral_Homotopies;

  procedure Pipelined_Polyhedral_Homotopies
              ( file,cfile,qfile : in file_type; nt : in integer32;
                misufile,contrep : in boolean;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mv : out natural32;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                qsols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Polynomial_Convertors;
    use QuadDobl_Complex_Solutions;
    use Drivers_for_Static_Lifting;

    stp : Standard_Complex_Laur_Systems.Laur_Sys(p'range)
        := QuadDobl_Complex_to_Standard_Laur_Sys(p);
    timer : Timing_Widget;
    nbequ : constant integer32 := p'last;
    nbpts : integer32 := 0;
    cnt,ind : Standard_Integer_Vectors.Vector(1..nbequ);
    sup : Standard_Integer_Vectors.Link_to_Vector;
    r : integer32;
    mtype,perm : Standard_Integer_Vectors.Link_to_Vector;
    lif : Link_to_Array_of_Lists;
    mcc : Mixed_Subdivision;

  begin
    tstart(timer);
    Extract_Supports(nbequ,stp,nbpts,ind,cnt,sup);
    if contrep then
      Reporting_Multitasking_Tracker
        (file,nt,nbequ,nbpts,ind,cnt,sup,r,mtype,perm,lif,mcc,mv,q,qsols);
    else
      Silent_Multitasking_Tracker
        (nt,nbequ,nbpts,ind,cnt,sup,r,mtype,perm,lif,mcc,mv,q,qsols);
    end if;
    tstop(timer);
    put(file,q'last,1); new_line(file);
    put(file,q);
    new_line(file);
    print_times(file,timer,"pipelined polyhedral path trackers");
    declare
      mix : constant Standard_Integer_Vectors.Vector
          := Pipelined_Labeled_Cells.Mixture(r,mtype);
    begin
      Floating_Volume_Computation(file,nbequ,mix,mcc,mv);
      if misufile
       then put(cfile,natural32(nbequ),mix,mcc);
      end if;
    end;
    put(qfile,q'last,1); new_line(qfile);
    put(qfile,q);
    if Length_Of(qsols) > 0 then
      new_line(qfile);
      put_line(qfile,"THE SOLUTIONS :");
      put(qfile,Length_Of(qsols),natural32(Head_Of(qsols).n),qsols);
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,Length_Of(qsols),natural32(Head_Of(qsols).n),qsols);
    end if;
    Standard_Integer_Vectors.Clear(sup);
    Clear(mcc);
  end Pipelined_Polyhedral_Homotopies;

-- WITH LIFTING BOUND FOR THE ARTIFICIAL ORIGIN :

  procedure Pipelined_Polyhedral_Homotopies
              ( nt : in integer32;
                stable : in boolean; stlb : in double_float;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                lq : out Standard_Complex_Laur_Systems.Laur_Sys;
                q : out Standard_Complex_Poly_Systems.Poly_Sys;
                qsols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Poly_Laur_Convertors;
    use Standard_Laur_Poly_Convertors;

    lp : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    nbequ : constant integer32 := p'last;
    nbpts : integer32 := 0;
    cnt,ind : Standard_Integer_Vectors.Vector(1..nbequ);
    sup : Standard_Integer_Vectors.Link_to_Vector;

  begin
    lp := Polynomial_to_Laurent_System(p);
    Extract_Supports(nbequ,lp,nbpts,ind,cnt,sup);
    Silent_Multitasking_Tracker
      (nt,nbequ,nbpts,ind,cnt,sup,stable,stlb,
       r,mtype,perm,lif,mcc,mv,lq,qsols);
    q := Laurent_to_Polynomial_System(lq);
    Standard_Complex_Laur_Systems.Clear(lp);
  end Pipelined_Polyhedral_Homotopies;

  procedure Pipelined_Polyhedral_Homotopies
              ( nt : in integer32;
                stable : in boolean; stlb : in double_float;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                lq : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                qsols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Poly_Laur_Convertors;
    use DoblDobl_Laur_Poly_Convertors;
    use DoblDobl_Polynomial_Convertors;

    lp : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range)
       := Polynomial_to_Laurent_System(p);
    stp : Standard_Complex_Laur_Systems.Laur_Sys(lp'range)
        := DoblDobl_Complex_to_Standard_Laur_Sys(lp);
    nbequ : constant integer32 := p'last;
    nbpts : integer32 := 0;
    cnt,ind : Standard_Integer_Vectors.Vector(1..nbequ);
    sup : Standard_Integer_Vectors.Link_to_Vector;

  begin
    Extract_Supports(nbequ,stp,nbpts,ind,cnt,sup);
    Silent_Multitasking_Tracker
      (nt,nbequ,nbpts,ind,cnt,sup,stable,stlb,
       r,mtype,perm,lif,mcc,mv,lq,qsols);
    q := Laurent_to_Polynomial_System(lq);
    DoblDobl_Complex_Laur_Systems.Clear(lp);
    Standard_Complex_Laur_Systems.Clear(stp);
  end Pipelined_Polyhedral_Homotopies;

  procedure Pipelined_Polyhedral_Homotopies
              ( nt : in integer32;
                stable : in boolean; stlb : in double_float;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                lq : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                qsols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Poly_Laur_Convertors;
    use QuadDobl_Laur_Poly_Convertors;
    use QuadDobl_Polynomial_Convertors;

    lp : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range)
       := Polynomial_to_Laurent_System(p);
    stp : Standard_Complex_Laur_Systems.Laur_Sys(lp'range)
        := QuadDobl_Complex_to_Standard_Laur_Sys(lp);
    nbequ : constant integer32 := p'last;
    nbpts : integer32 := 0;
    cnt,ind : Standard_Integer_Vectors.Vector(1..nbequ);
    sup : Standard_Integer_Vectors.Link_to_Vector;

  begin
    Extract_Supports(nbequ,stp,nbpts,ind,cnt,sup);
    Silent_Multitasking_Tracker
      (nt,nbequ,nbpts,ind,cnt,sup,stable,stlb,
       r,mtype,perm,lif,mcc,mv,lq,qsols);
    q := Laurent_to_Polynomial_System(lq);
    QuadDobl_Complex_Laur_Systems.Clear(lp);
    Standard_Complex_Laur_Systems.Clear(stp);
  end Pipelined_Polyhedral_Homotopies;

end Pipelined_Polyhedral_Drivers;
