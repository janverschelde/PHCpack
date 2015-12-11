with Timing_Package;                     use Timing_Package;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
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
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;     use Floating_Mixed_Subdivisions_io;
with MixedVol_Algorithm;                 use MixedVol_Algorithm;
with Drivers_for_Static_Lifting;
with Pipelined_Labeled_Cells;
with Pipelined_Polyhedral_Trackers;      use Pipelined_Polyhedral_Trackers;

package body Pipelined_Polyhedral_Drivers is

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
    stlb : double_float := 0.0;
    r : integer32;
    mtype,perm : Standard_Integer_Vectors.Link_to_Vector;
    mcc : Mixed_Subdivision;

  begin
    tstart(timer);
    Extract_Supports(nbequ,p,nbpts,ind,cnt,sup);
    if contrep then
      Reporting_Multitasking_Tracker
        (file,nt,nbequ,nbpts,ind,cnt,sup,stlb,r,mtype,perm,mcc,mv,q,qsols);
    else
      Silent_Multitasking_Tracker
        (nt,nbequ,nbpts,ind,cnt,sup,stlb,r,mtype,perm,mcc,mv,q,qsols);
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
    stlb : double_float := 0.0;
    r : integer32;
    mtype,perm : Standard_Integer_Vectors.Link_to_Vector;
    mcc : Mixed_Subdivision;

  begin
    tstart(timer);
    Extract_Supports(nbequ,stp,nbpts,ind,cnt,sup);
    if contrep then
      Reporting_Multitasking_Tracker
        (file,nt,nbequ,nbpts,ind,cnt,sup,stlb,r,mtype,perm,mcc,mv,q,qsols);
    else
      Silent_Multitasking_Tracker
        (nt,nbequ,nbpts,ind,cnt,sup,stlb,r,mtype,perm,mcc,mv,q,qsols);
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
    stlb : double_float := 0.0;
    r : integer32;
    mtype,perm : Standard_Integer_Vectors.Link_to_Vector;
    mcc : Mixed_Subdivision;

  begin
    tstart(timer);
    Extract_Supports(nbequ,stp,nbpts,ind,cnt,sup);
    if contrep then
      Reporting_Multitasking_Tracker
        (file,nt,nbequ,nbpts,ind,cnt,sup,stlb,r,mtype,perm,mcc,mv,q,qsols);
    else
      Silent_Multitasking_Tracker
        (nt,nbequ,nbpts,ind,cnt,sup,stlb,r,mtype,perm,mcc,mv,q,qsols);
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
  end Pipelined_Polyhedral_Homotopies;

end Pipelined_Polyhedral_Drivers;
