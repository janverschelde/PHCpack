with Ada.Calendar;
with Timing_Package,Time_Stamps;         use Timing_Package,Time_Stamps;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Multprec_Natural_Numbers_io;        use Multprec_Natural_Numbers_io;
with Characters_and_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with DoblDobl_Polynomial_Convertors;     use DoblDobl_Polynomial_Convertors;
with QuadDobl_Polynomial_Convertors;     use QuadDobl_Polynomial_Convertors;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Total_Degree_Start_Systems;         use Total_Degree_Start_Systems;
with m_Homogeneous_Bezout_Numbers;       use m_Homogeneous_Bezout_Numbers;
with m_Homogeneous_Start_Systems;        use m_Homogeneous_Start_Systems;
with Set_Structure;
with Standard_Linear_Product_System;
with Standard_Complex_Prod_Systems;      use Standard_Complex_Prod_Systems;
with Standard_Complex_Prod_Systems_io;   use Standard_Complex_Prod_Systems_io;
with Standard_Complex_Prod_Planes;
with Random_Product_Start_Systems;       use Random_Product_Start_Systems;
with Floating_Lifting_Functions;
with Floating_Mixed_Subdivisions_io;
with Induced_Permutations;
with Apply_Induced_Permutations;         use Apply_Induced_Permutations;
with Mixed_Volume_Computation;
with Black_Mixed_Volume_Computations;    use Black_Mixed_Volume_Computations;
with Stable_Polyhedral_Continuation;     use Stable_Polyhedral_Continuation;
with Black_Polyhedral_Continuations;     use Black_Polyhedral_Continuations;
with Root_Counters_Output;
with Pipelined_Labeled_Cells;
with Pipelined_Polyhedral_Drivers;

package body Black_Box_Root_Counters is

  chicken_mv : constant natural32 := 40;  -- arbitrary bound ...

  procedure Compute_Total_Degree
               ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 tdeg : out natural64; mptdeg : out Natural_Number ) is

  begin
    tdeg := Total_Degree(p);
  exception
    when others =>
      mptdeg := Total_Degree(p); tdeg := 0;
      put("Huge total degree : "); put(mptdeg); new_line;
  end Compute_Total_Degree;

  procedure Compute_Total_Degree
               ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 tdeg : out natural64; mptdeg : out Natural_Number ) is

    q : Standard_Complex_Poly_Systems.Poly_Sys(p'range)
      := QuadDobl_Complex_to_Standard_Poly_Sys(p);

  begin
    Compute_Total_Degree(q,tdeg,mptdeg);
  end Compute_Total_Degree;

  procedure Compute_Total_Degree
               ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 tdeg : out natural64; mptdeg : out Natural_Number ) is

    q : Standard_Complex_Poly_Systems.Poly_Sys(p'range)
      := DoblDobl_Complex_to_Standard_Poly_Sys(p);

  begin
    Compute_Total_Degree(q,tdeg,mptdeg);
  end Compute_Total_Degree;

  function Set_Structure_Bound
             ( p : Standard_Complex_Poly_Systems.Poly_Sys ) return natural64 is
  begin
    Random_Product_Start_Systems.Build_Set_Structure(p);
   -- put_line("The set structure :"); Set_Structure_io.put;
    Standard_Linear_Product_System.Init(natural32(p'last));
    Random_Product_Start_Systems.Build_Random_Product_System(natural32(p'last));
    return Standard_Linear_Product_System.Count_All_Solutions(1.0E-8);
   -- return Permanent(Degree_Sets_Tables.Create);
  end Set_Structure_Bound;

  function Set_Structure_Bound
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys ) return natural64 is

    q : Standard_Complex_Poly_Systems.Poly_Sys(p'range)
      := DoblDobl_Complex_to_Standard_Poly_Sys(p);
    res : constant natural64 := Set_Structure_Bound(q);

  begin
    Standard_Complex_Poly_Systems.Clear(q);
    return res;
  end Set_Structure_Bound;

  function Set_Structure_Bound
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys ) return natural64 is

    q : Standard_Complex_Poly_Systems.Poly_Sys(p'range)
      := QuadDobl_Complex_to_Standard_Poly_Sys(p);
    res : constant natural64 := Set_Structure_Bound(q);

  begin
    Standard_Complex_Poly_Systems.Clear(q);
    return res;
  end Set_Structure_Bound;

  procedure Count_Roots 
               ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean;
                 tode : out natural64; mptode : out Natural_Number;
                 mhbz,setb : out natural64;
                 mivo,stmv : out natural32;
                 zz : out Partition; nz : out natural32;
                 stlb : out double_float;
                 lifsup : out Link_to_Array_of_Lists;
                 mix,perm,iprm : out Link_to_Vector;
                 orgmcc,stbmcc : out Mixed_Subdivision;
                 rocotime : out duration; verbose : in integer32 := 0 ) is

    timer : timing_widget;
    n : constant natural32 := natural32(p'length);
    d,bz,bs : natural64 := 0;
    m,mv,smv,tmv : natural32 := 0;
    z : Partition(1..n);
    no_mv : constant boolean := (n > chicken_mv);
    mixsub : Mixed_Subdivision;
    orgcnt,stbcnt : natural32;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Count_Roots 1,");
      put_line("for polynomial systems in double precision ...");
    end if;
    tstart(timer);
    Compute_Total_Degree(p,d,mptode);
    if d = 0 and Equal(mptode,0)       -- patch for GNAT optimizers ...
     then mptode := Total_Degree(p);
    end if;
    if verbose > 0
     then put("-> total degree : "); put(d,1); new_line;
    end if;
    declare
    begin
      PB(p,bz,m,z);
      if verbose > 0 then
        put("-> "); put(m,1); put("-homogeneous Bezout number : ");
        put(bz,1); new_line;
      end if;
      bs := Set_Structure_Bound(p);
      if verbose > 0 then
        put("-> set structure Bezout bound : "); put(bs,1); new_line;
      end if;
    exception
      when others => m := 1; bz := 0; bs := 0;
    end;
    if m = 1
     then bz := d;
    end if;
    if not deg and not no_mv then
      declare -- problems with systems with one monomial equation
      begin
        Black_Box_Mixed_Volume_Computation
          (p,mix,perm,iprm,stlb,lifsup,mixsub,orgmcc,stbmcc,
           mv,smv,tmv,orgcnt,stbcnt,verbose-1);
      exception
        when others => mv := 0; smv := 0; tmv := 0;
      end;
    end if;
    tstop(timer);
    if deg or no_mv
     then mivo := 0; stmv := 0;
     else mivo := mv; stmv := smv;
    end if;
    tode := d; mhbz := bz; setb := bs;
    zz := z; nz := m;
    rocotime := Elapsed_User_Time(timer);
  end Count_Roots;

  procedure Count_Roots 
               ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean;
                 tode : out natural64; mptode : out Natural_Number;
                 mhbz,setb : out natural64;
                 mivo,stmv : out natural32;
                 zz : out Partition; nz : out natural32;
                 stlb : out double_float;
                 lifsup : out Link_to_Array_of_Lists;
                 mix,perm,iprm : out Link_to_Vector;
                 orgmcc,stbmcc : out Mixed_Subdivision;
                 rocotime : out duration; verbose : in integer32 := 0 ) is

    q : Standard_Complex_Poly_Systems.Poly_Sys(p'range)
      := DoblDobl_Complex_to_Standard_Poly_Sys(p);

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Count_Roots 2,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    Count_Roots(q,deg,tode,mptode,mhbz,setb,mivo,stmv,zz,nz,stlb,lifsup,
                mix,perm,iprm,orgmcc,stbmcc,rocotime,verbose-1);
    if iprm /= null then
     -- put("The permutation : "); put(iprm); new_line;
     -- put_line("The system before the permutation : "); put_line(p);
      Induced_Permutations.Permute(iprm.all,p);
     -- put_line("The system after the permutation : "); put_line(p);
    end if;
    Standard_Complex_Poly_Systems.Clear(q);
  end Count_Roots;

  procedure Count_Roots 
               ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean;
                 tode : out natural64; mptode : out Natural_Number;
                 mhbz,setb : out natural64;
                 mivo,stmv : out natural32;
                 zz : out Partition; nz : out natural32;
                 stlb : out double_float;
                 lifsup : out Link_to_Array_of_Lists;
                 mix,perm,iprm : out Link_to_Vector;
                 orgmcc,stbmcc : out Mixed_Subdivision;
                 rocotime : out duration; verbose : in integer32 := 0 ) is

    q : Standard_Complex_Poly_Systems.Poly_Sys(p'range)
      := QuadDobl_Complex_to_Standard_Poly_Sys(p);

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Count_Roots 3,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    Count_Roots(q,deg,tode,mptode,mhbz,setb,mivo,stmv,zz,nz,stlb,lifsup,
                mix,perm,iprm,orgmcc,stbmcc,rocotime,verbose-1);
    if iprm /= null
     then Induced_Permutations.Permute(iprm.all,p);
    end if;
    Standard_Complex_Poly_Systems.Clear(q);
  end Count_Roots;

  procedure Count_Roots 
               ( file : in file_type;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean;
                 tode : out natural64; mptode : out Natural_Number;
                 mhbz,setb : out natural64;
                 mivo,stmv : out natural32;
                 zz : out Partition; nz : out natural32;
                 stlb : out double_float;
                 lifsup : out Link_to_Array_of_Lists;
                 mix,perm,iprm : out Link_to_Vector;
                 orgmcc,stbmcc : out Mixed_Subdivision;
                 rocotime : out duration; verbose : in integer32 := 0 ) is

    timer : timing_widget;
    n : constant natural32 := natural32(p'length);
    d,bz,bs : natural64 := 0;
    m,mv,smv,tmv,orgcnt,stbcnt : natural32 := 0;
    z : Partition(1..n);
    no_mv : constant boolean := deg or (n > chicken_mv);
    mixsub : Mixed_Subdivision;

    use Root_Counters_Output;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Count_Roots 4,");
      put_line("for polynomial systems in double precision ...");
    end if;
    tstart(timer);
    declare
    begin
      d := Total_Degree(p);
    exception
      when others => 
        mptode := Total_Degree(p); d := 0;
        put("huge total degree : "); put(mptode); new_line;
    end;
    if d = 0 and Equal(mptode,0)      -- patch for GNAT optimizers ...
     then mptode := Total_Degree(p);
    end if;
    if verbose > 0
     then put("-> total degree : "); put(d,1); new_line;
    end if;
    declare
    begin
      PB(p,bz,m,z); 
      if verbose > 0 then
        put("-> "); put(m,1); put("-homogeneous Bezout number : ");
        put(bz,1); new_line;
      end if;
      bs := Set_Structure_Bound(p);
      if verbose > 0 then
        put("-> set structure Bezout bound : "); put(bs,1); new_line;
      end if;
    exception
      when others => bz := 0; m := 1; bs := 0;
    end;
    if m = 1
     then bz := d;
    end if;
    if deg or not no_mv then
      declare -- problems with systems with one monomial equation
      begin
        Black_Box_Mixed_Volume_Computation
          (p,mix,perm,iprm,stlb,lifsup,mixsub,orgmcc,stbmcc,
           mv,smv,tmv,orgcnt,stbcnt,verbose-1);
      exception
        when others => mv := 0; smv := 0; tmv := 0;
      end;
    end if;
    tstop(timer);
    Write_Root_Counts(file,no_mv,d,mptode,m,bz,bs,mv,smv,z);
    new_line(file);
    print_times(file,timer,"Root Counting");
    flush(file);
    tode := d; mhbz := bz; setb := bs;
    zz := z; nz := m;
    if no_mv
     then mivo := 0; stmv := 0;
     else mivo := mv; stmv := smv;
    end if;
    rocotime := Elapsed_User_Time(timer);
  end Count_Roots;

  procedure Count_Roots 
               ( file : in file_type;
                 p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean;
                 tode : out natural64; mptode : out Natural_Number;
                 mhbz,setb : out natural64;
                 mivo,stmv : out natural32;
                 zz : out Partition; nz : out natural32;
                 stlb : out double_float;
                 lifsup : out Link_to_Array_of_Lists;
                 mix,perm,iprm : out Link_to_Vector;
                 orgmcc,stbmcc : out Mixed_Subdivision;
                 rocotime : out duration; verbose : in integer32 := 0 ) is

    sp : Standard_Complex_Poly_Systems.Poly_Sys(p'range)
       := DoblDobl_Complex_to_Standard_Poly_Sys(p);

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Count_Roots 5,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    Count_Roots(file,sp,deg,tode,mptode,mhbz,setb,mivo,stmv,zz,nz,stlb,
                lifsup,mix,perm,iprm,orgmcc,stbmcc,rocotime,verbose-1);
    if iprm /= null then
     -- put("The permutation : "); put(iprm); new_line;
     -- put_line("The system before the permutation : "); put_line(p);
      Induced_Permutations.Permute(iprm.all,p);
     -- put_line("The system after the permutation : "); put_line(p);
    end if;
    Standard_Complex_Poly_Systems.Clear(sp);
  end Count_Roots;

  procedure Count_Roots 
               ( file : in file_type;
                 p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean;
                 tode : out natural64; mptode : out Natural_Number;
                 mhbz,setb : out natural64;
                 mivo,stmv : out natural32;
                 zz : out Partition; nz : out natural32;
                 stlb : out double_float;
                 lifsup : out Link_to_Array_of_Lists;
                 mix,perm,iprm : out Link_to_Vector;
                 orgmcc,stbmcc : out Mixed_Subdivision;
                 rocotime : out duration; verbose : in integer32 := 0 ) is

    sp : Standard_Complex_Poly_Systems.Poly_Sys(p'range)
       := QuadDobl_Complex_to_Standard_Poly_Sys(p);

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Count_Roots 6,");
      put_line("for polynomial systems in double precision ...");
    end if;
    Count_Roots(file,sp,deg,tode,mptode,mhbz,setb,mivo,stmv,zz,nz,stlb,
                lifsup,mix,perm,iprm,orgmcc,stbmcc,rocotime,verbose-1);
    if iprm /= null
     then Induced_Permutations.Permute(iprm.all,p);
    end if;
    Standard_Complex_Poly_Systems.Clear(sp);
  end Count_Roots;

  function Minimum ( a,b : natural64 ) return natural64 is

  -- DESCRIPTION :
  --   Returns the minimum of a and b.

  begin
    if a <= b
     then return a;
     else return b;
    end if;
  end Minimum;

  function Minimum ( a,b,c : natural64 ) return natural64 is

  -- DESCRIPTION :
  --   Returns the minimum of the three numbers a, b, and c.

     minbc : constant natural64 := Minimum(b,c);

  begin
    if a <= minbc
     then return a;
     else return minbc;
    end if;
  end Minimum;

  function Minimum ( a,b,c,d : natural64 ) return natural64 is

  -- DESCRIPTION :
  --   Returns the minimum of the four numbers.

    minbcd : constant natural64 := Minimum(b,c,d);

  begin
    if a <= minbcd
     then return a;
     else return minbcd;
    end if;
  end Minimum;

  procedure Construct_Start_System
              ( nt : in integer32;
                p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                d,bz,bs : in natural64;
                mv,smv : in natural32; z : in Partition;
                mix : in Link_to_Vector;
                stlb : in double_float; lifted : in Link_to_Array_of_Lists;
                orgmcc,stbmcc : in Mixed_Subdivision; roco : out natural64;
                q : in out Standard_Complex_Poly_Systems.Poly_Sys;
                qsols,qsols0 : in out Standard_Complex_Solutions.Solution_List;
                hocotime : out duration; verbose : in integer32 := 0 ) is

    timer : timing_widget;
    n : constant natural32 := natural32(p'length);
    nl : natural32;
    gap : constant natural32 := smv - mv;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Construct_Start_System 1,");
      put_line("for polynomial systems in double precision ...");
    end if;
    tstart(timer);
    if mv = 0 and smv = 0
     then roco := Minimum(d,bz,bs);
     else roco := Minimum(d,bz,bs,natural64(smv));
    end if;
    if gap = 0 then
      if roco = d then
        Start_System(p,q,qsols);
      elsif roco = bz then
        m_Homogeneous_Start_System(p,z,q,qsols);
      elsif roco = bs then
        Standard_Linear_Product_System.Init(n);
        Build_Random_Product_System(n);
        q := Standard_Linear_Product_System.Polynomial_System;
        Standard_Linear_Product_System.Solve(qsols,nl);
        Set_Structure.Clear;
        Standard_Linear_Product_System.Clear;
      else 
        Black_Box_Polyhedral_Continuation
          (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0,verbose-1);
      end if;
    else 
      Black_Box_Polyhedral_Continuation
        (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0,verbose-1);
    end if;
    tstop(timer);
    hocotime := Elapsed_User_Time(timer);
  exception
    when others => put_line("exception raised in construct start system");
                   put_line("the lifted supports : ");
                   Floating_Mixed_Subdivisions_io.put(lifted.all);
                   raise;
  end Construct_Start_System;

  procedure Construct_Start_System
              ( nt : in integer32;
                p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                d,bz,bs : in natural64;
                mv,smv : in natural32; z : in Partition;
                mix : in Link_to_Vector;
                stlb : in double_float; lifted : in Link_to_Array_of_Lists;
                orgmcc,stbmcc : in Mixed_Subdivision; roco : out natural64;
                q : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                qsols,qsols0 : in out DoblDobl_Complex_Solutions.Solution_List;
                hocotime : out duration; verbose : in integer32 := 0 ) is

    timer : timing_widget;
    n : constant natural32 := natural32(p'length);
    nl : natural32;
    gap : constant natural32 := smv - mv;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Construct_Start_System 2,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    tstart(timer);
    if mv = 0 and smv = 0
     then roco := Minimum(d,bz,bs);
     else roco := Minimum(d,bz,bs,natural64(smv));
    end if;
    if gap = 0 then
      if roco = d or roco = bz or roco = bs then
        declare
          sp : Standard_Complex_Poly_Systems.Poly_Sys(p'range)
             := DoblDobl_Complex_to_Standard_Poly_Sys(p);
          sq : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
          sqsols : Standard_Complex_Solutions.Solution_List;
        begin
          if roco = d then
            Start_System(sp,sq,sqsols);
          elsif roco = bz then
            m_Homogeneous_Start_System(sp,z,sq,sqsols);
          elsif roco = bs then
            Standard_Linear_Product_System.Init(n);
            Build_Random_Product_System(n);
            sq := Standard_Linear_Product_System.Polynomial_System;
            Standard_Linear_Product_System.Solve(sqsols,nl);
            Set_Structure.Clear;
            Standard_Linear_Product_System.Clear;
          end if;
          q := Standard_Poly_Sys_to_DoblDobl_Complex(sq);
          qsols := DoblDobl_Complex_Solutions.Create(sqsols);
          Standard_Complex_Poly_Systems.Clear(sp);
          Standard_Complex_Poly_Systems.Clear(sq);
          Standard_Complex_Solutions.Clear(sqsols);
        end;
      else 
        Black_Box_Polyhedral_Continuation
          (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0,verbose-1);
      end if;
    else 
      Black_Box_Polyhedral_Continuation
        (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0,verbose-1);
    end if;
    tstop(timer);
    hocotime := Elapsed_User_Time(timer);
  exception
    when others => put_line("exception raised in construct start system");
                   put_line("the lifted supports : ");
                   Floating_Mixed_Subdivisions_io.put(lifted.all);
                   raise;
  end Construct_Start_System;

  procedure Construct_Start_System
              ( nt : in integer32;
                p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                d,bz,bs : in natural64;
                mv,smv : in natural32; z : in Partition;
                mix : in Link_to_Vector;
                stlb : in double_float; lifted : in Link_to_Array_of_Lists;
                orgmcc,stbmcc : in Mixed_Subdivision; roco : out natural64;
                q : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                qsols,qsols0 : in out QuadDobl_Complex_Solutions.Solution_List;
                hocotime : out duration; verbose : in integer32 := 0 ) is

    timer : timing_widget;
    n : constant natural32 := natural32(p'length);
    nl : natural32;
    gap : constant natural32 := smv - mv;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Construct_Start_System 3,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    tstart(timer);
    if mv = 0 and smv = 0
     then roco := Minimum(d,bz,bs);
     else roco := Minimum(d,bz,bs,natural64(smv));
    end if;
    if gap = 0 then
      if roco = d or roco = bz or roco = bs then
        declare
          sp : Standard_Complex_Poly_Systems.Poly_Sys(p'range)
             := QuadDobl_Complex_to_Standard_Poly_Sys(p);
          sq : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
          sqsols : Standard_Complex_Solutions.Solution_List;
        begin
          if roco = d then
            Start_System(sp,sq,sqsols);
          elsif roco = bz then
            m_Homogeneous_Start_System(sp,z,sq,sqsols);
          elsif roco = bs then
            Standard_Linear_Product_System.Init(n);
            Build_Random_Product_System(n);
            sq := Standard_Linear_Product_System.Polynomial_System;
            Standard_Linear_Product_System.Solve(sqsols,nl);
            Set_Structure.Clear;
            Standard_Linear_Product_System.Clear;
          end if;
          q := Standard_Poly_Sys_to_QuadDobl_Complex(sq);
          qsols := QuadDobl_Complex_Solutions.Create(sqsols);
          Standard_Complex_Poly_Systems.Clear(sp);
          Standard_Complex_Poly_Systems.Clear(sq);
          Standard_Complex_Solutions.Clear(sqsols);
        end;
      else 
        Black_Box_Polyhedral_Continuation
          (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0,verbose-1);
      end if;
    else 
      Black_Box_Polyhedral_Continuation
        (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0,verbose-1);
    end if;
    tstop(timer);
    hocotime := Elapsed_User_Time(timer);
  exception
    when others => put_line("exception raised in construct start system");
                   put_line("the lifted supports : ");
                   Floating_Mixed_Subdivisions_io.put(lifted.all);
                   raise;
  end Construct_Start_System;

  procedure Construct_Start_System
              ( file : in file_type; nt : in integer32;
                p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                d,bz,bs : in natural64; mv,smv : in natural32;
                z : in Partition; mix : in Link_to_Vector;
                stlb : in double_float; lifted : in Link_to_Array_of_Lists;
                orgmcc,stbmcc : in Mixed_Subdivision; roco : out natural64;
                q : in out Standard_Complex_Poly_Systems.Poly_Sys;
                qsols,qsols0 : in out Standard_Complex_Solutions.Solution_List;
                hocotime : out duration; verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    timer : timing_widget;
    n : constant natural32 := natural32(p'length);
    rq : Prod_Sys(p'range);
    nl : natural32;
    prodform : boolean := false;
    wsols : Solution_List;
    gap : constant natural32 := smv - mv;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Construct_Start_System 4,");
      put_line("for polynomial systems in double precision ...");
    end if;
    new_line(file);
    tstart(timer);
    if mv = 0 and smv = 0
     then roco := Minimum(d,bz,bs);
     else roco := Minimum(d,bz,bs,natural64(smv));
    end if;
    if gap = 0 then
      if roco = d then
        put_line(file,"START SYSTEM BASED ON TOTAL DEGREE :");
        Start_System(p,q,qsols);
      elsif roco = bz then
        put(file,natural32(z'length),1);
        put_line(file,"-HOMOGENEOUS START SYSTEM :");
        m_Homogeneous_Start_System(p,z,q,rq,qsols);
        prodform := true;
      elsif roco = bs then
        put_line(file,"LINEAR-PRODUCT START SYSTEM : ");
        Standard_Linear_Product_System.Init(n);
        Build_Random_Product_System(n);
        q := Standard_Linear_Product_System.Polynomial_System;
        rq := Standard_Complex_Prod_Planes.Create;
        prodform := true;
        Standard_Linear_Product_System.Solve(qsols,nl);
        Set_Structure.Clear;
        Standard_Linear_Product_System.Clear;
      else
        put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
        Black_Box_Polyhedral_Continuation
          (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0,verbose-1);
      end if;
    else
      put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
      Black_Box_Polyhedral_Continuation
        (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0,verbose-1);
    end if;
    tstop(timer);
    new_line(file);
    if prodform
     then put_line(file,natural32(rq'last),rq); Deep_Clear(rq);
     else put_line(file,q);
    end if;
    new_line(file);
    put_line(file,"START SOLUTIONS : ");
    new_line(file);
    if Is_Null(qsols0) then
      put(file,Length_Of(qsols),natural32(q'last),qsols);
    else
      Push(qsols,wsols);
      Push(qsols0,wsols);
      put(file,Length_Of(wsols),natural32(q'last),wsols);
    end if;
    new_line(file);
    print_times(file,timer,"Construction of Start System");
    flush(file);
    hocotime := Elapsed_User_Time(timer);
  end Construct_Start_System;

  procedure Construct_Start_System
              ( file : in file_type; nt : in integer32;
                p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                d,bz,bs : in natural64; mv,smv : in natural32;
                z : in Partition; mix : in Link_to_Vector;
                stlb : in double_float; lifted : in Link_to_Array_of_Lists;
                orgmcc,stbmcc : in Mixed_Subdivision; roco : out natural64;
                q : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                qsols,qsols0 : in out DoblDobl_Complex_Solutions.Solution_List;
                hocotime : out duration; verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;

    timer : timing_widget;
    n : constant natural32 := natural32(p'length);
    rq : Prod_Sys(p'range);
    nl : natural32;
    prodform : boolean := false;
    wsols : Solution_List;
    gap : constant natural32 := smv - mv;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Construct_Start_System 5,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    new_line(file);
    tstart(timer);
    if mv = 0 and smv = 0
     then roco := Minimum(d,bz,bs);
     else roco := Minimum(d,bz,bs,natural64(smv));
    end if;
    if gap = 0 then
      if roco = d or roco = bz or roco = bs then
        declare
          sp : Standard_Complex_Poly_Systems.Poly_Sys(p'range)
             := DoblDobl_Complex_to_Standard_Poly_Sys(p);
          sq : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
          sqsols : Standard_Complex_Solutions.Solution_List;
        begin
          if roco = d then
            put_line(file,"START SYSTEM BASED ON TOTAL DEGREE :");
            Start_System(sp,sq,sqsols);
          elsif roco = bz then
            put(file,natural32(z'length),1);
            put_line(file,"-HOMOGENEOUS START SYSTEM :");
            m_Homogeneous_Start_System(sp,z,sq,rq,sqsols);
            prodform := true;
          elsif roco = bs then
            put_line(file,"LINEAR-PRODUCT START SYSTEM : ");
            Standard_Linear_Product_System.Init(n);
            Build_Random_Product_System(n);
            sq := Standard_Linear_Product_System.Polynomial_System;
            rq := Standard_Complex_Prod_Planes.Create;
            prodform := true;
            Standard_Linear_Product_System.Solve(sqsols,nl);
            Set_Structure.Clear;
            Standard_Linear_Product_System.Clear;
          end if;
          q := Standard_Poly_Sys_to_DoblDobl_Complex(sq);
          qsols := DoblDobl_Complex_Solutions.Create(sqsols);
          Standard_Complex_Poly_Systems.Clear(sp);
          Standard_Complex_Poly_Systems.Clear(sq);
          Standard_Complex_Solutions.Clear(sqsols);
        end;
      else
        put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
        Black_Box_Polyhedral_Continuation
          (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0,verbose-1);
      end if;
    else
      put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
      Black_Box_Polyhedral_Continuation
        (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0,verbose-1);
    end if;
    tstop(timer);
    new_line(file);
    if prodform
     then put_line(file,natural32(rq'last),rq); Deep_Clear(rq);
     else put_line(file,q);
    end if;
    new_line(file);
    put_line(file,"START SOLUTIONS : ");
    new_line(file);
    if Is_Null(qsols0) then
      put(file,Length_Of(qsols),natural32(q'last),qsols);
    else
      Push(qsols,wsols);
      Push(qsols0,wsols);
      put(file,Length_Of(wsols),natural32(q'last),wsols);
    end if;
    new_line(file);
    print_times(file,timer,"Construction of Start System");
    flush(file);
    hocotime := Elapsed_User_Time(timer);
  end Construct_Start_System;

  procedure Construct_Start_System
              ( file : in file_type; nt : in integer32;
                p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                d,bz,bs : in natural64; mv,smv : in natural32;
                z : in Partition; mix : in Link_to_Vector;
                stlb : in double_float; lifted : in Link_to_Array_of_Lists;
                orgmcc,stbmcc : in Mixed_Subdivision; roco : out natural64;
                q : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                qsols,qsols0 : in out QuadDobl_Complex_Solutions.Solution_List;
                hocotime : out duration; verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;

    timer : timing_widget;
    n : constant natural32 := natural32(p'length);
    rq : Prod_Sys(p'range);
    nl : natural32;
    prodform : boolean := false;
    wsols : Solution_List;
    gap : constant natural32 := smv - mv;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Construct_Start_System 6,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    new_line(file);
    tstart(timer);
    if mv = 0 and smv = 0
     then roco := Minimum(d,bz,bs);
     else roco := Minimum(d,bz,bs,natural64(smv));
    end if;
    if gap = 0 then
      if roco = d or roco = bz or roco = bs then
        declare
          sp : Standard_Complex_Poly_Systems.Poly_Sys(p'range)
             := QuadDobl_Complex_to_Standard_Poly_Sys(p);
          sq : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
          sqsols : Standard_Complex_Solutions.Solution_List;
        begin
          if roco = d then
            put_line(file,"START SYSTEM BASED ON TOTAL DEGREE :");
            Start_System(sp,sq,sqsols);
          elsif roco = bz then
            put(file,natural32(z'length),1);
            put_line(file,"-HOMOGENEOUS START SYSTEM :");
            m_Homogeneous_Start_System(sp,z,sq,rq,sqsols);
            prodform := true;
          elsif roco = bs then
            put_line(file,"LINEAR-PRODUCT START SYSTEM : ");
            Standard_Linear_Product_System.Init(n);
            Build_Random_Product_System(n);
            sq := Standard_Linear_Product_System.Polynomial_System;
            rq := Standard_Complex_Prod_Planes.Create;
            prodform := true;
            Standard_Linear_Product_System.Solve(sqsols,nl);
            Set_Structure.Clear;
            Standard_Linear_Product_System.Clear;
          end if;
          q := Standard_Poly_Sys_to_QuadDobl_Complex(sq);
          qsols := QuadDobl_Complex_Solutions.Create(sqsols);
          Standard_Complex_Poly_Systems.Clear(sp);
          Standard_Complex_Poly_Systems.Clear(sq);
          Standard_Complex_Solutions.Clear(sqsols);
        end;
      else
        put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
        Black_Box_Polyhedral_Continuation
          (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0,verbose-1);
      end if;
    else
      put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
      Black_Box_Polyhedral_Continuation
        (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0,verbose-1);
    end if;
    tstop(timer);
    new_line(file);
    if prodform
     then put_line(file,natural32(rq'last),rq); Deep_Clear(rq);
     else put_line(file,q);
    end if;
    new_line(file);
    put_line(file,"START SOLUTIONS : ");
    new_line(file);
    if Is_Null(qsols0) then
      put(file,Length_Of(qsols),natural32(q'last),qsols);
    else
      Push(qsols,wsols);
      Push(qsols0,wsols);
      put(file,Length_Of(wsols),natural32(q'last),wsols);
    end if;
    new_line(file);
    print_times(file,timer,"Construction of Start System");
    flush(file);
    hocotime := Elapsed_User_Time(timer);
  end Construct_Start_System;

  procedure Black_Box_Root_Counting
               ( nt : in integer32; silent : in boolean;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    d,bz,bs,wrc : natural64;
    mv,smv : natural32;
    z : partition(1..natural32(p'last));
    nz : natural32;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    orgmcc,stbmcc : Mixed_Subdivision;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists;
    no_mv : constant boolean := deg or (natural32(p'last) > chicken_mv);
    mptdeg : Natural_Number;

    use Root_Counters_Output;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Black_Box_Root_Counting 1,");
      put_line("for polynomial systems in double precision ...");
    end if;
    Count_Roots(p,deg,d,mptdeg,bz,bs,mv,smv,z,nz,
                stlb,lifsup,mix,perm,iprm,orgmcc,stbmcc,rocotime,verbose-1);
    if (not silent) or (verbose > 0)
     then Write_Root_Counts(standard_output,no_mv,d,mptdeg,nz,bz,bs,mv,smv,z);
    end if;
    Construct_Start_System
      (nt,p,d,bz,bs,mv,smv,z(1..nz),mix,stlb,lifsup,orgmcc,stbmcc,wrc,
       q,qsols,qsols0,hocotime,verbose-1);
    rc := natural32(wrc);
    Clear(z);
    Clear(mix);
    Deep_Clear(lifsup);
    Deep_Clear(orgmcc);
    Deep_Clear(stbmcc);
  exception
    when others => put_line("exception raised in black box root counting");
                   raise;
  end Black_Box_Root_Counting;

  procedure Black_Box_Root_Counting
               ( nt : in integer32; silent : in boolean;
                 p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out DoblDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    d,bz,bs,wrc : natural64;
    mv,smv : natural32;
    z : partition(1..natural32(p'last));
    nz : natural32;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    orgmcc,stbmcc : Mixed_Subdivision;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists;
    no_mv : constant boolean := deg or (natural32(p'last) > chicken_mv);
    mptdeg : Natural_Number;

    use Root_Counters_Output;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Black_Box_Root_Counting 2,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    Count_Roots(p,deg,d,mptdeg,bz,bs,mv,smv,z,nz,
                stlb,lifsup,mix,perm,iprm,orgmcc,stbmcc,rocotime,verbose-1);
    if (not silent) or (verbose > 0)
     then Write_Root_Counts(standard_output,no_mv,d,mptdeg,nz,bz,bs,mv,smv,z);
    end if;
    Construct_Start_System
      (nt,p,d,bz,bs,mv,smv,z(1..nz),mix,stlb,lifsup,orgmcc,stbmcc,wrc,
       q,qsols,qsols0,hocotime,verbose-1);
    rc := natural32(wrc);
    Clear(z);
    Clear(mix);
    Deep_Clear(lifsup);
    Deep_Clear(orgmcc);
    Deep_Clear(stbmcc);
  exception
    when others => put_line("exception raised in black box root counting");
                   raise;
  end Black_Box_Root_Counting;

  procedure Black_Box_Root_Counting
               ( nt : in integer32; silent : in boolean;
                 p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out QuadDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    d,bz,bs,wrc : natural64;
    mv,smv : natural32;
    z : partition(1..natural32(p'last));
    nz : natural32;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    orgmcc,stbmcc : Mixed_Subdivision;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists;
    no_mv : constant boolean := deg or (natural32(p'last) > chicken_mv);
    mptdeg : Natural_Number;

    use Root_Counters_Output;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Black_Box_Root_Counting 3,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    Count_Roots(p,deg,d,mptdeg,bz,bs,mv,smv,z,nz,
                stlb,lifsup,mix,perm,iprm,orgmcc,stbmcc,rocotime,verbose-1);
    if (not silent) or (verbose > 0)
     then Write_Root_Counts(standard_output,no_mv,d,mptdeg,nz,bz,bs,mv,smv,z);
    end if;
    Construct_Start_System
      (nt,p,d,bz,bs,mv,smv,z(1..nz),mix,stlb,lifsup,orgmcc,stbmcc,wrc,
       q,qsols,qsols0,hocotime,verbose-1);
    rc := natural32(wrc);
    Clear(z);
    Clear(mix);
    Deep_Clear(lifsup);
    Deep_Clear(orgmcc);
    Deep_Clear(stbmcc);
  exception
    when others => put_line("exception raised in black box root counting");
                   raise;
  end Black_Box_Root_Counting;

  procedure Black_Box_Root_Counting 
               ( nt : in integer32;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 rocos : out Link_to_String;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    d,bz,bs,wrc : natural64;
    mv,smv : natural32;
    z : partition(1..natural32(p'last));
    nz : natural32;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    orgmcc,stbmcc : Mixed_Subdivision;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists;
    no_mv : constant boolean := deg or (natural32(p'last) > chicken_mv);
    mptdeg : Natural_Number;

    use Root_Counters_Output;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Black_Box_Root_Counting 4,");
      put_line("for polynomial systems in double precision ...");
    end if;
    Count_Roots(p,deg,d,mptdeg,bz,bs,mv,smv,z,nz,
                stlb,lifsup,mix,perm,iprm,orgmcc,stbmcc,rocotime,verbose-1);
    declare
      rcs : constant string
          := Root_Counts_to_String(no_mv,d,mptdeg,nz,bz,bs,mv,smv,z);
    begin
      rocos := new string'(rcs);
      if verbose > 0
       then put_line(rcs);
      end if;
    end;
    Construct_Start_System
      (nt,p,d,bz,bs,mv,smv,z(1..nz),mix,stlb,lifsup,orgmcc,stbmcc,wrc,
       q,qsols,qsols0,hocotime,verbose-1);
    rc := natural32(wrc);
    Clear(z);
    Clear(mix);
    Deep_Clear(lifsup);
    Deep_Clear(orgmcc);
    Deep_Clear(stbmcc);
  exception
    when others => put_line("exception raised in black box root counting");
                   raise;
  end Black_Box_Root_Counting;

  procedure Black_Box_Root_Counting
               ( nt : in integer32;
                 p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 rocos : out Link_to_String;
                 q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out DoblDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    d,bz,bs,wrc : natural64;
    mv,smv : natural32;
    z : partition(1..natural32(p'last));
    nz : natural32;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    orgmcc,stbmcc : Mixed_Subdivision;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists;
    no_mv : constant boolean := deg or (natural32(p'last) > chicken_mv);
    mptdeg : Natural_Number;

    use Root_Counters_Output;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Black_Box_Root_Counting 5,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    Count_Roots(p,deg,d,mptdeg,bz,bs,mv,smv,z,nz,
                stlb,lifsup,mix,perm,iprm,orgmcc,stbmcc,rocotime,verbose-1);
    declare
      rcs : constant string
          := Root_Counts_to_String(no_mv,d,mptdeg,nz,bz,bs,mv,smv,z);
    begin
      rocos := new string'(rcs);
      if verbose > 0
       then put_line(rcs);
      end if;
    end;
    Construct_Start_System
      (nt,p,d,bz,bs,mv,smv,z(1..nz),mix,stlb,lifsup,orgmcc,stbmcc,wrc,
       q,qsols,qsols0,hocotime,verbose-1);
    rc := natural32(wrc);
    Clear(z);
    Clear(mix);
    Deep_Clear(lifsup);
    Deep_Clear(orgmcc);
    Deep_Clear(stbmcc);
  exception
    when others => put_line("exception raised in black box root counting");
                   raise;
  end Black_Box_Root_Counting;

  procedure Black_Box_Root_Counting
               ( nt : in integer32;
                 p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 rocos : out Link_to_String;
                 q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out QuadDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    d,bz,bs,wrc : natural64;
    mv,smv : natural32;
    z : partition(1..natural32(p'last));
    nz : natural32;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    orgmcc,stbmcc : Mixed_Subdivision;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists;
    no_mv : constant boolean := deg or (natural32(p'last) > chicken_mv);
    mptdeg : Natural_Number;

    use Root_Counters_Output;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Black_Box_Root_Counting 6,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    Count_Roots(p,deg,d,mptdeg,bz,bs,mv,smv,z,nz,
                stlb,lifsup,mix,perm,iprm,orgmcc,stbmcc,rocotime,verbose-1);
    declare
      rcs : constant string
          := Root_Counts_to_String(no_mv,d,mptdeg,nz,bz,bs,mv,smv,z);
    begin
      rocos := new string'(rcs);
      if verbose > 0
       then put_line(rcs);
      end if;
    end;
    Construct_Start_System
      (nt,p,d,bz,bs,mv,smv,z(1..nz),mix,stlb,lifsup,orgmcc,stbmcc,wrc,
       q,qsols,qsols0,hocotime,verbose-1);
    rc := natural32(wrc);
    Clear(z);
    Clear(mix);
    Deep_Clear(lifsup);
    Deep_Clear(orgmcc);
    Deep_Clear(stbmcc);
  exception
    when others => put_line("exception raised in black box root counting");
                   raise;
  end Black_Box_Root_Counting;

  procedure Black_Box_Root_Counting
               ( file : in file_type; nt : in integer32;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    d,bz,bs,wrc : natural64;
    mv,smv : natural32;
    z : partition(1..natural32(p'last));
    nz : natural32;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    orgmcc,stbmcc : Mixed_Subdivision;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists;
    mptdeg : Natural_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Black_Box_Root_Counting 7,");
      put_line("for polynomial systems in double precision ...");
    end if;
    Count_Roots(file,p,deg,d,mptdeg,bz,bs,mv,smv,z,nz,
                stlb,lifsup,mix,perm,iprm,orgmcc,stbmcc,rocotime,verbose-1);
    Construct_Start_System
      (file,nt,p,d,bz,bs,mv,smv,z(1..nz),mix,stlb,lifsup,orgmcc,stbmcc,
       wrc,q,qsols,qsols0,hocotime,verbose-1);
    rc := natural32(wrc);
    Clear(z);
    Clear(mix);
    Deep_Clear(lifsup);
    Deep_Clear(orgmcc);
    Deep_Clear(stbmcc);
  end Black_Box_Root_Counting;

  procedure Black_Box_Root_Counting
               ( file : in file_type; nt : in integer32;
                 p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out DoblDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    d,bz,bs,wrc : natural64;
    mv,smv : natural32;
    z : partition(1..natural32(p'last));
    nz : natural32;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    orgmcc,stbmcc : Mixed_Subdivision;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists;
    mptdeg : Natural_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Black_Box_Root_Counting 8,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    Count_Roots(file,p,deg,d,mptdeg,bz,bs,mv,smv,z,nz,
                stlb,lifsup,mix,perm,iprm,orgmcc,stbmcc,rocotime,verbose-1);
    Construct_Start_System
      (file,nt,p,d,bz,bs,mv,smv,z(1..nz),mix,stlb,lifsup,orgmcc,stbmcc,
       wrc,q,qsols,qsols0,hocotime,verbose-1);
    rc := natural32(wrc);
    Clear(z);
    Clear(mix);
    Deep_Clear(lifsup);
    Deep_Clear(orgmcc);
    Deep_Clear(stbmcc);
  end Black_Box_Root_Counting;

  procedure Black_Box_Root_Counting
               ( file : in file_type; nt : in integer32;
                 p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out QuadDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    d,bz,bs,wrc : natural64;
    mv,smv : natural32;
    z : partition(1..natural32(p'last));
    nz : natural32;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    orgmcc,stbmcc : Mixed_Subdivision;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists;
    mptdeg : Natural_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Black_Box_Root_Counting 9,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    Count_Roots(file,p,deg,d,mptdeg,bz,bs,mv,smv,z,nz,
                stlb,lifsup,mix,perm,iprm,orgmcc,stbmcc,rocotime,verbose-1);
    Construct_Start_System
      (file,nt,p,d,bz,bs,mv,smv,z(1..nz),mix,stlb,lifsup,orgmcc,stbmcc,
       wrc,q,qsols,qsols0,hocotime,verbose-1);
    rc := natural32(wrc);
    Clear(z);
    Clear(mix);
    Deep_Clear(lifsup);
    Deep_Clear(orgmcc);
    Deep_Clear(stbmcc);
  end Black_Box_Root_Counting;

-- PIPELINED BLACK BOX ROOT COUNTING :

  procedure Pipelined_Stable_Continuation
              ( silent : in boolean; r : in integer32;
                mtype : in Standard_Integer_Vectors.Link_to_Vector;
                stlb : in double_float;
                lifsup: in Link_to_Array_of_Lists;
                mcc : in Mixed_Subdivision; tmv : in natural32;
                lq : in Standard_Complex_Laur_Systems.Laur_Sys;
                mv,smv : out natural32;
                qsols0 : out Standard_Complex_Solutions.Solution_List ) is

   -- start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
   -- ended_moment : Ada.Calendar.Time;
    dim : constant integer32 := lq'last;
    orgmcc,stbmcc : Mixed_Subdivision;
    orgcnt,stbcnt : natural32;

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

  begin
    Split_Original_Cells(mcc,stlb,orgmcc,stbmcc,orgcnt,stbcnt);
    if stbcnt = 0 then
      if not silent then
        put("mixed volume : "); put(tmv,1); new_line;
        put("stable mixed volume : "); put(tmv,1); new_line;
      end if;
      mv := tmv; smv := tmv;
    else
      declare
        mix : constant Standard_Integer_Vectors.Vector
            := Pipelined_Labeled_Cells.Mixture(r,mtype);
        mtp : Standard_Integer_Vectors.Link_to_Vector;
      begin
        Mixed_Volume_Computation.Mixed_Volume(dim,mix,orgmcc,mv);
        if not silent then
          put("mixed volume : "); put(mv,1); new_line;
        end if;
        Mixed_Volume_Computation.Mixed_Volume(dim,mix,stbmcc,smv);
        smv := smv + mv;
        if not silent then
          put("stable mixed volume : "); put(smv,1); new_line;
        end if;
        if not Is_Null(stbmcc) then
          mtp := new Standard_Integer_Vectors.Vector'(mix);
          Induced_Permutations.Remove_Artificial_Origin(lifsup.all,stlb);
          Silent_Polyhedral_Continuation
            (lq,stlb,mtp,lifsup.all,stbmcc,qsols0);
          Set_Continuation_Parameter(qsols0,Create(0.0));
        end if;
      end;
    end if;
   -- ended_moment := Ada.Calendar.Clock;
   -- if not silent
   --  then Write_Elapsed_Time(standard_output,start_moment,ended_moment);
   -- end if;
  end Pipelined_Stable_Continuation;

  procedure Pipelined_Stable_Continuation
              ( silent : in boolean; r : in integer32;
                mtype : in Standard_Integer_Vectors.Link_to_Vector;
                stlb : in double_float;
                lifsup: in Link_to_Array_of_Lists;
                mcc : in Mixed_Subdivision; tmv : in natural32;
                lq : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mv,smv : out natural32;
                qsols0 : out DoblDobl_Complex_Solutions.Solution_List ) is

   -- start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
   -- ended_moment : Ada.Calendar.Time;
    dim : constant integer32 := lq'last;
    orgmcc,stbmcc : Mixed_Subdivision;
    orgcnt,stbcnt : natural32;
    zero : constant double_double := create(0.0);

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

  begin
    Split_Original_Cells(mcc,stlb,orgmcc,stbmcc,orgcnt,stbcnt);
    if stbcnt = 0 then
      if not silent then
        put("mixed volume : "); put(tmv,1); new_line;
        put("stable mixed volume : "); put(tmv,1); new_line;
      end if;
      mv := tmv; smv := tmv;
    else
      declare
        mix : constant Standard_Integer_Vectors.Vector
            := Pipelined_Labeled_Cells.Mixture(r,mtype);
        mtp : Standard_Integer_Vectors.Link_to_Vector;
      begin
        Mixed_Volume_Computation.Mixed_Volume(dim,mix,orgmcc,mv);
        if not silent then
          put("mixed volume : "); put(mv,1); new_line;
        end if;
        Mixed_Volume_Computation.Mixed_Volume(dim,mix,stbmcc,smv);
        smv := smv + mv;
        if not silent then
          put("stable mixed volume : "); put(smv,1); new_line;
        end if;
        if not Is_Null(stbmcc) then
          mtp := new Standard_Integer_Vectors.Vector'(mix);
          Induced_Permutations.Remove_Artificial_Origin(lifsup.all,stlb);
          Silent_Polyhedral_Continuation
            (lq,stlb,mtp,lifsup.all,stbmcc,qsols0);
          Set_Continuation_Parameter(qsols0,Create(zero));
        end if;
      end;
    end if;
   -- ended_moment := Ada.Calendar.Clock;
   -- if not silent
   --  then Write_Elapsed_Time(standard_output,start_moment,ended_moment);
   -- end if;
  end Pipelined_Stable_Continuation;

  procedure Pipelined_Stable_Continuation
              ( silent : in boolean; r : in integer32;
                mtype : in Standard_Integer_Vectors.Link_to_Vector;
                stlb : in double_float;
                lifsup: in Link_to_Array_of_Lists;
                mcc : in Mixed_Subdivision; tmv : in natural32;
                lq : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mv,smv : out natural32;
                qsols0 : out QuadDobl_Complex_Solutions.Solution_List ) is

   -- start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
   -- ended_moment : Ada.Calendar.Time;
    dim : constant integer32 := lq'last;
    orgmcc,stbmcc : Mixed_Subdivision;
    orgcnt,stbcnt : natural32;
    zero : constant quad_double := create(0.0);

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

  begin
    Split_Original_Cells(mcc,stlb,orgmcc,stbmcc,orgcnt,stbcnt);
    if stbcnt = 0 then
      if not silent then
        put("mixed volume : "); put(tmv,1); new_line;
        put("stable mixed volume : "); put(tmv,1); new_line;
      end if;
      mv := tmv; smv := tmv;
    else
      declare
        mix : constant Standard_Integer_Vectors.Vector
            := Pipelined_Labeled_Cells.Mixture(r,mtype);
        mtp : Standard_Integer_Vectors.Link_to_Vector;
      begin
        Mixed_Volume_Computation.Mixed_Volume(dim,mix,orgmcc,mv);
        if not silent then
          put("mixed volume : "); put(mv,1); new_line;
        end if;
        Mixed_Volume_Computation.Mixed_Volume(dim,mix,stbmcc,smv);
        smv := smv + mv;
        if not silent then
          put("stable mixed volume : "); put(smv,1); new_line;
        end if;
        if not Is_Null(stbmcc) then
          mtp := new Standard_Integer_Vectors.Vector'(mix);
          Induced_Permutations.Remove_Artificial_Origin(lifsup.all,stlb);
          Silent_Polyhedral_Continuation
            (lq,stlb,mtp,lifsup.all,stbmcc,qsols0);
          Set_Continuation_Parameter(qsols0,Create(zero));
        end if;
      end;
    end if;
   -- ended_moment := Ada.Calendar.Clock;
   -- if not silent
   --  then Write_Elapsed_Time(standard_output,start_moment,ended_moment);
   -- end if;
  end Pipelined_Stable_Continuation;

  procedure Pipelined_Root_Counting
               ( nt : in integer32; silent : in boolean;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    d : natural64;
    mptdeg : Natural_Number;
    mv,smv,tmv : natural32;
    r : integer32;
    mtype,perm : Standard_Integer_Vectors.Link_to_Vector;
    mcc : Mixed_Subdivision;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists;
    no_mv : constant boolean := deg or (natural32(p'last) > chicken_mv);
    timer : Timing_Widget;
    lq : Standard_Complex_Laur_Systems.Laur_Sys(p'range);

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;
    use Root_Counters_Output;
    use Pipelined_Polyhedral_Drivers;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Pipelined_Count_Rooting 1,");
      put_line("for polynomial systems in standard double precision ...");
    end if;
    if nt < 2 or no_mv then
      Black_Box_Root_Counting
        (nt,silent,p,deg,rc,q,qsols,qsols0,rocotime,hocotime);
    else
      Compute_Total_Degree(p,d,mptdeg);
      if not silent
       then Write_Total_Degree(standard_output,d,mptdeg);
      end if;
      stlb := Floating_Lifting_Functions.Lifting_Bound(p);
      tstart(timer);
      Pipelined_Polyhedral_Homotopies
        (nt,true,stlb,p,r,mtype,perm,lifsup,mcc,tmv,lq,q,qsols);
      Apply_Induced_Permutation(p,stlb,r,perm,mtype,mcc);
      Set_Continuation_Parameter(qsols,Create(0.0));
      Pipelined_Stable_Continuation
        (silent,r,mtype,stlb,lifsup,mcc,tmv,lq,mv,smv,qsols0);
      rc := smv;
      tstop(timer);
      rocotime := 0.0;
      hocotime := Elapsed_User_Time(timer);
    end if;
    ended_moment := Ada.Calendar.Clock;
    if not silent then
      new_line;
      print_times(standard_output,timer,"pipelined polyhedral homotopies");
      new_line;
      Write_Elapsed_Time(standard_output,start_moment,ended_moment);
    end if;
  end Pipelined_Root_Counting;

  procedure Pipelined_Root_Counting
               ( nt : in integer32; silent : in boolean;
                 p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out DoblDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    d : natural64;
    mptdeg : Natural_Number;
    mv,smv,tmv : natural32;
    r : integer32;
    mtype,perm : Standard_Integer_Vectors.Link_to_Vector;
    mcc : Mixed_Subdivision;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists;
    no_mv : constant boolean := deg or (natural32(p'last) > chicken_mv);
    timer : Timing_Widget;
    lq : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    zero : constant double_double := create(0.0);

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;
    use Root_Counters_Output;
    use Pipelined_Polyhedral_Drivers;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Pipelined_Count_Rooting 2,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    if nt < 2 or no_mv then
      Black_Box_Root_Counting
        (nt,silent,p,deg,rc,q,qsols,qsols0,rocotime,hocotime);
    else
      Compute_Total_Degree(p,d,mptdeg);
      if not silent
       then Write_Total_Degree(standard_output,d,mptdeg);
      end if;
      stlb := Floating_Lifting_Functions.Lifting_Bound(p);
      tstart(timer);
      Pipelined_Polyhedral_Homotopies
        (nt,true,stlb,p,r,mtype,perm,lifsup,mcc,tmv,lq,q,qsols);
      Apply_Induced_Permutation(p,stlb,r,perm,mtype,mcc);
      Set_Continuation_Parameter(qsols,Create(zero));
      Pipelined_Stable_Continuation
        (silent,r,mtype,stlb,lifsup,mcc,tmv,lq,mv,smv,qsols0);
      rc := smv;
      tstop(timer);
      rocotime := 0.0;
      hocotime := Elapsed_User_Time(timer);
    end if;
    ended_moment := Ada.Calendar.Clock;
    if not silent then
      new_line;
      print_times(standard_output,timer,"pipelined polyhedral homotopies");
      new_line;
      Write_Elapsed_Time(standard_output,start_moment,ended_moment);
    end if;
  end Pipelined_Root_Counting;

  procedure Pipelined_Root_Counting
               ( nt : in integer32; silent : in boolean;
                 p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out QuadDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    d : natural64;
    mptdeg : Natural_Number;
    mv,smv,tmv : natural32;
    r : integer32;
    mtype,perm : Standard_Integer_Vectors.Link_to_Vector;
    mcc : Mixed_Subdivision;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists;
    no_mv : constant boolean := deg or (natural32(p'last) > chicken_mv);
    timer : Timing_Widget;
    lq : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    zero : constant quad_double := create(0.0);

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;
    use Root_Counters_Output;
    use Pipelined_Polyhedral_Drivers;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Pipelined_Count_Rooting 3,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    if nt < 2 or no_mv then
      Black_Box_Root_Counting
        (nt,silent,p,deg,rc,q,qsols,qsols0,rocotime,hocotime);
    else
      Compute_Total_Degree(p,d,mptdeg);
      if not silent
       then Write_Total_Degree(standard_output,d,mptdeg);
      end if;
      stlb := Floating_Lifting_Functions.Lifting_Bound(p);
      tstart(timer);
      Pipelined_Polyhedral_Homotopies
        (nt,true,stlb,p,r,mtype,perm,lifsup,mcc,tmv,lq,q,qsols);
      Apply_Induced_Permutation(p,stlb,r,perm,mtype,mcc);
      Set_Continuation_Parameter(qsols,Create(zero));
      Pipelined_Stable_Continuation
        (silent,r,mtype,stlb,lifsup,mcc,tmv,lq,mv,smv,qsols0);
      rc := smv;
      tstop(timer);
      rocotime := 0.0;
      hocotime := Elapsed_User_Time(timer);
    end if;
    ended_moment := Ada.Calendar.Clock;
    if not silent then
      new_line;
      print_times(standard_output,timer,"pipelined polyhedral homotopies");
      new_line;
      Write_Elapsed_Time(standard_output,start_moment,ended_moment);
    end if;
  end Pipelined_Root_Counting;

  procedure Pipelined_Root_Counting 
               ( nt : in integer32;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 rocos : out Link_to_String;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    d : natural64;
    mptdeg : Natural_Number;
    mv,smv,tmv : natural32;
    r : integer32;
    mtype,perm : Standard_Integer_Vectors.Link_to_Vector;
    mcc : Mixed_Subdivision;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists;
    no_mv : constant boolean := deg or (natural32(p'last) > chicken_mv);
    timer : Timing_Widget;
    lq : Standard_Complex_Laur_Systems.Laur_Sys(p'range);

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;
    use Root_Counters_Output;
    use Pipelined_Polyhedral_Drivers;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Pipelined_Count_Rooting 4,");
      put_line("for polynomial systems in standard double precision ...");
    end if;
    if nt < 2 or no_mv then
      Black_Box_Root_Counting
        (nt,p,deg,rc,rocos,q,qsols,qsols0,rocotime,hocotime);
    else
      Compute_Total_Degree(p,d,mptdeg);
      stlb := Floating_Lifting_Functions.Lifting_Bound(p);
      tstart(timer);
      Pipelined_Polyhedral_Homotopies
        (nt,true,stlb,p,r,mtype,perm,lifsup,mcc,tmv,lq,q,qsols);
      Set_Continuation_Parameter(qsols0,Create(0.0));
      Pipelined_Stable_Continuation
        (true,r,mtype,stlb,lifsup,mcc,tmv,lq,mv,smv,qsols0);
      Apply_Induced_Permutation(p,stlb,r,perm,mtype,mcc);
      rc := smv;
      tstop(timer);
      rocotime := 0.0;
      hocotime := Elapsed_User_Time(timer);
      declare
        rcs : constant string := Mixed_Volumes_to_String(d,mv,smv);
      begin
        rocos := new string'(rcs);
      end;
    end if;
  end Pipelined_Root_Counting;

  procedure Pipelined_Root_Counting 
               ( nt : in integer32;
                 p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 rocos : out Link_to_String;
                 q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out DoblDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    d : natural64;
    mptdeg : Natural_Number;
    mv,smv,tmv : natural32;
    r : integer32;
    mtype,perm : Standard_Integer_Vectors.Link_to_Vector;
    mcc : Mixed_Subdivision;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists;
    no_mv : constant boolean := deg or (natural32(p'last) > chicken_mv);
    timer : Timing_Widget;
    zero : constant double_double := create(0.0);
    lq : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;
    use Root_Counters_Output;
    use Pipelined_Polyhedral_Drivers;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Pipelined_Count_Rooting 5,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    if nt < 2 or no_mv then
      Black_Box_Root_Counting
        (nt,p,deg,rc,rocos,q,qsols,qsols0,rocotime,hocotime);
    else
      Compute_Total_Degree(p,d,mptdeg);
      stlb := Floating_Lifting_Functions.Lifting_Bound(p);
      tstart(timer);
      Pipelined_Polyhedral_Homotopies
        (nt,true,stlb,p,r,mtype,perm,lifsup,mcc,tmv,lq,q,qsols);
      Set_Continuation_Parameter(qsols,Create(zero));
      Pipelined_Stable_Continuation
        (true,r,mtype,stlb,lifsup,mcc,tmv,lq,mv,smv,qsols0);
      Apply_Induced_Permutation(p,stlb,r,perm,mtype,mcc);
      rc := smv;
      tstop(timer);
      rocotime := 0.0;
      hocotime := Elapsed_User_Time(timer);
      declare
        rcs : constant string := Mixed_Volumes_to_String(d,mv,smv);
      begin
        rocos := new string'(rcs);
      end;
    end if;
  end Pipelined_Root_Counting;

  procedure Pipelined_Root_Counting 
               ( nt : in integer32;
                 p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 rocos : out Link_to_String;
                 q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out QuadDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    d : natural64;
    mptdeg : Natural_Number;
    mv,smv,tmv : natural32;
    r : integer32;
    mtype,perm : Standard_Integer_Vectors.Link_to_Vector;
    mcc : Mixed_Subdivision;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists;
    no_mv : constant boolean := deg or (natural32(p'last) > chicken_mv);
    timer : Timing_Widget;
    lq : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    zero : constant quad_double := create(0.0);

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;
    use Root_Counters_Output;
    use Pipelined_Polyhedral_Drivers;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Pipelined_Count_Rooting 6,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    if nt < 2 or no_mv then
      Black_Box_Root_Counting
        (nt,p,deg,rc,rocos,q,qsols,qsols0,rocotime,hocotime);
    else
      Compute_Total_Degree(p,d,mptdeg);
      stlb := Floating_Lifting_Functions.Lifting_Bound(p);
      tstart(timer);
      Pipelined_Polyhedral_Homotopies
        (nt,true,stlb,p,r,mtype,perm,lifsup,mcc,tmv,lq,q,qsols);
      Apply_Induced_Permutation(p,stlb,r,perm,mtype,mcc);
      Set_Continuation_Parameter(qsols,Create(zero));
      Pipelined_Stable_Continuation
        (true,r,mtype,stlb,lifsup,mcc,tmv,lq,mv,smv,qsols0);
      rc := smv;
      tstop(timer);
      rocotime := 0.0;
      hocotime := Elapsed_User_Time(timer);
      declare
        rcs : constant string := Mixed_Volumes_to_String(d,mv,smv);
      begin
        rocos := new string'(rcs);
      end;
      rc := mv;
    end if;
  end Pipelined_Root_Counting;

  procedure Pipelined_Root_Counting
               ( file : in file_type; nt : in integer32;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    d : natural64;
    mptdeg : Natural_Number;
    mv,smv,tmv : natural32;
    r : integer32;
    mtype,perm : Standard_Integer_Vectors.Link_to_Vector;
    mcc : Mixed_Subdivision;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists;
    timer : Timing_Widget;
    lq : Standard_Complex_Laur_Systems.Laur_Sys(p'range);

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;
    use Root_Counters_Output;
    use Pipelined_Polyhedral_Drivers;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Pipelined_Count_Rooting 7,");
      put_line("for polynomial systems in standard double precision ...");
    end if;
    if nt < 2 or deg then
      Black_Box_Root_Counting
        (file,nt,p,deg,rc,q,qsols,qsols0,rocotime,hocotime);
    else
      Compute_Total_Degree(p,d,mptdeg);
      Write_Total_Degree(file,d,mptdeg);
      flush(file);
      stlb := Floating_Lifting_Functions.Lifting_Bound(p);
      tstart(timer);
      Pipelined_Polyhedral_Homotopies
        (nt,true,stlb,p,r,mtype,perm,lifsup,mcc,tmv,lq,q,qsols);
      Apply_Induced_Permutation(p,stlb,r,perm,mtype,mcc);
      Set_Continuation_Parameter(qsols,create(0.0));
      Pipelined_Stable_Continuation
        (true,r,mtype,stlb,lifsup,mcc,tmv,lq,mv,smv,qsols0);
      Set_Continuation_Parameter(qsols0,create(0.0));
      put(file,"mixed volume : "); put(file,mv,1); new_line(file);
      put(file,"stable mixed volume : "); put(file,smv,1); new_line(file);
      new_line(file);
      put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
      put_line(file,q);
      new_line(file);
      put_line(file,"START SOLUTIONS :");
      put(file,Length_Of(qsols),natural32(q'last),qsols);
      if not Is_Null(qsols0) then
        new_line(file);
        put_line(file,"START SOLUTIONS with zero coordinates :");
        put(file,Length_Of(qsols0),natural32(q'last),qsols0);
      end if;
      tstop(timer);
      rc := smv;
      rocotime := 0.0;
      hocotime := Elapsed_User_Time(timer);
      new_line(file);
      print_times(file,timer,"pipelined polyhedral continuation");
    end if;
    ended_moment := Ada.Calendar.Clock;
    new_line(file);
    Write_Elapsed_Time(file,start_moment,ended_moment);
    flush(file);
  end Pipelined_Root_Counting;

  procedure Pipelined_Root_Counting
               ( file : in file_type; nt : in integer32;
                 p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out DoblDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    d : natural64;
    mptdeg : Natural_Number;
    mv,smv,tmv : natural32;
    r : integer32;
    mtype,perm : Standard_Integer_Vectors.Link_to_Vector;
    mcc : Mixed_Subdivision;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists;
    timer : Timing_Widget;
    lq : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    zero : constant double_double := create(0.0);

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;
    use Root_Counters_Output;
    use Pipelined_Polyhedral_Drivers;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Pipelined_Count_Rooting 8,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    if nt < 2 or deg then
      Black_Box_Root_Counting
        (file,nt,p,deg,rc,q,qsols,qsols0,rocotime,hocotime);
    else
      Compute_Total_Degree(p,d,mptdeg);
      Write_Total_Degree(file,d,mptdeg);
      flush(file);
      stlb := Floating_Lifting_Functions.Lifting_Bound(p);
      tstart(timer);
      Pipelined_Polyhedral_Homotopies
        (nt,true,stlb,p,r,mtype,perm,lifsup,mcc,tmv,lq,q,qsols);
      Apply_Induced_Permutation(p,stlb,r,perm,mtype,mcc);
      Set_Continuation_Parameter(qsols,create(zero));
      Pipelined_Stable_Continuation
        (true,r,mtype,stlb,lifsup,mcc,tmv,lq,mv,smv,qsols0);
      Set_Continuation_Parameter(qsols0,create(zero));
      put(file,"mixed volume : "); put(file,mv,1); new_line(file);
      put(file,"stable mixed volume : "); put(file,smv,1); new_line(file);
      new_line(file);
      put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
      put_line(file,q);
      new_line(file);
      put_line(file,"START SOLUTIONS :");
      put(file,Length_Of(qsols),natural32(q'last),qsols);
      if not Is_Null(qsols0) then
        new_line(file);
        put_line(file,"START SOLUTIONS with zero coordinates :");
        put(file,Length_Of(qsols0),natural32(q'last),qsols0);
      end if;
      tstop(timer);
      rc := smv;
      rocotime := 0.0;
      hocotime := Elapsed_User_Time(timer);
      new_line(file);
      print_times(file,timer,"pipelined polyhedral continuation");
    end if;
    ended_moment := Ada.Calendar.Clock;
    new_line(file);
    Write_Elapsed_Time(file,start_moment,ended_moment);
    flush(file);
  end Pipelined_Root_Counting;

  procedure Pipelined_Root_Counting
               ( file : in file_type; nt : in integer32;
                 p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out QuadDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    d : natural64;
    mptdeg : Natural_Number;
    mv,smv,tmv : natural32;
    r : integer32;
    mtype,perm : Standard_Integer_Vectors.Link_to_Vector;
    mcc : Mixed_Subdivision;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists;
    timer : Timing_Widget;
    lq : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    zero : constant quad_double := create(0.0);

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;
    use Root_Counters_Output;
    use Pipelined_Polyhedral_Drivers;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Pipelined_Count_Rooting 9,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    if nt < 2 or deg then
      Black_Box_Root_Counting
        (file,nt,p,deg,rc,q,qsols,qsols0,rocotime,hocotime);
    else
      Compute_Total_Degree(p,d,mptdeg);
      Write_Total_Degree(file,d,mptdeg);
      flush(file);
      stlb := Floating_Lifting_Functions.Lifting_Bound(p);
      tstart(timer);
      Pipelined_Polyhedral_Homotopies
        (nt,true,stlb,p,r,mtype,perm,lifsup,mcc,tmv,lq,q,qsols);
      Apply_Induced_Permutation(p,stlb,r,perm,mtype,mcc);
      Set_Continuation_Parameter(qsols,create(zero));
      Pipelined_Stable_Continuation
        (true,r,mtype,stlb,lifsup,mcc,tmv,lq,mv,smv,qsols0);
      Set_Continuation_Parameter(qsols0,create(zero));
      put(file,"mixed volume : "); put(file,mv,1); new_line(file);
      put(file,"stable mixed volume : "); put(file,smv,1); new_line(file);
      new_line(file);
      put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
      put_line(file,q);
      new_line(file);
      put_line(file,"START SOLUTIONS :");
      put(file,Length_Of(qsols),natural32(q'last),qsols);
      if not Is_Null(qsols0) then
        new_line(file);
        put_line(file,"START SOLUTIONS with zero coordinates :");
        put(file,Length_Of(qsols0),natural32(q'last),qsols0);
      end if;
      tstop(timer);
      rc := smv;
      rocotime := 0.0;
      hocotime := Elapsed_User_Time(timer);
      new_line(file);
      print_times(file,timer,"pipelined polyhedral continuation");
    end if;
    ended_moment := Ada.Calendar.Clock;
    new_line(file);
    Write_Elapsed_Time(file,start_moment,ended_moment);
    flush(file);
  end Pipelined_Root_Counting;

-- FOR LAURENT SYSTEMS :

  procedure Black_Box_Root_Counting 
               ( nt : in integer32; silent : in boolean;
                 p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out Standard_Complex_Laur_Systems.Laur_Sys;
                 qsols : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    timer : Timing_Widget;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    lifsup : Link_to_Array_of_Lists;
    mixsub : Mixed_Subdivision;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Black_Box_Root_Counting 10,");
      put_line("for Laurent systems in double precision ...");
    end if;
    tstart(timer);
    Black_Box_Mixed_Volume_Computation
      (p,mix,perm,iprm,lifsup,mixsub,rc,verbose-1);
    tstop(timer);
    if not silent
     then put("mixed volume : "); put(rc,1); new_line;
    end if;
    rocotime := Elapsed_User_Time(timer);
    tstart(timer);
    Black_Box_Polyhedral_Continuation
      (nt,p,mix,lifsup.all,mixsub,q,qsols,verbose-1);
    tstop(timer);
    hocotime := Elapsed_User_Time(timer);
  end Black_Box_Root_Counting;

  procedure Black_Box_Root_Counting 
               ( nt : in integer32; silent : in boolean;
                 p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out DoblDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    timer : Timing_Widget;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    lifsup : Link_to_Array_of_Lists;
    mixsub : Mixed_Subdivision;
    sp : Standard_Complex_Laur_Systems.Laur_Sys(p'range)
       := DoblDobl_Complex_to_Standard_Laur_Sys(p);

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Black_Box_Root_Counting 11,");
      put_line("for Laurent systems in double double precision ...");
    end if;
    tstart(timer);
    Black_Box_Mixed_Volume_Computation
      (sp,mix,perm,iprm,lifsup,mixsub,rc,verbose-1);
    tstop(timer);
    if not silent
     then put("mixed volume : "); put(rc,1); new_line;
    end if;
    rocotime := Elapsed_User_Time(timer);
    tstart(timer);
    Black_Box_Polyhedral_Continuation
      (nt,p,mix,lifsup.all,mixsub,q,qsols,verbose-1);
    tstop(timer);
    Standard_Complex_Laur_Systems.Clear(sp);
    hocotime := Elapsed_User_Time(timer);
  end Black_Box_Root_Counting;

  procedure Black_Box_Root_Counting 
               ( nt : in integer32; silent : in boolean;
                 p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out QuadDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    timer : Timing_Widget;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    lifsup : Link_to_Array_of_Lists;
    mixsub : Mixed_Subdivision;
    sp : Standard_Complex_Laur_Systems.Laur_Sys(p'range)
       := QuadDobl_Complex_to_Standard_Laur_Sys(p);

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Black_Box_Root_Counting 12,");
      put_line("for Laurent systems in quad double precision ...");
    end if;
    tstart(timer);
    Black_Box_Mixed_Volume_Computation
      (sp,mix,perm,iprm,lifsup,mixsub,rc,verbose-1);
    tstop(timer);
    if not silent
     then put("mixed volume : "); put(rc,1); new_line;
    end if;
    rocotime := Elapsed_User_Time(timer);
    tstart(timer);
    Black_Box_Polyhedral_Continuation
      (nt,p,mix,lifsup.all,mixsub,q,qsols,verbose-1);
    tstop(timer);
    Standard_Complex_Laur_Systems.Clear(sp);
    hocotime := Elapsed_User_Time(timer);
  end Black_Box_Root_Counting;

  procedure Black_Box_Root_Counting 
               ( nt : in integer32;
                 p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32; rocos : out Link_to_String;
                 q : out Standard_Complex_Laur_Systems.Laur_Sys;
                 qsols : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    timer : Timing_Widget;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    lifsup : Link_to_Array_of_Lists;
    mixsub : Mixed_Subdivision;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Black_Box_Root_Counting 13,");
      put_line("for Laurent systems in double precision ...");
    end if;
    tstart(timer);
    Black_Box_Mixed_Volume_Computation
      (p,mix,perm,iprm,lifsup,mixsub,rc,verbose-1);
    tstop(timer);
    Append(rocos,"mixed volume : ");
    Append(rocos,Characters_and_Numbers.nConvert(rc));
    rocotime := Elapsed_User_Time(timer);
    tstart(timer);
    Black_Box_Polyhedral_Continuation
      (nt,p,mix,lifsup.all,mixsub,q,qsols,verbose-1);
    tstop(timer);
    hocotime := Elapsed_User_Time(timer);
  end Black_Box_Root_Counting;

  procedure Black_Box_Root_Counting 
               ( nt : in integer32;
                 p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32; rocos : out Link_to_String;
                 q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out DoblDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    timer : Timing_Widget;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    lifsup : Link_to_Array_of_Lists;
    mixsub : Mixed_Subdivision;
    sp : Standard_Complex_Laur_Systems.Laur_Sys(p'range)
       := DoblDobl_Complex_to_Standard_Laur_Sys(p);

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Black_Box_Root_Counting 14,");
      put_line("for Laurent systems in double double precision ...");
    end if;
    tstart(timer);
    Black_Box_Mixed_Volume_Computation
      (sp,mix,perm,iprm,lifsup,mixsub,rc,verbose-1);
    tstop(timer);
    Append(rocos,"mixed volume : ");
    Append(rocos,Characters_and_Numbers.nConvert(rc));
    rocotime := Elapsed_User_Time(timer);
    tstart(timer);
    Black_Box_Polyhedral_Continuation
      (nt,p,mix,lifsup.all,mixsub,q,qsols,verbose-1);
    tstop(timer);
    Standard_Complex_Laur_Systems.Clear(sp);
    hocotime := Elapsed_User_Time(timer);
  end Black_Box_Root_Counting;

  procedure Black_Box_Root_Counting 
               ( nt : in integer32;
                 p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32; rocos : out Link_to_String;
                 q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out QuadDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    timer : Timing_Widget;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    lifsup : Link_to_Array_of_Lists;
    mixsub : Mixed_Subdivision;
    sp : Standard_Complex_Laur_Systems.Laur_Sys(p'range)
       := QuadDobl_Complex_to_Standard_Laur_Sys(p);

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Black_Box_Root_Counting 15,");
      put_line("for Laurent systems in quad double precision ...");
    end if;
    tstart(timer);
    Black_Box_Mixed_Volume_Computation
      (sp,mix,perm,iprm,lifsup,mixsub,rc,verbose-1);
    tstop(timer);
    Append(rocos,"mixed volume : ");
    Append(rocos,Characters_and_Numbers.nConvert(rc));
    rocotime := Elapsed_User_Time(timer);
    tstart(timer);
    Black_Box_Polyhedral_Continuation
      (nt,p,mix,lifsup.all,mixsub,q,qsols,verbose-1);
    tstop(timer);
    Standard_Complex_Laur_Systems.Clear(sp);
    hocotime := Elapsed_User_Time(timer);
  end Black_Box_Root_Counting;

  procedure Black_Box_Root_Counting 
               ( file : in file_type; nt : in integer32;
                 p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out Standard_Complex_Laur_Systems.Laur_Sys;
                 qsols : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    lifsup : Link_to_Array_of_Lists;
    mixsub : Mixed_Subdivision;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Black_Box_Root_Counting 16,");
      put_line("for Laurent systems in double precision ...");
    end if;
    tstart(timer);
    Black_Box_Mixed_Volume_Computation
      (p,mix,perm,iprm,lifsup,mixsub,rc,verbose-1);
    tstop(timer);
    new_line(file);
    put(file,"mixed volume : "); put(file,rc,1); new_line(file);
    new_line(file);
    print_times(file,timer,"mixed-volume computation");
    flush(file);
    rocotime := Elapsed_User_Time(timer);
    put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
    tstart(timer);
    Black_Box_Polyhedral_Continuation
      (nt,p,mix,lifsup.all,mixsub,q,qsols,verbose-1);
    tstop(timer);
    hocotime := Elapsed_User_Time(timer);
    new_line(file);
    put_line(file,q);
    new_line(file);
    put_line(file,"START SOLUTIONS : ");
    new_line(file);
    if not Is_Null(qsols)
     then put(file,Length_Of(qsols),natural32(q'last),qsols);
    end if;
    new_line(file);
    print_times(file,timer,"polyhedral homotopy continuation");
    flush(file);
  end Black_Box_Root_Counting;

  procedure Black_Box_Root_Counting 
               ( file : in file_type; nt : in integer32;
                 p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out DoblDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;

    timer : Timing_Widget;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    lifsup : Link_to_Array_of_Lists;
    mixsub : Mixed_Subdivision;
    sp : Standard_Complex_Laur_Systems.Laur_Sys(p'range)
       := DoblDobl_Complex_to_Standard_Laur_Sys(p);

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Black_Box_Root_Counting 17,");
      put_line("for Laurent systems in double double precision ...");
    end if;
    tstart(timer);
    Black_Box_Mixed_Volume_Computation
      (sp,mix,perm,iprm,lifsup,mixsub,rc,verbose-1);
    tstop(timer);
    new_line(file);
    put(file,"mixed volume : "); put(file,rc,1); new_line(file);
    new_line(file);
    print_times(file,timer,"mixed-volume computation");
    flush(file);
    rocotime := Elapsed_User_Time(timer);
    put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
   -- careful to apply the permutation on p!
    tstart(timer);
    Black_Box_Polyhedral_Continuation
      (nt,p,mix,lifsup.all,mixsub,q,qsols,verbose-1);
    tstop(timer);
    hocotime := Elapsed_User_Time(timer);
    new_line(file);
    put_line(file,q);
    new_line(file);
    put_line(file,"START SOLUTIONS : ");
    new_line(file);
    if not Is_Null(qsols)
     then put(file,Length_Of(qsols),natural32(q'last),qsols);
    end if;
    new_line(file);
    print_times(file,timer,"polyhedral homotopy continuation");
    flush(file);
    Standard_Complex_Laur_Systems.Clear(sp);
  end Black_Box_Root_Counting;

  procedure Black_Box_Root_Counting 
               ( file : in file_type; nt : in integer32;
                 p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out QuadDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;

    timer : Timing_Widget;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    lifsup : Link_to_Array_of_Lists;
    mixsub : Mixed_Subdivision;
    sp : Standard_Complex_Laur_Systems.Laur_Sys(p'range)
       := QuadDobl_Complex_to_Standard_Laur_Sys(p);

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Black_Box_Root_Counting 18,");
      put_line("for Laurent systems in quad double precision ...");
    end if;
    tstart(timer);
    Black_Box_Mixed_Volume_Computation
      (sp,mix,perm,iprm,lifsup,mixsub,rc,verbose-1);
    tstop(timer);
    new_line(file);
    put(file,"mixed volume : "); put(file,rc,1); new_line(file);
    new_line(file);
    print_times(file,timer,"mixed-volume computation");
    flush(file);
    rocotime := Elapsed_User_Time(timer);
    put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
    tstart(timer);
    Black_Box_Polyhedral_Continuation
      (nt,p,mix,lifsup.all,mixsub,q,qsols,verbose-1);
    tstop(timer);
    hocotime := Elapsed_User_Time(timer);
    new_line(file);
    put_line(file,q);
    new_line(file);
    put_line(file,"START SOLUTIONS : ");
    new_line(file);
    if not Is_Null(qsols)
     then put(file,Length_Of(qsols),natural32(q'last),qsols);
    end if;
    new_line(file);
    print_times(file,timer,"polyhedral homotopy continuation");
    flush(file);
    Standard_Complex_Laur_Systems.Clear(sp);
  end Black_Box_Root_Counting;

-- PIPELINING FOR LAURENT SYSTEMS :

  procedure Pipelined_Root_Counting 
               ( nt : in integer32; silent : in boolean;
                 p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out Standard_Complex_Laur_Systems.Laur_Sys;
                 qsols : out Standard_Complex_Solutions.Solution_List;
                 elaptime : out duration; verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
    use Pipelined_Polyhedral_Drivers;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    timer : Timing_Widget;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Pipelined_Count_Rooting 10,");
      put_line("for Laurent systems in standard double precision ...");
    end if;
    tstart(timer);
    Pipelined_Polyhedral_Homotopies(nt,p,rc,q,qsols);
    tstop(timer);
    if not silent
     then put("mixed volume : "); put(rc,1); new_line;
    end if;
    elaptime := Elapsed_User_Time(timer);
    ended_moment := Ada.Calendar.Clock;
    if not silent then
      new_line;
      print_times(standard_output,timer,"pipelined polyhedral homotopies");
      new_line;
      Write_Elapsed_Time(standard_output,start_moment,ended_moment);
    end if;
  end Pipelined_Root_Counting;

  procedure Pipelined_Root_Counting 
               ( nt : in integer32; silent : in boolean;
                 p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out DoblDobl_Complex_Solutions.Solution_List;
                 elaptime : out duration; verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;
    use Pipelined_Polyhedral_Drivers;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    timer : Timing_Widget;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Pipelined_Count_Rooting 11,");
      put_line("for Laurent systems in double double precision ...");
    end if;
    tstart(timer);
    Pipelined_Polyhedral_Homotopies(nt,p,rc,q,qsols);
    tstop(timer);
    if not silent
     then put("mixed volume : "); put(rc,1); new_line;
    end if;
    elaptime := Elapsed_User_Time(timer);
    ended_moment := Ada.Calendar.Clock;
    if not silent then
      new_line;
      print_times(standard_output,timer,"pipelined polyhedral homotopies");
      new_line;
      Write_Elapsed_Time(standard_output,start_moment,ended_moment);
    end if;
  end Pipelined_Root_Counting;

  procedure Pipelined_Root_Counting 
               ( nt : in integer32; silent : in boolean;
                 p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out QuadDobl_Complex_Solutions.Solution_List;
                 elaptime : out duration; verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;
    use Pipelined_Polyhedral_Drivers;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    timer : Timing_Widget;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Pipelined_Count_Rooting 12,");
      put_line("for Laurent systems in quad double precision ...");
    end if;
    tstart(timer);
    Pipelined_Polyhedral_Homotopies(nt,p,rc,q,qsols);
    tstop(timer);
    if not silent
     then put("mixed volume : "); put(rc,1); new_line;
    end if;
    elaptime := Elapsed_User_Time(timer);
    ended_moment := Ada.Calendar.Clock;
    if not silent then
      new_line;
      print_times(standard_output,timer,"pipelined polyhedral homotopies");
      new_line;
      Write_Elapsed_Time(standard_output,start_moment,ended_moment);
    end if;
  end Pipelined_Root_Counting;

  procedure Pipelined_Root_Counting 
               ( nt : in integer32;
                 p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32; rocos : out Link_to_String;
                 q : out Standard_Complex_Laur_Systems.Laur_Sys;
                 qsols : out Standard_Complex_Solutions.Solution_List;
                 elaptime : out duration; verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
    use Pipelined_Polyhedral_Drivers;

    timer : Timing_Widget;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Pipelined_Count_Rooting 13,");
      put_line("for Laurent systems in standard double precision ...");
    end if;
    tstart(timer);
    Pipelined_Polyhedral_Homotopies(nt,p,rc,q,qsols);
    tstop(timer);
    Append(rocos,"mixed volume : ");
    Append(rocos,Characters_and_Numbers.nConvert(rc));
    elaptime := Elapsed_User_Time(timer);
  end Pipelined_Root_Counting;

  procedure Pipelined_Root_Counting 
               ( nt : in integer32;
                 p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32; rocos : out Link_to_String;
                 q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out DoblDobl_Complex_Solutions.Solution_List;
                 elaptime : out duration; verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;
    use Pipelined_Polyhedral_Drivers;

    timer : Timing_Widget;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Pipelined_Count_Rooting 14,");
      put_line("for Laurent systems in double double precision ...");
    end if;
    tstart(timer);
    Pipelined_Polyhedral_Homotopies(nt,p,rc,q,qsols);
    tstop(timer);
    Append(rocos,"mixed volume : ");
    Append(rocos,Characters_and_Numbers.nConvert(rc));
    elaptime := Elapsed_User_Time(timer);
  end Pipelined_Root_Counting;

  procedure Pipelined_Root_Counting 
               ( nt : in integer32;
                 p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32; rocos : out Link_to_String;
                 q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out QuadDobl_Complex_Solutions.Solution_List;
                 elaptime : out duration; verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;
    use Pipelined_Polyhedral_Drivers;

    timer : Timing_Widget;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Pipelined_Count_Rooting 15,");
      put_line("for Laurent systems in quad double precision ...");
    end if;
    tstart(timer);
    Pipelined_Polyhedral_Homotopies(nt,p,rc,q,qsols);
    tstop(timer);
    Append(rocos,"mixed volume : ");
    Append(rocos,Characters_and_Numbers.nConvert(rc));
    elaptime := Elapsed_User_Time(timer);
  end Pipelined_Root_Counting;

  procedure Pipelined_Root_Counting 
               ( file : in file_type; nt : in integer32;
                 p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out Standard_Complex_Laur_Systems.Laur_Sys;
                 qsols : out Standard_Complex_Solutions.Solution_List;
                 elaptime : out duration; verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
    use Pipelined_Polyhedral_Drivers;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    timer : Timing_Widget;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Pipelined_Count_Rooting 16,");
      put_line("for Laurent systems in standard double precision ...");
    end if;
    tstart(timer);
    Pipelined_Polyhedral_Homotopies(nt,p,rc,q,qsols);
    tstop(timer);
    new_line(file);
    put(file,"mixed volume : "); put(file,rc,1); new_line(file);
    new_line(file);
    put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
    put_line(file,q);
    new_line(file);
    put_line(file,"START SOLUTIONS : ");
    new_line(file);
    if not Is_Null(qsols)
     then put(file,Length_Of(qsols),natural32(q'last),qsols);
    end if;
    new_line(file);
    print_times(file,timer,"pipelined polyhedral homotopy continuation");
    elaptime := Elapsed_User_Time(timer);
    ended_moment := Ada.Calendar.Clock;
    new_line(file);
    Write_Elapsed_Time(file,start_moment,ended_moment);
    flush(file);
  end Pipelined_Root_Counting;

  procedure Pipelined_Root_Counting 
               ( file : in file_type; nt : in integer32;
                 p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out DoblDobl_Complex_Solutions.Solution_List;
                 elaptime : out duration; verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;
    use Pipelined_Polyhedral_Drivers;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    timer : Timing_Widget;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Pipelined_Count_Rooting 17,");
      put_line("for Laurent systems in double double precision ...");
    end if;
    tstart(timer);
    Pipelined_Polyhedral_Homotopies(nt,p,rc,q,qsols);
    tstop(timer);
    new_line(file);
    put(file,"mixed volume : "); put(file,rc,1); new_line(file);
    new_line(file);
    put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
    put_line(file,q);
    new_line(file);
    put_line(file,"START SOLUTIONS : ");
    new_line(file);
    if not Is_Null(qsols)
     then put(file,Length_Of(qsols),natural32(q'last),qsols);
    end if;
    new_line(file);
    print_times(file,timer,"pipelined polyhedral homotopy continuation");
    elaptime := Elapsed_User_Time(timer);
    ended_moment := Ada.Calendar.Clock;
    new_line(file);
    Write_Elapsed_Time(file,start_moment,ended_moment);
    flush(file);
  end Pipelined_Root_Counting;

  procedure Pipelined_Root_Counting 
               ( file : in file_type; nt : in integer32;
                 p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out QuadDobl_Complex_Solutions.Solution_List;
                 elaptime : out duration; verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;
    use Pipelined_Polyhedral_Drivers;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    timer : Timing_Widget;

  begin
    if verbose > 0 then
      put_line("-> in black_box_root_counters.Pipelined_Count_Rooting 18,");
      put_line("for Laurent systems in quad double precision ...");
    end if;
    tstart(timer);
    Pipelined_Polyhedral_Homotopies(nt,p,rc,q,qsols);
    tstop(timer);
    new_line(file);
    put(file,"mixed volume : "); put(file,rc,1); new_line(file);
    new_line(file);
    put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
    put_line(file,q);
    new_line(file);
    put_line(file,"START SOLUTIONS : ");
    new_line(file);
    if not Is_Null(qsols)
     then put(file,Length_Of(qsols),natural32(q'last),qsols);
    end if;
    new_line(file);
    print_times(file,timer,"pipelined polyhedral homotopy continuation");
    elaptime := Elapsed_User_Time(timer);
    ended_moment := Ada.Calendar.Clock;
    new_line(file);
    Write_Elapsed_Time(file,start_moment,ended_moment);
    flush(file);
  end Pipelined_Root_Counting;

  procedure Main ( nt : in natural32; infilename,outfilename : in string;
                   verbose : in integer32 := 0 ) is

    use Standard_Complex_Poly_Systems;

    procedure Read_System
                ( file : in out file_type; filename : in string;
                  lp : out Link_to_Poly_Sys ) is
    begin
      if filename /= "" then
        Open_Input_File(file,filename);
        get(file,lp);
      end if;
    exception
      when others => put_line("Something is wrong with argument file...");
                     lp := null; return;
    end Read_System;

    lp,lq : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    infile,outfile : file_type;
    rc : natural32;
    roco,poco : duration;
    qsols,qsols0 : Standard_Complex_Solutions.Solution_List;

  begin
    if verbose > 0 then
      put("At verbose level "); put(verbose,1);
      put_line(", in black_box_root_counters.Main ...");
    end if;
    Read_System(infile,infilename,lp);
    if lp = null then
      new_line;
      get(lp);
    end if;
    Create_Output_File(outfile,outfilename);
    put(outfile,lp.all);
    lq := new Poly_Sys(lp'range);
    Black_Box_Root_Counting
      (outfile,integer32(nt),lp.all,false,rc,lq.all,qsols,qsols0,roco,poco);
  end Main;

end Black_Box_Root_Counters;
