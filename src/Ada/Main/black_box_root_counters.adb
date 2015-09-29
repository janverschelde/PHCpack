with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Multprec_Natural_Numbers_io;        use Multprec_Natural_Numbers_io;
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
with Partitions_of_Sets_Of_Unknowns_io;  use Partitions_of_Sets_of_Unknowns_io;
with m_Homogeneous_Bezout_Numbers;       use m_Homogeneous_Bezout_Numbers;
with m_Homogeneous_Start_Systems;        use m_Homogeneous_Start_Systems;
with Set_Structure,Set_Structure_io;
with Standard_Linear_Product_System;
with Standard_Complex_Prod_Systems;      use Standard_Complex_Prod_Systems;
with Standard_Complex_Prod_Systems_io;   use Standard_Complex_Prod_Systems_io;
with Standard_Complex_Prod_Planes;
with Random_Product_Start_Systems;       use Random_Product_Start_Systems;
with Floating_Mixed_Subdivisions_io;
with Induced_Permutations;
with Black_Mixed_Volume_Computations;    use Black_Mixed_Volume_Computations;
with Black_Polyhedral_Continuations;     use Black_Polyhedral_Continuations;

--with Standard_Integer_Vectors_io;
-- use Standard_Integer_Vectors_io;

package body Black_Box_Root_Counters is

  chicken_mv : constant natural32 := 40;  -- arbitrary bound ...

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
    res : natural64 := Set_Structure_Bound(q);

  begin
    Standard_Complex_Poly_Systems.Clear(q);
    return res;
  end Set_Structure_Bound;

  function Set_Structure_Bound
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys ) return natural64 is

    q : Standard_Complex_Poly_Systems.Poly_Sys(p'range)
      := QuadDobl_Complex_to_Standard_Poly_Sys(p);
    res : natural64 := Set_Structure_Bound(q);

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
                 rocotime : out duration ) is

    timer : timing_widget;
    n : constant natural32 := natural32(p'length);
    d,bz,bs : natural64 := 0;
    m,mv,smv,tmv : natural32 := 0;
    z : Partition(1..n);
    no_mv : constant boolean := (n > chicken_mv);
    mixsub : Mixed_Subdivision;
    orgcnt,stbcnt : natural32;

  begin
    tstart(timer);
    declare
    begin
      d := Total_Degree(p);
    exception
      when others =>
        mptode := Total_Degree(p); d := 0;
        put("Huge total degree : "); put(mptode); new_line;
    end;
    if d = 0 and Equal(mptode,0)       -- patch for GNAT optimizers ...
     then mptode := Total_Degree(p);
    end if;
    declare
    begin
      PB(p,bz,m,z);
      bs := Set_Structure_Bound(p);
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
           mv,smv,tmv,orgcnt,stbcnt);
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
                 rocotime : out duration ) is

    q : Standard_Complex_Poly_Systems.Poly_Sys(p'range)
      := DoblDobl_Complex_to_Standard_Poly_Sys(p);

  begin
    Count_Roots(q,deg,tode,mptode,mhbz,setb,mivo,stmv,zz,nz,stlb,lifsup,
                mix,perm,iprm,orgmcc,stbmcc,rocotime);
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
                 rocotime : out duration ) is

    q : Standard_Complex_Poly_Systems.Poly_Sys(p'range)
      := QuadDobl_Complex_to_Standard_Poly_Sys(p);

  begin
    Count_Roots(q,deg,tode,mptode,mhbz,setb,mivo,stmv,zz,nz,stlb,lifsup,
                mix,perm,iprm,orgmcc,stbmcc,rocotime);
    if iprm /= null
     then Induced_Permutations.Permute(iprm.all,p);
    end if;
    Standard_Complex_Poly_Systems.Clear(q);
  end Count_Roots;

  procedure Write_Root_Counts
               ( file : in file_type; no_mv : in boolean;
                 d : in natural64; mp_d : in Natural_Number;
                 m : in natural32; bz,bs : in natural64;
                 mv,smv : in natural32; z : in Partition ) is
  begin
    new_line(file);
    put_line(file,"ROOT COUNTS :");
    new_line(file);
    put(file,"total degree : ");
    if d = 0
     then put(file,mp_d,1); new_line(file);
     else put(file,d,1); new_line(file);
    end if;
    if m > 1 then
      put(file,m,1); put(file,"-homogeneous Bezout number : ");
      put(file,bz,1); new_line(file);
      put(file,"  with partition : "); put(file,z(1..m)); new_line(file);
    end if;
    put(file,"general linear-product Bezout number : ");
    put(file,bs,1); new_line(file);
    if bs > 0 then
      put_line(file,"  based on the set structure :");
      Set_Structure_io.put(file);
    end if;
    if not no_mv then
      put(file,"mixed volume : "); put(file,mv,1); new_line(file);
      put(file,"stable mixed volume : "); put(file,smv,1); new_line(file);
    end if;
  end Write_Root_Counts;

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
                 rocotime : out duration ) is

    timer : timing_widget;
    n : constant natural32 := natural32(p'length);
    d,bz,bs : natural64 := 0;
    m,mv,smv,tmv,orgcnt,stbcnt : natural32 := 0;
    z : Partition(1..n);
    no_mv : constant boolean := deg or (n > chicken_mv);
    mixsub : Mixed_Subdivision;

  begin
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
   -- put("total degree : "); put(d,1); new_line;
    declare
    begin
      PB(p,bz,m,z); 
     -- put(m,1); put("-homogeneous Bezout number : "); put(bz,1); new_line;
      bs := Set_Structure_Bound(p);
     -- put("set structure Bezout bound : "); put(bs,1); new_line;
    exception
      when others => bz := 0; m := 1; bs := 0;
    end;
    if m = 1
     then bz := d;
    end if;
    if not no_mv then
      declare -- problems with systems with one monomial equation
      begin
        Black_Box_Mixed_Volume_Computation
          (p,mix,perm,iprm,stlb,lifsup,mixsub,orgmcc,stbmcc,
           mv,smv,tmv,orgcnt,stbcnt);
       -- put("the mixed volume : "); put(mv,1); new_line;
       -- put("the stable mixed volume : "); put(smv,1); new_line;
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
                 rocotime : out duration ) is

    sp : Standard_Complex_Poly_Systems.Poly_Sys(p'range)
       := DoblDobl_Complex_to_Standard_Poly_Sys(p);

  begin
    Count_Roots(file,sp,deg,tode,mptode,mhbz,setb,mivo,stmv,
                zz,nz,stlb,lifsup,mix,perm,iprm,orgmcc,stbmcc,rocotime);
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
                 rocotime : out duration ) is

    sp : Standard_Complex_Poly_Systems.Poly_Sys(p'range)
       := QuadDobl_Complex_to_Standard_Poly_Sys(p);

  begin
    Count_Roots(file,sp,deg,tode,mptode,mhbz,setb,mivo,stmv,
                zz,nz,stlb,lifsup,mix,perm,iprm,orgmcc,stbmcc,rocotime);
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
                hocotime : out duration ) is

    timer : timing_widget;
    n : constant natural32 := natural32(p'length);
    nl : natural32;
    gap : constant natural32 := smv - mv;

  begin
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
          (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0);
      end if;
    else 
      Black_Box_Polyhedral_Continuation
        (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0);
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
                hocotime : out duration ) is

    timer : timing_widget;
    n : constant natural32 := natural32(p'length);
    nl : natural32;
    gap : constant natural32 := smv - mv;

  begin
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
          (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0);
      end if;
    else 
      Black_Box_Polyhedral_Continuation
        (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0);
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
                hocotime : out duration ) is

    timer : timing_widget;
    n : constant natural32 := natural32(p'length);
    nl : natural32;
    gap : constant natural32 := smv - mv;

  begin
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
          (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0);
      end if;
    else 
      Black_Box_Polyhedral_Continuation
        (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0);
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
                hocotime : out duration ) is

    use Standard_Complex_Solutions;

    timer : timing_widget;
    n : constant natural32 := natural32(p'length);
    rq : Prod_Sys(p'range);
    nl : natural32;
    prodform : boolean := false;
    wsols : Solution_List;
    gap : constant natural32 := smv - mv;

  begin
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
          (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0);
      end if;
    else
      put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
      Black_Box_Polyhedral_Continuation
        (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0);
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
                hocotime : out duration ) is

    use DoblDobl_Complex_Solutions;

    timer : timing_widget;
    n : constant natural32 := natural32(p'length);
    rq : Prod_Sys(p'range);
    nl : natural32;
    prodform : boolean := false;
    wsols : Solution_List;
    gap : constant natural32 := smv - mv;

  begin
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
          (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0);
      end if;
    else
      put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
      Black_Box_Polyhedral_Continuation
        (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0);
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
                hocotime : out duration ) is

    use QuadDobl_Complex_Solutions;

    timer : timing_widget;
    n : constant natural32 := natural32(p'length);
    rq : Prod_Sys(p'range);
    nl : natural32;
    prodform : boolean := false;
    wsols : Solution_List;
    gap : constant natural32 := smv - mv;

  begin
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
          (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0);
      end if;
    else
      put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
      Black_Box_Polyhedral_Continuation
        (nt,p,mix,stlb,lifted.all,orgmcc,stbmcc,q,qsols,qsols0);
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
                 rocotime,hocotime : out duration ) is

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

  begin
    Count_Roots(p,deg,d,mptdeg,bz,bs,mv,smv,z,nz,
                stlb,lifsup,mix,perm,iprm,orgmcc,stbmcc,rocotime);
    if not silent
     then Write_Root_Counts(standard_output,no_mv,d,mptdeg,nz,bz,bs,mv,smv,z);
    end if;
    Construct_Start_System
      (nt,p,d,bz,bs,mv,smv,z(1..nz),mix,stlb,lifsup,orgmcc,stbmcc,wrc,
       q,qsols,qsols0,hocotime);
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
                 rocotime,hocotime : out duration ) is

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

  begin
    Count_Roots(p,deg,d,mptdeg,bz,bs,mv,smv,z,nz,
                stlb,lifsup,mix,perm,iprm,orgmcc,stbmcc,rocotime);
    if not silent
     then Write_Root_Counts(standard_output,no_mv,d,mptdeg,nz,bz,bs,mv,smv,z);
    end if;
    Construct_Start_System
      (nt,p,d,bz,bs,mv,smv,z(1..nz),mix,stlb,lifsup,orgmcc,stbmcc,wrc,
       q,qsols,qsols0,hocotime);
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
                 rocotime,hocotime : out duration ) is

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

  begin
    Count_Roots(p,deg,d,mptdeg,bz,bs,mv,smv,z,nz,
                stlb,lifsup,mix,perm,iprm,orgmcc,stbmcc,rocotime);
    if not silent
     then Write_Root_Counts(standard_output,no_mv,d,mptdeg,nz,bz,bs,mv,smv,z);
    end if;
    Construct_Start_System
      (nt,p,d,bz,bs,mv,smv,z(1..nz),mix,stlb,lifsup,orgmcc,stbmcc,wrc,
       q,qsols,qsols0,hocotime);
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
                 rocotime,hocotime : out duration ) is

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
    Count_Roots(file,p,deg,d,mptdeg,bz,bs,mv,smv,z,nz,
                stlb,lifsup,mix,perm,iprm,orgmcc,stbmcc,rocotime);
    Construct_Start_System
      (file,nt,p,d,bz,bs,mv,smv,z(1..nz),mix,stlb,lifsup,orgmcc,stbmcc,
       wrc,q,qsols,qsols0,hocotime);
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
                 rocotime,hocotime : out duration ) is

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
    Count_Roots(file,p,deg,d,mptdeg,bz,bs,mv,smv,z,nz,
                stlb,lifsup,mix,perm,iprm,orgmcc,stbmcc,rocotime);
    Construct_Start_System
      (file,nt,p,d,bz,bs,mv,smv,z(1..nz),mix,stlb,lifsup,orgmcc,stbmcc,
       wrc,q,qsols,qsols0,hocotime);
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
                 rocotime,hocotime : out duration ) is

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
    Count_Roots(file,p,deg,d,mptdeg,bz,bs,mv,smv,z,nz,
                stlb,lifsup,mix,perm,iprm,orgmcc,stbmcc,rocotime);
    Construct_Start_System
      (file,nt,p,d,bz,bs,mv,smv,z(1..nz),mix,stlb,lifsup,orgmcc,stbmcc,
       wrc,q,qsols,qsols0,hocotime);
    rc := natural32(wrc);
    Clear(z);
    Clear(mix);
    Deep_Clear(lifsup);
    Deep_Clear(orgmcc);
    Deep_Clear(stbmcc);
  end Black_Box_Root_Counting;

  procedure Black_Box_Root_Counting 
               ( nt : in integer32; silent : in boolean;
                 p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out Standard_Complex_Laur_Systems.Laur_Sys;
                 qsols : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration ) is

    timer : Timing_Widget;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    lifsup : Link_to_Array_of_Lists;
    mixsub : Mixed_Subdivision;

  begin
    tstart(timer);
    Black_Box_Mixed_Volume_Computation(p,mix,perm,iprm,lifsup,mixsub,rc);
    tstop(timer);
    if not silent
     then put("mixed volume : "); put(rc,1); new_line;
    end if;
    rocotime := Elapsed_User_Time(timer);
    tstart(timer);
    Black_Box_Polyhedral_Continuation(nt,p,mix,lifsup.all,mixsub,q,qsols);
    tstop(timer);
    hocotime := Elapsed_User_Time(timer);
  end Black_Box_Root_Counting;

  procedure Black_Box_Root_Counting 
               ( nt : in integer32; silent : in boolean;
                 p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out DoblDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration ) is

    timer : Timing_Widget;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    lifsup : Link_to_Array_of_Lists;
    mixsub : Mixed_Subdivision;
    sp : Standard_Complex_Laur_Systems.Laur_Sys(p'range)
       := DoblDobl_Complex_to_Standard_Laur_Sys(p);

  begin
    tstart(timer);
    Black_Box_Mixed_Volume_Computation(sp,mix,perm,iprm,lifsup,mixsub,rc);
    tstop(timer);
    if not silent
     then put("mixed volume : "); put(rc,1); new_line;
    end if;
    rocotime := Elapsed_User_Time(timer);
    tstart(timer);
    Black_Box_Polyhedral_Continuation(nt,p,mix,lifsup.all,mixsub,q,qsols);
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
                 rocotime,hocotime : out duration ) is

    timer : Timing_Widget;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    lifsup : Link_to_Array_of_Lists;
    mixsub : Mixed_Subdivision;
    sp : Standard_Complex_Laur_Systems.Laur_Sys(p'range)
       := QuadDobl_Complex_to_Standard_Laur_Sys(p);

  begin
    tstart(timer);
    Black_Box_Mixed_Volume_Computation(sp,mix,perm,iprm,lifsup,mixsub,rc);
    tstop(timer);
    if not silent
     then put("mixed volume : "); put(rc,1); new_line;
    end if;
    rocotime := Elapsed_User_Time(timer);
    tstart(timer);
    Black_Box_Polyhedral_Continuation(nt,p,mix,lifsup.all,mixsub,q,qsols);
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
                 rocotime,hocotime : out duration ) is

    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    lifsup : Link_to_Array_of_Lists;
    mixsub : Mixed_Subdivision;

  begin
    tstart(timer);
    Black_Box_Mixed_Volume_Computation(p,mix,perm,iprm,lifsup,mixsub,rc);
    tstop(timer);
    new_line(file);
    put(file,"mixed volume : "); put(file,rc,1); new_line(file);
    new_line(file);
    print_times(file,timer,"mixed-volume computation");
    flush(file);
    rocotime := Elapsed_User_Time(timer);
    put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
    tstart(timer);
    Black_Box_Polyhedral_Continuation(nt,p,mix,lifsup.all,mixsub,q,qsols);
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
                 rocotime,hocotime : out duration ) is

    use DoblDobl_Complex_Solutions;

    timer : Timing_Widget;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    lifsup : Link_to_Array_of_Lists;
    mixsub : Mixed_Subdivision;
    sp : Standard_Complex_Laur_Systems.Laur_Sys(p'range)
       := DoblDobl_Complex_to_Standard_Laur_Sys(p);

  begin
    tstart(timer);
    Black_Box_Mixed_Volume_Computation(sp,mix,perm,iprm,lifsup,mixsub,rc);
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
    Black_Box_Polyhedral_Continuation(nt,p,mix,lifsup.all,mixsub,q,qsols);
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
                 rocotime,hocotime : out duration ) is

    use QuadDobl_Complex_Solutions;

    timer : Timing_Widget;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    lifsup : Link_to_Array_of_Lists;
    mixsub : Mixed_Subdivision;
    sp : Standard_Complex_Laur_Systems.Laur_Sys(p'range)
       := QuadDobl_Complex_to_Standard_Laur_Sys(p);

  begin
    tstart(timer);
    Black_Box_Mixed_Volume_Computation(sp,mix,perm,iprm,lifsup,mixsub,rc);
    tstop(timer);
    new_line(file);
    put(file,"mixed volume : "); put(file,rc,1); new_line(file);
    new_line(file);
    print_times(file,timer,"mixed-volume computation");
    flush(file);
    rocotime := Elapsed_User_Time(timer);
    put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
    tstart(timer);
    Black_Box_Polyhedral_Continuation(nt,p,mix,lifsup.all,mixsub,q,qsols);
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

end Black_Box_Root_Counters;
