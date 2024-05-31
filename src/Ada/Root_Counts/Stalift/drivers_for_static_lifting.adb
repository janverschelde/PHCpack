with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_VecVecs;
with Standard_Complex_VecVecs;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Functions;    use Standard_Complex_Laur_Functions;
with Standard_Complex_Laur_SysFun;       use Standard_Complex_Laur_SysFun;
with Standard_Complex_Laur_Jacomats;     use Standard_Complex_Laur_Jacomats;
with Exponent_Vectors;                   use Exponent_Vectors;
with Standard_Poly_Laur_Convertors;      use Standard_Poly_Laur_Convertors;
with Standard_Laur_Poly_Convertors;      use Standard_Laur_Poly_Convertors;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Lists_of_Integer_Vectors_io;        use Lists_of_Integer_Vectors_io;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Floating_Integer_Convertors;
with Integer_Faces_of_Polytope;
with Floating_Faces_of_Polytope;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Integer_Mixed_Subdivisions_io;      use Integer_Mixed_Subdivisions_io;
with Floating_Mixed_Subdivisions_io;     use Floating_Mixed_Subdivisions_io;
with Mixed_Volume_Computation;           use Mixed_Volume_Computation;
with Integer_Lifting_Utilities;          use Integer_Lifting_Utilities;
with Floating_Lifting_Functions;         use Floating_Lifting_Functions;
with Floating_Lifting_Utilities;         use Floating_Lifting_Utilities;
with Integer_Pruning_Methods;            use Integer_Pruning_Methods;
with Floating_Pruning_Methods;           use Floating_Pruning_Methods;
with Driver_for_Criterion;
with Main_Lifting_Functions;
with Pruning_Statistics;
with Integer_Polyhedral_Continuation;    use Integer_Polyhedral_Continuation;
with Floating_Polyhedral_Continuation;   use Floating_Polyhedral_Continuation;
with Stable_Polyhedral_Continuation;     use Stable_Polyhedral_Continuation;
with Drivers_for_Coefficient_Systems;    use Drivers_for_Coefficient_Systems;
with Multitasking_Polyhedral_Trackers;   use Multitasking_Polyhedral_Trackers;

package body Drivers_for_Static_Lifting is

  procedure Static_Lifting_Info is

    i : array(1..11) of string(1..65);

  begin
    i( 1):="  Static lifting is a  general  procedure  to  construct  regular";
    i( 2):="mixed  subdivisions  of  a tuple of polytopes.   For mixed-volume";
    i( 3):="computation, only those cells that are  spanned  by  a  tuple  of";
    i( 4):="edges  contribute  to  the mixed volume.  These cells are the so-";
    i( 5):="called mixed cells in the subdivision.  The collection  of  mixed";
    i( 6):="cells  is  computed  efficiently by pruning in the tree of lifted";
    i( 7):="edge-edge combinations.                                          ";
    i( 8):="  These mixed cells provide the start systems in  the  polyhedral";
    i( 9):="homotopy methods used to solve a random coefficient start system.";
    i(10):="Recursion is applied in case the lifting does not induce at  once";
    i(11):="a fine mixed subdivision.                                        ";
    for k in i'range loop
      put_line(i(k));
    end loop;
  end Static_Lifting_Info;

  procedure Compute_Mixture 
              ( file : in file_type; n : in integer32; compmix : in boolean;
                sup : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                mix : in out Standard_Integer_Vectors.Link_to_Vector;
                permsys : in out Poly_Sys;
                vrblvl : in integer32 := 0 ) is

    perm : Standard_Integer_Vectors.Link_to_Vector;

  begin
    if vrblvl > 0
     then put_line("-> in drivers_for_static_lifting.Compute_Mixture 1 ...");
    end if;
    if compmix
     then Compute_Mixture(sup,mix,perm);
     else perm := Compute_Permutation(n,mix.all,sup);
    end if;
    permsys := Permute(permsys,perm);
    new_line(file);
    put(file,"TYPE OF MIXTURE : "); put(file,"#supports : ");
    put(file,mix'last,1);
    put(file,"  occurrences : "); put(file,mix);
    if mix'last /= n then
      new_line(file);
      put(file,"  the permutation : "); put(file,perm);
      new_line(file);
    end if;
    new_line(file);
    Clear(perm);
  end Compute_Mixture;

  procedure Compute_Mixture 
              ( file : in file_type; n : in integer32; compmix : in boolean;
                sup : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                mix : in out Standard_Integer_Vectors.Link_to_Vector;
                permsys : in out Laur_Sys;
                vrblvl : in integer32 := 0 ) is

    perm : Standard_Integer_Vectors.Link_to_Vector;

  begin
    if vrblvl > 0
     then put_line("-> in drivers_for_static_lifting.Compute_Mixture 2 ...");
    end if;
    if compmix
     then Compute_Mixture(sup,mix,perm);
     else perm := Compute_Permutation(n,mix.all,sup);
    end if;
    permsys := Permute(permsys,perm);
    new_line(file);
    put(file,"TYPE OF MIXTURE : "); put(file,"#supports : ");
    put(file,mix'last,1);
    put(file,"  occurrences : "); put(file,mix);
    if mix'last /= n then
      new_line(file);
      put(file,"  the permutation : "); put(file,perm);
      new_line(file);
    end if;
    new_line(file);
    Clear(perm);
  end Compute_Mixture;

  function Expand ( mix : Vector;
                    points : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return Arrays_of_Integer_Vector_Lists.Array_of_Lists is

  -- DESCRIPTION :
  --   Returns a tuple of expanded lists, according to the type of mixture.

    sum : constant integer32 := Standard_Integer_Vectors.Sum(mix);
    res : Arrays_of_Integer_Vector_Lists.Array_of_Lists(1..sum);
    cnt : integer32 := 0;

  begin
    for i in mix'range loop
      for j in 1..mix(i) loop
        cnt := cnt + 1;
        res(cnt) := points(i);
      end loop;
    end loop;
    return res;
  end Expand;

  procedure Write_Cardinalities
               ( file : in file_type;
                 mix,card : in Standard_Integer_Vectors.Vector ) is
 
  begin
    new_line(file);
    put_line(file,"CARDINALITIES OF THE LIFTED FACES :");
    new_line(file);
    for i in card'range loop
      put(file,"  #lifted "); put(file,mix(i),1);
      put(file,"-faces of polytope "); put(file,i,1);
      put(file," : "); put(file,card(i),1); new_line(file);
    end loop;
  end Write_Cardinalities;

  procedure Integer_Create_Mixed_Cells
              ( n : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                lifted : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                mixsub : in out Integer_Mixed_Subdivisions.Mixed_Subdivision;
                vrblvl : in integer32 := 0 ) is

    use Integer_Faces_of_Polytope;
    afa : Array_of_Faces(mix'range);
    nbsucc,nbfail : Standard_Floating_Vectors.Vector(mix'range)
                  := (mix'range => 0.0);
    timer : timing_widget;

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_static_lifting.");
      put_line("Integer_Create_Mixed_Cells 1 ...");
    end if;
    for i in afa'range loop
      afa(i) := Create_Lower(mix(i),n+1,lifted(i));
    end loop;
    Create_CS(n,mix,afa,lifted,nbsucc,nbfail,mixsub);
  end Integer_Create_Mixed_Cells;

  procedure Integer_Create_Mixed_Cells
              ( file : in file_type; n : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                report : in boolean;
                lifted : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                mixsub : in out Integer_Mixed_Subdivisions.Mixed_Subdivision;
                vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   The pruning algorithm will be applied to compute the mixed cells.

    use Integer_Faces_of_Polytope;
    afa : Array_of_Faces(mix'range);
    cardafa : Standard_Integer_Vectors.Vector(mix'range);
    nbsucc,nbfail : Standard_Floating_Vectors.Vector(mix'range)
                  := (mix'range => 0.0);
    timer : timing_widget;
    tmv,nbcells : natural32;

    procedure Write_Cell ( mic : in Integer_Mixed_Subdivisions.Mixed_Cell;
                           continue : out boolean ) is

      vol : natural32;

    begin
      nbcells := nbcells + 1;
      put(file,"Cell no. "); put(file,nbcells,1); put_line(file," : ");
      put(file," normal to cell : "); put(file,mic.nor); new_line(file);
      put_line(file," the points in the cell : ");
      for k in mic.pts'range loop
        put(file,"  component "); put(file,k,1); put(file," with ");
        put(file,Length_Of(mic.pts(k)),1); put_line(file," points :");
        put(file,mic.pts(k));
      end loop;
      vol := Mixed_Volume(n,mix,mic);
      put(file," with volume addition : ");
      put(file,tmv,1); put(file," + "); put(file,vol,1);
      tmv := tmv + vol; put(file," = "); put(file,tmv,1); new_line(file);
      continue := true;
    end Write_Cell;
    procedure Report_and_Create1_CS is new Gen1_Create_CS(Write_Cell);

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_static_lifting.");
      put_line("Integer_Create_Mixed_Cells 2 ...");
    end if;
    tstart(timer);
    for i in afa'range loop
      afa(i) := Create_Lower(mix(i),n+1,lifted(i));
    end loop;
    tstop(timer);
    for i in afa'range loop
      cardafa(i) := integer32(Length_Of(afa(i)));
    end loop;
    Write_Cardinalities(file,mix,cardafa);
    new_line(file);
    print_times(file,timer,"Creation of the faces of lower hull");
    new_line(file);
    put_line(file,"PRUNING TO EXTRACT THE MIXED CELLS :");
    tstart(timer);
    if report
     then nbcells := 0; tmv := 0;
          Report_and_Create1_CS(n,mix,afa,lifted,nbsucc,nbfail,mixsub);
     else Create_CS(n,mix,afa,lifted,nbsucc,nbfail,mixsub);
    end if;
    tstop(timer);
    Pruning_Statistics(file,nbsucc,nbfail);
    new_line(file);
    print_times(file,timer,"Pruning for Mixed Cells");
  end Integer_Create_Mixed_Cells;

  procedure Floating_Create_Mixed_Cells
              ( file : in file_type;
                n : in integer32; mix : in Standard_Integer_Vectors.Vector;
                fltsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                lilifu : in Standard_Floating_VecVecs.Link_to_VecVec;
                lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                fltsub : in out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                vrblvl : in integer32 := 0 ) is

    use Standard_Floating_VecVecs;
    use Floating_Faces_of_Polytope;
    timer : timing_widget;
    tol : constant double_float := 1.0E-12;
    supfa,liffaces : Array_of_Faces(mix'range);
    cardafa : Standard_Integer_Vectors.Vector(mix'range);
    nbsucc,nbfail : Standard_Floating_Vectors.Vector(mix'range)
                  := (mix'range => 0.0);

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_static_lifting.");
      put_line("Floating_Create_Mixed_Cells ...");
    end if;
    tstart(timer);
    if lilifu /= null then
      for i in supfa'range loop
        supfa(i) := Create(mix(i),n,fltsup(i),tol);
        liffaces(i) := Linear_Lift(supfa(i),lilifu(i).all);
      end loop;
    elsif mix'length > 1 then
      for i in liffaces'range loop
        liffaces(i) := Create_Lower(mix(i),n+1,lifsup(i),tol);
      end loop;
    else
      liffaces(1) := Create_Lower_Facets(n+1,lifsup(1),tol);
    end if;
    tstop(timer);
    for i in liffaces'range loop
      cardafa(i) := integer32(Length_Of(liffaces(i)));
    end loop;
    Write_Cardinalities(file,mix,cardafa);
    new_line(file);
    print_times(file,timer,"Creation of the faces of lower hull");
    new_line(file);
    put_line(file,"PRUNING TO EXTRACT THE MIXED CELLS :");
    tstart(timer);
    Create(n,mix,liffaces,lifsup,tol,nbsucc,nbfail,fltsub);
    tstop(timer);
    Pruning_Statistics(file,nbsucc,nbfail);
    new_line(file);
    print_times(file,timer,"Pruning for Mixed Cells.");
  end Floating_Create_Mixed_Cells;

  procedure Integer_Volume_Computation
              ( n : in integer32; mix : in Standard_Integer_Vectors.Vector;
                compmisu : in boolean;
                lifpts : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                mixsub : in out Integer_Mixed_Subdivisions.Mixed_Subdivision;
                mv : out natural32; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put("-> in drivers_for_static_lifting.");
      put_line("Integer_Volume_Computation 1 ...");
    end if;
    if not compmisu then
      Mixed_Volume(n,mix,mixsub,mv);
    else
      mixsub := Integer_Mixed_Subdivisions.Create(lifpts,mixsub);
      Mixed_Volume(n,mix,mixsub,mv);
    end if;
  end Integer_Volume_Computation;

  procedure Integer_Volume_Computation
              ( file : in file_type;
                n : in integer32; mix : in Standard_Integer_Vectors.Vector;
                compmisu : in boolean;
                lifpts : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                mixsub : in out Integer_Mixed_Subdivisions.Mixed_Subdivision;
                mv : out natural32; vrblvl : in integer32 := 0 ) is

    timer : timing_widget;
    mixvol : natural32;

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_static_lifting.");
      put_line("Integer_Volume_Computation 2 ...");
    end if;
    new_line(file);
    put_line(file,"VOLUMES OF MIXED CELLS :");
    new_line(file);
    tstart(timer);
    put(file,natural32(n),mix,mixsub,mixvol);
    tstop(timer);
    put(file,"The total mixed volume equals ");
    put(file,mixvol,1); new_line(file);
    new_line(file);
    print_times(file,timer,"Volume computation of mixed cells");
    if compmisu then
      mixsub := Integer_Mixed_Subdivisions.Create(lifpts,mixsub);
      new_line(file);
      put_line(file,"VOLUMES OF MIXED CELLS, AFTER REFINEMENT :");
      new_line(file);
      tstart(timer);
      put(file,natural32(n),mix,mixsub,mixvol);
      tstop(timer);
      put(file,"The total mixed volume equals ");
      put(file,mixvol,1); new_line(file);
      new_line(file);
      print_times(file,timer,"Volume computation of mixed cells");
    end if;
    mv := mixvol;
  end Integer_Volume_Computation;

  procedure Floating_Volume_Computation
              ( n : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                mixsub : in out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                mv : out natural32;
                multprec_hermite : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    res : natural32 := 0;
    use Floating_Mixed_Subdivisions;
    tmp : Mixed_Subdivision := mixsub;
    mic : Mixed_Cell;
    vol : natural32;

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_static_lifting.");
      put_line("Floating_Volume_Computation 1 ...");
    end if;
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      Mixed_Volume(n,mix,mic,vol,multprec_hermite);
      res := res + vol;
      tmp := Tail_Of(tmp);
    end loop;
    mv := res;
  end Floating_Volume_Computation;

  procedure Floating_Volume_Computation
              ( file : in file_type; n : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                mixsub : in out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                mv : out natural32;
                multprec_hermite : in boolean := false;
                vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   Wraps a timer around the computation of the volumes of all cells.
  --   In case stlb is nonzero, the stable mixed volume is also computed.

    timer : timing_widget;

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_static_lifting.");
      put_line("Floating_Volume_Computation 2 ...");
    end if;
    new_line(file);
    put_line(file,"THE MIXED SUBDIVISION : ");
    new_line(file);
    tstart(timer);
    Floating_Mixed_Subdivisions_io.put
      (file,natural32(n),mix,mixsub,mv,multprec_hermite);
    tstop(timer);
    put(file,"the mixed volume : ");
    put(file,mv,1); new_line(file);
    new_line(file);
    print_times(file,timer,"Volume computation of mixed cells");
  end Floating_Volume_Computation;

  procedure Floating_Volume_Computation
              ( n : in integer32; stlb : in double_float;
                mix : in Standard_Integer_Vectors.Vector;
                mixsub : in out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                mv,smv,tmv : out natural32;
                multprec_hermite : in boolean := false;
                vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   Computes the mixed volume of the given mixed-cell configuration.
  --   In case stlb is nonzero, the stable mixed volume is also computed.

    use Floating_Mixed_Subdivisions;
    tmp : Mixed_Subdivision := mixsub;
    mic : Mixed_Cell;
    vol : natural32;

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_static_lifting.");
      put_line("Floating_Volume_Computation 3 ...");
    end if;
    if stlb = 0.0 then
      Floating_Volume_Computation(n,mix,mixsub,mv);
      smv := 0; tmv := 0;
    else
      mv := 0; smv := 0; tmv := 0;
      while not Is_Null(tmp) loop
        mic := Head_Of(tmp);
        Mixed_Volume(n,mix,mic,vol,multprec_hermite);
        if Is_Original(mic,stlb) then
	  mv := mv + vol;
	  smv := smv + vol;
        elsif Is_Stable(mic.nor.all,stlb,mic.pts.all) then
          smv := smv + vol;
        else
          tmv := tmv + vol;
        end if;
        tmp := Tail_Of(tmp);
      end loop;
    end if;
  end Floating_Volume_Computation;

  procedure Floating_Volume_Computation
              ( file : in file_type;
                n : in integer32; stlb : in double_float;
                mix : in Standard_Integer_Vectors.Vector;
                mixsub : in out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                mv,smv,tmv : out natural32;
                multprec_hermite : in boolean := false;
                vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   Wraps a timer around the computation of the volumes of all cells.
  --   In case stlb is nonzero, the stable mixed volume is also computed.

    timer : timing_widget;

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_static_lifting.");
      put_line("Floating_Volume_Computation 4 ...");
    end if;
    new_line(file);
    put_line(file,"THE MIXED SUBDIVISION : ");
    new_line(file);
    tstart(timer);
    if stlb /= 0.0 then
      Floating_Mixed_Subdivisions_io.put
        (file,natural32(n),stlb,mix,mixsub,mv,smv,tmv,multprec_hermite);
    else
      Floating_Mixed_Subdivisions_io.put
        (file,natural32(n),mix,mixsub,mv,multprec_hermite);
      smv := 0; tmv := 0;
    end if;
    tstop(timer);
    put(file,"common mixed volume : ");
    put(file,mv,1); new_line(file);
    if smv /= 0 then
      put(file,"stable mixed volume : "); put(file,smv,1); new_line(file);
      put(file," total mixed volume : "); put(file,tmv,1); new_line(file);
    end if;
    new_line(file);
    print_times(file,timer,"Volume computation of mixed cells");
  end Floating_Volume_Computation;

  procedure Integer_Polyhedral_Homotopy_Continuation
               ( file : in file_type; contrep : in boolean;
                 n : in integer32; mix : in Standard_Integer_Vectors.Vector;
                 q : in out Poly_Sys; qsols : in out Solution_List;
                 lifted : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 mixsub : in Integer_Mixed_Subdivisions.Mixed_Subdivision;
                 vrblvl : in integer32 := 0 ) is

    timer : timing_widget;
    lifted_lq,lq : Laur_Sys(q'range);
    h : Eval_Coeff_Laur_Sys(q'range);
    c : Standard_Complex_VecVecs.VecVec(h'range);
    e : Exponent_Vectors_Array(h'range);
    j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
    m : Mult_Factors(j'range(1),j'range(2));

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_static_lifting.");
      put_line("Integer_Polyhedral_Homotopy_Continuation 1 ...");
    end if;
    new_line(file);
    put_line(file,"POLYHEDRAL HOMOTOPY CONTINUATION :");
    lq := Polynomial_to_Laurent_System(q);
    lifted_lq := Perform_Lifting(n,mix,lifted,lq);
    Clear(lq);
    lq := Eval(lifted_lq,Create(1.0),n+1);
    q := Laurent_to_Polynomial_System(lq);
    h := Create(lq);
    for i in c'range loop
      declare
        coeff_lq : constant Standard_Complex_Vectors.Vector := Coeff(lq(i));
      begin
        c(i) := new Standard_Complex_Vectors.Vector(coeff_lq'range);
        for k in coeff_lq'range loop
          c(i)(k) := coeff_lq(k);
        end loop;
      end;
    end loop;
    e := Create(lq);
    Create(lq,j,m);
    tstart(timer);
    if contrep
     then --Mixed_Solve(file,lifted_lq,mix,mixsub,qsols);
          Mixed_Solve(file,lifted_lq,lifted,h,c,e,j,m,mix,mixsub,qsols);
     else --Mixed_Solve(lifted_lq,mix,mixsub,qsols);
          Mixed_Solve(lifted_lq,lifted,h,c,e,j,m,mix,mixsub,qsols);
    end if;
    Clear(h); Clear(j); Clear(m); Clear(e);
    Standard_Complex_VecVecs.Clear(c);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Polyhedral Continuation");
  end Integer_Polyhedral_Homotopy_Continuation;

  procedure Integer_Polyhedral_Homotopy_Continuation
              ( file : in file_type; contrep : in boolean;
                n : in integer32; mix : in Standard_Integer_Vectors.Vector;
                q : in out Laur_Sys; qsols : in out Solution_List;
                lifted : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                mixsub : in Integer_Mixed_Subdivisions.Mixed_Subdivision;
                vrblvl : in integer32 := 0 ) is

    timer : timing_widget;
    lifted_q : Laur_Sys(q'range);
    h : Eval_Coeff_Laur_Sys(q'range);
    c : Standard_Complex_VecVecs.VecVec(h'range);
    e : Exponent_Vectors_Array(h'range);
    j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
    m : Mult_Factors(j'range(1),j'range(2));

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_static_lifting.");
      put_line("Integer_Polyhedral_Homotopy_Continuation 2 ...");
    end if;
    new_line(file);
    put_line(file,"POLYHEDRAL HOMOTOPY CONTINUATION :");
    lifted_q := Perform_Lifting(n,mix,lifted,q);
    Clear(q);
    q := Eval(lifted_q,Create(1.0),n+1);
    h := Create(q);
    for i in c'range loop
      declare
        coeff_lq : constant Standard_Complex_Vectors.Vector := Coeff(q(i));
      begin
        c(i) := new Standard_Complex_Vectors.Vector(coeff_lq'range);
        for k in coeff_lq'range loop
          c(i)(k) := coeff_lq(k);
        end loop;
      end;
    end loop;
    e := Create(q);
    Create(q,j,m);
    tstart(timer);
    if contrep
     then --Mixed_Solve(file,lifted_q,mix,mixsub,qsols);
          Mixed_Solve(file,lifted_q,lifted,h,c,e,j,m,mix,mixsub,qsols);
     else --Mixed_Solve(lifted_q,mix,mixsub,qsols);
          Mixed_Solve(lifted_q,lifted,h,c,e,j,m,mix,mixsub,qsols);
    end if;
    Clear(h); Clear(j); Clear(m); Clear(e);
    Standard_Complex_VecVecs.Clear(c);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Polyhedral Continuation");
  end Integer_Polyhedral_Homotopy_Continuation;

  procedure Floating_Polyhedral_Homotopy_Continuation
              ( file : in file_type; nt : in integer32; contrep : in boolean;
                n : in integer32; mix : in Standard_Integer_Vectors.Vector;
                q : in Poly_Sys; qsols : in out Solution_List;
                lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                fltsub : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                vrblvl : in integer32 := 0 ) is

    lq : constant Laur_Sys(q'range) := Polynomial_to_Laurent_System(q);
    h : Eval_Coeff_Laur_Sys(q'range);
    c : Standard_Complex_VecVecs.VecVec(h'range);
    e : Exponent_Vectors_Array(h'range);
    j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
    m : Mult_Factors(j'range(1),j'range(2));
    timer : timing_widget;

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_static_lifting.");
      put_line("Floating_Polyhedral_Homotopy_Continuation 1 ...");
    end if;
    new_line(file);
    if nt = 0 then
      put_line(file,"POLYHEDRAL HOMOTOPY CONTINUATION without tasking :");
    else
      put(file,"POLYHEDRAL HOMOTOPY CONTINUATION with ");
      put(file,nt,1); put_line(file," tasks :");
    end if;
    h := Create(lq);
    for i in c'range loop
      declare
        coeff_lq : constant Standard_Complex_Vectors.Vector := Coeff(lq(i));
      begin
        c(i) := new Standard_Complex_Vectors.Vector(coeff_lq'range);
        for k in coeff_lq'range loop
          c(i)(k) := coeff_lq(k);
        end loop;
      end;
    end loop;
    e := Create(lq);
    Create(lq,j,m);
    tstart(timer);
    if nt > 0 then
      Silent_Multitasking_Path_Tracker
        (lq,nt,n,mix'last,mix,lifsup,fltsub,h,c,e,j,m,qsols);
    else
      if contrep
       then Mixed_Solve(file,lq,lifsup,h,c,e,j,m,mix,fltsub,qsols);
       else Mixed_Solve(lq,lifsup,h,c,e,j,m,mix,fltsub,qsols);
      end if;
    end if;
    Clear(h); Clear(j); Clear(m); Clear(e);
    Standard_Complex_VecVecs.Clear(c);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Polyhedral Homotopy Continuation");
  end Floating_Polyhedral_Homotopy_Continuation;

  procedure Floating_Polyhedral_Homotopy_Continuation
              ( file : in file_type; nt : in integer32; contrep : in boolean;
                n : in integer32; mix : in Standard_Integer_Vectors.Vector;
                q : in Laur_Sys; qsols : in out Solution_List;
                lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                fltsub : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                vrblvl : in integer32 := 0 ) is

    h : Eval_Coeff_Laur_Sys(q'range);
    c : Standard_Complex_VecVecs.VecVec(h'range);
    e : Exponent_Vectors_Array(h'range);
    j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
    m : Mult_Factors(j'range(1),j'range(2));
    timer : timing_widget;

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_static_lifting.");
      put_line("Floating_Polyhedral_Homotopy_Continuation 2 ...");
    end if;
    new_line(file);
    if nt = 0 then
      put_line(file,"POLYHEDRAL HOMOTOPY CONTINUATION without tasking :");
    else
      put(file,"POLYHEDRAL HOMOTOPY CONTINUATION with ");
      put(file,nt,1); put_line(file," tasks :");
    end if;
    h := Create(q);
    for i in c'range loop
      declare
        coeff_lq : constant Standard_Complex_Vectors.Vector := Coeff(q(i));
      begin
        c(i) := new Standard_Complex_Vectors.Vector(coeff_lq'range);
        for k in coeff_lq'range loop
          c(i)(k) := coeff_lq(k);
        end loop;
      end;
    end loop;
    e := Create(q);
    Create(q,j,m);
    tstart(timer);
    if nt > 0 then
      Silent_Multitasking_Path_Tracker
        (q,nt,n,mix'last,mix,lifsup,fltsub,h,c,e,j,m,qsols);
    else
      if contrep
       then Mixed_Solve(file,q,lifsup,h,c,e,j,m,mix,fltsub,qsols);
       else Mixed_Solve(q,lifsup,h,c,e,j,m,mix,fltsub,qsols);
      end if;
    end if;
    Clear(h); Clear(j); Clear(m); Clear(e);
    Standard_Complex_VecVecs.Clear(c);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Polyhedral Homotopy Continuation");
  end Floating_Polyhedral_Homotopy_Continuation;

  procedure Data_Management
              ( file : in file_type;
                mix : in out Standard_Integer_Vectors.Link_to_Vector;
                imixsub : out Integer_Mixed_Subdivisions.Mixed_Subdivision;
                fmixsub : out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                compmisu,compmix,fltlif : out boolean;
                vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   This procedure allows to use previously computed mixed subdivisions.

    ans : character;
    m : natural32 := 0;

  begin
    if vrblvl > 0
     then put_line("-> in drivers_for_static_lifting.Data_Management ...");
    end if;
    new_line;
    put("Do you already have a mixed subdivision ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put("Induced by integer or floating-point lifting (i/f) ");
      Ask_Alternative(ans,"if");
      fltlif := (ans = 'f');
      declare
        insubft : file_type;
        nn : natural32 := 0;
      begin
        put_line("Reading the name of the input file.");
        Read_Name_and_Open_File(insubft);
        if ans = 'f'
         then get(insubft,nn,m,mix,fmixsub);
         else get(insubft,nn,m,mix,imixsub);
        end if;
        Close(insubft);
        new_line(file);
        put_line(file,"Mixed subdivision supplied by user.");
        new_line(file);
        compmisu := false; compmix := false;
      exception
        when DATA_ERROR
          => put_line("Data not in correct format.  Will ignore it...");
             Close(insubft);
      end;
    else
      compmisu := true;
      put("Do you want to enforce a type mixture ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        put("Give number of different supports : "); get(m);
        put("Give vector indicating occurrences : ");
        get(m,mix);
        compmix := false;
      else
        compmix := true;
      end if;
    end if;
  end Data_Management;

  procedure Write_Results
               ( file : in file_type; n : in integer32;
                 ranstart : in boolean; qft,solsft : in out file_type;
                 q : in Poly_Sys; qsols,qsols0 : in Solution_List ) is

  -- DESCRIPTION :
  --   Writes the results of the polyhedral continuation to the file,
  --   and to the other separate files qft and solft.

  -- REQUIRED : all needed files must be opened for output.

  -- ON ENTRY :
  --   file      to write diagnostics on;
  --   n         length of the solution vectors;
  --   ranstart  if true, then a random coefficient system was created
  --             and it will be written to the separate file qft,
  --             otherwise only to the solutions will be written to
  --             the separate file solsft;
  --   qft       to write the random coefficient system to,
  --   solsft    to write the start solutions to;
  --   q         random coefficient start system;
  --   qsols     solutions of q;
  --   qsols0    solutions of q with zero components.

  -- ON RETURN :
  --   qft       closed when q was written to it;
  --   solsft    closed when qsols was written to it.

  begin
    new_line(file);
    put_line(file,"THE RANDOM COEFFICIENT START SYSTEM :");
    new_line(file);
    put_line(file,q);
    new_line(file);
    put_line(file,"THE START SOLUTIONS :");
    new_line(file);
    put(file,Length_Of(qsols),natural32(q'length),qsols);
    if ranstart then
      new_line(qft);
      put_line(qft,"THE SOLUTIONS : ");
      new_line(qft);
      put(qft,Length_Of(qsols),natural32(n),qsols);
      if not Is_Null(qsols0) then
        new_line(qft);
        put_line(qft,"THE SOLUTIONS WITH ZERO COMPONENTS : ");
        new_line(qft);
        put(qft,Length_Of(qsols0),natural32(n),qsols0);
      end if;
      Close(qft);
    else
      put(solsft,Length_Of(qsols),natural32(n),qsols);
      if not Is_Null(qsols0)
       then put(solsft,Length_Of(qsols0),natural32(n),qsols0);
      end if;
      Close(solsft);
    end if;
  end Write_Results;

  procedure Write_Results
               ( file : in file_type; n : in integer32;
                 ranstart : in boolean; qft,solsft : in out file_type;
                 q : in Laur_Sys; qsols,qsols0 : in Solution_List ) is

  -- DESCRIPTION :
  --   Writes the results of the polyhedral continuation to the file,
  --   and to the other separate files qft and solft.

  -- REQUIRED : all needed files must be opened for output.

  -- ON ENTRY :
  --   file      to write diagnostics on;
  --   n         length of the solution vectors;
  --   ranstart  if true, then a random coefficient system was created
  --             and it will be written to the separate file qft,
  --             otherwise only to the solutions will be written to
  --             the separate file solsft;
  --   qft       to write the random coefficient system to,
  --   solsft    to write the start solutions to;
  --   q         random coefficient start system;
  --   qsols     solutions of q;
  --   qsols0    solutions of q with zero components.

  -- ON RETURN :
  --   qft       closed when q was written to it;
  --   solsft    closed when qsols was written to it.

  begin
    new_line(file);
    put_line(file,"THE RANDOM COEFFICIENT START SYSTEM :");
    new_line(file);
    put_line(file,q);
    new_line(file);
    put_line(file,"THE START SOLUTIONS :");
    new_line(file);
    put(file,Length_Of(qsols),natural32(q'length),qsols);
    if ranstart then
      new_line(qft);
      put_line(qft,"THE SOLUTIONS : ");
      new_line(qft);
      put(qft,Length_Of(qsols),natural32(n),qsols);
      if not Is_Null(qsols0) then
        new_line(qft);
        put_line(qft,"THE SOLUTIONS WITH ZERO COMPONENTS : ");
        new_line(qft);
        put(qft,Length_Of(qsols0),natural32(n),qsols0);
      end if;
      Close(qft);
    else
      put(solsft,Length_Of(qsols),natural32(n),qsols);
      if not Is_Null(qsols0)
       then put(solsft,Length_Of(qsols0),natural32(n),qsols0);
      end if;
      Close(solsft);
    end if;
  end Write_Results;

  procedure Prompt_for_File
               ( misufile : out boolean;
                 outsubft : in out file_type ) is

  -- DESCRIPTION :
  --   Prompts the user regarding the output file before the
  --   computation of a regular mixed-cell configuration.

  -- ON RETURN :
  --   misufile  indicates whether the mixed-cell configuration has
  --             to be written on a separate file or not;
  --   outsubft  if misufile, then on return, this file pointer
  --             has been assigned to some output file.

    ans : character;

  begin
    put("Do you want the mixed cells on separate file ? (y/n) ");
    Ask_Yes_or_No(ans);
    misufile := (ans = 'y');
    if misufile then
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(outsubft);
    end if;
  end Prompt_for_File;

  procedure Prompt_for_Options
               ( report,misufile : out boolean;
                 outsubft : in out file_type ) is

  -- DESCRIPTION :
  --   Prompts the user for options during the computation
  --   of a regular mixed-cell configuration.

  -- ON RETURN :
  --   report    indicates whether intermediate output during the
  --             computation of the mixed cells is requested;
  --   misufile  indicates whether the mixed-cell configuration has
  --             to be written on a separate file or not;
  --   outsubft  if misufile, then on return, this file pointer
  --             has been assigned to some output file.

    ans : character;

  begin
    new_line;
    put("Do you want intermediate output on file, during creation ? (y/n) ");
    Ask_Yes_or_No(ans);
    report := (ans = 'y');
    Prompt_for_File(misufile,outsubft);
  end Prompt_for_Options;

  procedure Polyhedral_Homotopies_with_Integer_Lifting
               ( file : in file_type; compmisu,byebye : in boolean;
                 n : in integer32; p : in Poly_Sys;
                 mix : in Standard_Integer_Vectors.Vector;
                 points : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                -- mixpts : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 lif : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 mcc : in out Integer_Mixed_Subdivisions.Mixed_Subdivision;
                 q : out Poly_Sys; qsols : out Solution_List;
                 gft,solsft : in out file_type;
                 tosolve,ranstart,contrep : out boolean;
                 mv : out natural32; vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   After rearranging the input supports to compute the mixture,
  --   a mixed-cell configuration induced by the given integer lifting
  --   will be computed and a random coefficient system will be solved,
  --   depending on the interaction with the user.

    misufile,report : boolean;
    outsubft : file_type;

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_static_lifting.");
      put_line("Polyhedral_Homotopies_with_Integer_Lifting 1 ...");
    end if;
    if compmisu
     then Prompt_for_Options(report,misufile,outsubft);
     else lif := Induced_Lifting(n,mix,points,mcc);
    end if;
    Driver_for_Coefficient_System
      (file,p,0,byebye,q,gft,solsft,tosolve,ranstart,contrep);
    if compmisu
     -- then Integer_Create_Mixed_Cells(file,n,mix,report,mixpts,lif,mcc);
     then Integer_Create_Mixed_Cells(file,n,mix,report,lif,mcc,vrblvl-1);
    end if;
    if not Integer_Mixed_Subdivisions.Is_Null(mcc) then
      Integer_Volume_Computation(file,n,mix,compmisu,lif,mcc,mv,vrblvl-1);
      if compmisu and then misufile
       then put(outsubft,natural32(n),mix,mcc);
      end if;
      if tosolve
       then Integer_Polyhedral_Homotopy_Continuation
               (file,contrep,n,mix,q,qsols,lif,mcc,vrblvl-1);
      end if;
    end if;
  end Polyhedral_Homotopies_with_Integer_Lifting;

  procedure Polyhedral_Homotopies_with_Integer_Lifting
               ( file : in file_type; compmisu,byebye : in boolean;
                 n : in integer32; p : in Laur_Sys;
                 mix : in Standard_Integer_Vectors.Vector;
                 points : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                -- mixpts : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 lif : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 mcc : in out Integer_Mixed_Subdivisions.Mixed_Subdivision;
                 q : out Laur_Sys; qsols : out Solution_List;
                 gft,solsft : in out file_type;
                 tosolve,ranstart,contrep : out boolean;
                 mv : out natural32; vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   After rearranging the input supports to compute the mixture,
  --   a mixed-cell configuration induced by the given integer lifting
  --   will be computed and a random coefficient system will be solved,
  --   depending on the interaction with the user.
  --   This version allows Laurent polynomial systems.

    misufile,report : boolean;
    outsubft : file_type;

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_static_lifting.");
      put_line("Polyhedral_Homotopies_with_Integer_Lifting 2 ...");
    end if;
    if compmisu
     then Prompt_for_Options(report,misufile,outsubft);
     else lif := Induced_Lifting(n,mix,points,mcc);
    end if;
    Driver_for_Coefficient_System
      (file,p,0,byebye,q,gft,solsft,tosolve,ranstart,contrep);
    if compmisu
     -- then Integer_Create_Mixed_Cells(file,n,mix,report,mixpts,lif,mcc);
     then Integer_Create_Mixed_Cells(file,n,mix,report,lif,mcc,vrblvl-1);
    end if;
    if not Integer_Mixed_Subdivisions.Is_Null(mcc) then
      Integer_Volume_Computation(file,n,mix,compmisu,lif,mcc,mv,vrblvl-1);
      if compmisu and then misufile
       then put(outsubft,natural32(n),mix,mcc);
      end if;
      if tosolve
       then Integer_Polyhedral_Homotopy_Continuation
              (file,contrep,n,mix,q,qsols,lif,mcc,vrblvl-1);
      end if;
    end if;
  end Polyhedral_Homotopies_with_Integer_Lifting;

  procedure Polyhedral_Homotopies_with_Float_Lifting
              ( file : in file_type; nt : in integer32;
                compmisu,byebye : in boolean;
                n : in integer32; p : in Poly_Sys; stlb : in double_float;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                points : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mixpts : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                flili : in Standard_Floating_VecVecs.Link_to_VecVec;
                lif : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mcc : in out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                q : out Poly_Sys; qsols,qsols0 : out Solution_List;          
                gft,solsft : in out file_type;
                tosolve,ranstart,contrep : out boolean;
                mv,smv,tmv : out natural32; vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   This is the analogue driver for floating-point lifting functions.
  --   The user is prompted to make some choices before the computations.

    misufile : boolean;
    outsubft : file_type;
   -- sp : Poly_Sys(p'range);
    lq : Laur_Sys(q'range);
    bnd,bnd2 : double_float;
    lifted : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mixpts'range);
    orgmcc,stbmcc : Floating_Mixed_Subdivisions.Mixed_Subdivision;
    orgcnt,stbcnt : natural32;
    use Floating_Mixed_Subdivisions;

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_static_lifting.");
      put_line("Polyhedral_Homotopies_with_Float_Lifting 1 ...");
    end if;
    if compmisu then
      new_line;
      Prompt_for_File(misufile,outsubft);
      bnd := stlb;
    else 
      lif := Occurred_Lifting(n,mix.all,points,mcc);
      if stlb = 0.0 then
        lifted := Lifted_Supports(mix'last,mcc);
        bnd := Lifting_Bound(lifted);
        put("The extracted lifting bound :"); put(bnd); new_line;
        bnd2 := Lifting_Bound(p);
        put(" the computed lifting bound :"); put(bnd2); new_line;
        if bnd < bnd2 then
          bnd := 0.0;  -- reject extracted stable lifting bound
        end if;
      else
        bnd := stlb;   -- use variable bnd to hold stable lifting bound
      end if;
     -- sp := Select_Lifted(p,mix.all,lif);
    end if;
    Driver_for_Coefficient_System
      (file,p,0,byebye,q,gft,solsft,tosolve,ranstart,contrep);
    if compmisu then
      Floating_Create_Mixed_Cells
        (file,n,mix.all,points,flili,lif,mcc,vrblvl-1);
      if misufile
       then put(outsubft,natural32(n),mix.all,mcc);
      end if;
    end if;
    if not Floating_Mixed_Subdivisions.Is_Null(mcc) then
      Floating_Volume_Computation
        (file,n,bnd,mix.all,mcc,mv,smv,tmv,vrblvl=>vrblvl-1);
      if tosolve then
        if bnd = 0.0 then
          Floating_Polyhedral_Homotopy_Continuation
              (file,nt,contrep,n,mix.all,q,qsols,lif,mcc,vrblvl-1);
        else
          Split_Original_Cells(mcc,bnd,orgmcc,stbmcc,orgcnt,stbcnt);
          new_line(file);
          put(file," Number of cells without artificial origin : ");
          put(file,orgcnt,1); new_line(file);
          put(file,"#extra stable cells with artificial origin : ");
          put(file,stbcnt,1); new_line(file);
          Floating_Polyhedral_Homotopy_Continuation
              (file,nt,contrep,n,mix.all,q,qsols,lif,orgmcc,vrblvl-1);
          lq := Polynomial_to_Laurent_System(q);
          if contrep then
            Reporting_Polyhedral_Continuation
              (file,lq,bnd,mix,lif,stbmcc,qsols0);
	  else
	    Silent_Polyhedral_Continuation(lq,bnd,mix,lif,stbmcc,qsols0);
	  end if;
        end if;
      end if;
    end if;
  end Polyhedral_Homotopies_with_Float_Lifting;

  procedure Polyhedral_Homotopies_with_Float_Lifting
              ( file : in file_type; nt : in integer32;
                compmisu,byebye : in boolean;
                n : in integer32; p : in Laur_Sys; stlb : in double_float;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                points : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mixpts : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                flili : in Standard_Floating_VecVecs.Link_to_VecVec;
                lif : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mcc : in out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                q : out Laur_Sys; qsols,qsols0 : out Solution_List;          
                gft,solsft : in out file_type;
                tosolve,ranstart,contrep : out boolean;
                mv,smv,tmv : out natural32; vrblvl : in integer32 := 0 ) is

   -- DESCRIPTION :
   --   This is the analogue driver for floating-point lifting functions,
   --   but for Laurent polynomial systems.
   --   The user is prompted to make some choices before the computations.

    misufile : boolean;
    outsubft : file_type;
   -- sp : Laur_Sys(p'range);
    lifted : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mixpts'range);
    bnd,bnd2 : double_float;
    orgmcc,stbmcc : Floating_Mixed_Subdivisions.Mixed_Subdivision;
    orgcnt,stbcnt : natural32;
    use Floating_Mixed_Subdivisions;

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_static_lifting.");
      put_line("Polyhedral_Homotopies_with_Float_Lifting 2 ...");
    end if;
    if compmisu then
      new_line;
      Prompt_for_File(misufile,outsubft);
      bnd := stlb;
    else 
      lif := Occurred_Lifting(n,mix.all,points,mcc);
      if Is_Genuine_Laurent(p) then
        bnd := 0.0;
      elsif stlb = 0.0 then
        lifted := Lifted_Supports(mix'last,mcc);
        bnd := Lifting_Bound(lifted);
        put("The lifting bound : "); put(bnd); new_line;
        put("The extracted lifting bound :"); put(bnd); new_line;
        bnd2 := Lifting_Bound(p);
        put(" the computed lifting bound :"); put(bnd2); new_line;
        if bnd < bnd2
         then bnd := 0.0; -- reject extracted lifting bound
        end if;
      else
        bnd := stlb; -- hold in bnd the value of stable lifting bound
      end if;
     -- sp := Select_Lifted(p,mix.all,lif);
    end if;
    Driver_for_Coefficient_System
      (file,p,0,byebye,q,gft,solsft,tosolve,ranstart,contrep);
    if compmisu then
      Floating_Create_Mixed_Cells(file,n,mix.all,points,flili,lif,mcc,vrblvl-1);
      if misufile
       then put(outsubft,natural32(n),mix.all,mcc);
      end if;
    end if;
    if not Floating_Mixed_Subdivisions.Is_Null(mcc) then
      Floating_Volume_Computation
        (file,n,bnd,mix.all,mcc,mv,smv,tmv,vrblvl=>vrblvl-1);
      if tosolve then
        if bnd = 0.0 then
          Floating_Polyhedral_Homotopy_Continuation
              (file,nt,contrep,n,mix.all,q,qsols,lif,mcc,vrblvl-1);
        else
          Split_Original_Cells(mcc,bnd,orgmcc,stbmcc,orgcnt,stbcnt);
          new_line(file);
          put(file," Number of cells without artificial origin : ");
          put(file,orgcnt,1); new_line(file);
          put(file,"#extra stable cells with artificial origin : ");
          put(file,stbcnt,1); new_line(file);
          Floating_Polyhedral_Homotopy_Continuation
            (file,nt,contrep,n,mix.all,q,qsols,lif,orgmcc,vrblvl-1);
          if contrep then
            Reporting_Polyhedral_Continuation
              (file,q,bnd,mix,lif,stbmcc,qsols0);
	  else
	    Silent_Polyhedral_Continuation(q,bnd,mix,lif,stbmcc,qsols0);
	  end if;
        end if;
      end if;
    end if;
  end Polyhedral_Homotopies_with_Float_Lifting;

  procedure Prepare_Supports 
              ( file : in file_type; n : in integer32; p : in Poly_Sys;
                mix : in out Standard_Integer_Vectors.Link_to_Vector;
                sup : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                mixpts,ilifpts
                 : out Arrays_of_Integer_Vector_Lists.Link_to_Array_of_Lists;
                fpts,flifpts
                 : out Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
                ilili : out Standard_Integer_VecVecs.Link_to_VecVec;
                flili : out Standard_Floating_VecVecs.Link_to_VecVec;
                fltlif : in out boolean; stlb : out double_float;
                compmisu,compmix : in boolean; sp : out Poly_Sys;
                vrblvl : in integer32 := 0 ) is

  -- DESRIPTION :
  --   To prepare the supports for mixed-volume computation, the type
  --   of mixture is computed and superfluous points are omitted.

    permp : Poly_Sys(p'range) := p;
    mixpts1 : Arrays_of_Integer_Vector_Lists.Link_to_Array_of_Lists;

  begin
    if vrblvl > 0
     then put_line("-> in drivers_for_static_lifting.Prepare_Supports 1 ...");
    end if;
    if compmisu then
      Compute_Mixture(file,n,compmix,sup,mix,permp);
      mixpts1 := new Arrays_of_Integer_Vector_Lists.
                     Array_of_Lists'(Typed_Lists(mix.all,sup));
      Driver_for_Criterion(file,mixpts1.all);
      mixpts := new Arrays_of_Integer_Vector_Lists.
                    Array_of_Lists'(Expand(mix.all,mixpts1.all));
      Compute_Mixture(file,n,compmix,mixpts.all,mix,permp);
      mixpts := new Arrays_of_Integer_Vector_Lists.
                    Array_of_Lists'(Typed_Lists(mix.all,mixpts.all));

      ilifpts := new Arrays_of_Integer_Vector_Lists.
                     Array_of_Lists(mix'range);
      fpts := new Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
      flifpts
        := new Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);

      sp := Select_Terms(permp,mix.all,mixpts.all);
      new_line;
      Main_Lifting_Functions.Main_Polynomial
        (file,sp,mixpts.all,fltlif,stlb,fpts.all,ilifpts.all,
              flifpts.all,ilili,flili);
    else
      mixpts := new Arrays_of_Integer_Vector_Lists.
                    Array_of_Lists'(Typed_Lists(mix.all,sup));
      ilifpts := new Arrays_of_Integer_Vector_Lists.
                     Array_of_Lists(mix'range);
      fpts := new Arrays_of_Floating_Vector_Lists.
                  Array_of_Lists(sup'range);
      fpts.all := Floating_Integer_Convertors.Convert(sup);
      flifpts := new Arrays_of_Floating_Vector_Lists.
                     Array_of_Lists(mix'range);
      sp := p;
      stlb := 0.0;
    end if;
  end Prepare_Supports;

  procedure Prepare_Supports 
              ( file : in file_type; n : in integer32; p : in Laur_Sys;
                mix : in out Standard_Integer_Vectors.Link_to_Vector;
                sup : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                mixpts,ilifpts
                 : out Arrays_of_Integer_Vector_Lists.Link_to_Array_of_Lists;
                fpts,flifpts
                 : out Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
                ilili : out Standard_Integer_VecVecs.Link_to_VecVec;
                flili : out Standard_Floating_VecVecs.Link_to_VecVec;
                fltlif : in out boolean; stlb : out double_float;
                compmisu,compmix : in boolean; sp : out Laur_Sys;
                vrblvl : in integer32 := 0 ) is

  -- DESRIPTION :
  --   To prepare the supports for mixed-volume computation, the type
  --   of mixture is computed and superfluous points are omitted.
  --   This version allows Laurent systems on input.

    permp : Laur_Sys(p'range) := p;
    mixpts1 : Arrays_of_Integer_Vector_Lists.Link_to_Array_of_Lists;

  begin
    if vrblvl > 0
     then put_line("-> in drivers_for_static_lifting.Prepare_Supports 2 ...");
    end if;
    if compmisu then
      Compute_Mixture(file,n,compmix,sup,mix,permp);
      mixpts1 := new Arrays_of_Integer_Vector_Lists.
                     Array_of_Lists'(Typed_Lists(mix.all,sup));
      Driver_for_Criterion(file,mixpts1.all);
      mixpts := new Arrays_of_Integer_Vector_Lists.
                    Array_of_Lists'(Expand(mix.all,mixpts1.all));
      Compute_Mixture(file,n,compmix,mixpts.all,mix,permp);
      mixpts := new Arrays_of_Integer_Vector_Lists.
                    Array_of_Lists'(Typed_Lists(mix.all,mixpts.all));

      ilifpts := new Arrays_of_Integer_Vector_Lists.
                     Array_of_Lists(mix'range);
      fpts := new Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
      flifpts
        := new Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);

      sp := Select_Terms(permp,mix.all,mixpts.all);
      new_line;
      Main_Lifting_Functions.Main_Laurent
        (file,sp,mixpts.all,fltlif,stlb,fpts.all,ilifpts.all,
              flifpts.all,ilili,flili);
    else
      mixpts := new Arrays_of_Integer_Vector_Lists.
                    Array_of_Lists'(Typed_Lists(mix.all,sup));
      ilifpts := new Arrays_of_Integer_Vector_Lists.
                     Array_of_Lists(mix'range);
      fpts := new Arrays_of_Floating_Vector_Lists.
                  Array_of_Lists(sup'range);
      fpts.all := Floating_Integer_Convertors.Convert(sup);
      flifpts := new Arrays_of_Floating_Vector_Lists.
                     Array_of_Lists(mix'range);
      sp := p;
      stlb := 0.0;
    end if;
  end Prepare_Supports;

  procedure Driver_for_Mixed_Volume_Computation 
               ( file : in file_type; nt : in integer32;
                 p : in Poly_Sys; byebye : in boolean;
                 q : out Poly_Sys; qsols,qsols0 : out Solution_List;
                 mv,smv,tmv : out natural32;
                 vrblvl : in integer32 := 0 ) is

    welcome : constant string := "Mixed-Volume Computation by Static Lifting";
    mix : Standard_Integer_Vectors.Link_to_Vector;
    imixsub : Integer_Mixed_Subdivisions.Mixed_Subdivision;
    fmixsub : Floating_Mixed_Subdivisions.Mixed_Subdivision;
    compmix   : boolean;   -- compute the type of mixture
    compmisu  : boolean;   -- compute the subdivision
    tosolve   : boolean;   -- solve the system 
    ranstart  : boolean;   -- construct random coefficient start system
    contrep   : boolean;   -- intermediate output during continuation
    fltlif    : boolean;   -- floating valued lifting
    stlb : double_float;
    solsft,gft : file_type;
    n : constant integer32 := p'length;
    sp : Poly_Sys(p'range);
    totaltimer : timing_widget;
    points : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);
    ilili : Standard_Integer_VecVecs.Link_to_VecVec;
    flili : Standard_Floating_VecVecs.Link_to_VecVec;
    mixpts,ilifpts : Arrays_of_Integer_Vector_Lists.Link_to_Array_of_Lists;
    fpts,flifpts : Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_static_lifting.");
      put_line("Driver_for_Mixed_Volume_Computation 1 ...");
    end if;
    new_line; put_line(welcome);
    tstart(totaltimer);
    points := Create(p);
    Data_Management(file,mix,imixsub,fmixsub,compmisu,compmix,fltlif,vrblvl-1);
    Prepare_Supports(file,n,p,mix,points,mixpts,ilifpts,fpts,flifpts,
                     ilili,flili,fltlif,stlb,compmisu,compmix,sp,vrblvl-1);
    if fltlif then
      Polyhedral_Homotopies_with_Float_Lifting
        (file,nt,compmisu,byebye,n,sp,stlb,mix,fpts.all,mixpts.all,flili,
         flifpts.all,fmixsub,q,qsols,qsols0,gft,solsft,
         tosolve,ranstart,contrep,mv,smv,tmv,vrblvl-1);
    else
      Polyhedral_Homotopies_with_Integer_Lifting
       -- (file,compmisu,byebye,n,sp,mix.all,points,mixpts.all,ilifpts.all,
        (file,compmisu,byebye,n,sp,mix.all,points,ilifpts.all,
         imixsub,q,qsols,gft,solsft,tosolve,ranstart,contrep,mv,vrblvl-1);
     end if;
    if tosolve
     then Write_Results(file,n,ranstart,gft,solsft,q,qsols,qsols0);
    end if;
    tstop(totaltimer);
    new_line(file);
    print_times(file,totaltimer,"All Computations");
  end Driver_for_Mixed_Volume_Computation;

  procedure Driver_for_Mixed_Volume_Computation 
               ( file : in file_type; nt : in integer32;
                 p : in Laur_Sys; byebye : in boolean;
                 q : out Laur_Sys; qsols,qsols0 : out Solution_List;
                 mv,smv,tmv : out natural32;
                 vrblvl : in integer32 := 0 ) is

    welcome : constant string := "Mixed-Volume Computation by Static Lifting";
    mix : Standard_Integer_Vectors.Link_to_Vector;
    imixsub : Integer_Mixed_Subdivisions.Mixed_Subdivision;
    fmixsub : Floating_Mixed_Subdivisions.Mixed_Subdivision;
    compmix   : boolean;   -- compute the type of mixture
    compmisu  : boolean;   -- compute the subdivision
    tosolve   : boolean;   -- solve the system 
    ranstart  : boolean;   -- construct random coefficient start system
    contrep   : boolean;   -- intermediate output during continuation
    fltlif    : boolean;   -- floating valued lifting function
    stlb : double_float;
    solsft,gft : file_type;
    n : constant integer32 := p'length;
    sp : Laur_Sys(p'range);
    totaltimer : timing_widget;
    points : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);
    ilili : Standard_Integer_VecVecs.Link_to_VecVec;
    flili : Standard_Floating_VecVecs.Link_to_VecVec;
    mixpts,ilifpts : Arrays_of_Integer_Vector_Lists.Link_to_Array_of_Lists;
    fpts,flifpts : Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;

  begin
    if vrblvl > 0 then
      put("-> in drivers_for_static_lifting.");
      put_line("Driver_for_Mixed_Volume_Computation 2 ...");
    end if;
    new_line; put_line(welcome);
    tstart(totaltimer);
    points := Create(p);
    Data_Management(file,mix,imixsub,fmixsub,compmisu,compmix,fltlif,vrblvl-1);
    Prepare_Supports(file,n,p,mix,points,mixpts,ilifpts,fpts,flifpts,
                     ilili,flili,fltlif,stlb,compmisu,compmix,sp,vrblvl-1);
    if fltlif then
      Polyhedral_Homotopies_with_Float_Lifting
        (file,nt,compmisu,byebye,n,sp,stlb,mix,fpts.all,mixpts.all,flili,
         flifpts.all,fmixsub,q,qsols,qsols0,gft,solsft,
         tosolve,ranstart,contrep,mv,smv,tmv,vrblvl-1);
    else
      Polyhedral_Homotopies_with_Integer_Lifting
       -- (file,compmisu,byebye,n,sp,mix.all,points,mixpts.all,ilifpts.all,
        (file,compmisu,byebye,n,sp,mix.all,points,ilifpts.all,
         imixsub,q,qsols,gft,solsft,tosolve,ranstart,contrep,mv,vrblvl-1);
     end if;
    if tosolve
     then Write_Results(file,n,ranstart,gft,solsft,q,qsols,qsols0);
    end if;
    tstop(totaltimer);
    new_line(file);
    print_times(file,totaltimer,"All Computations");
  end Driver_for_Mixed_Volume_Computation;

end Drivers_for_Static_Lifting;
