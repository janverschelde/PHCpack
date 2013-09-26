with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Poly_Randomizers; 
with Standard_Complex_Laur_Randomizers;  
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Functions;    use Standard_Complex_Laur_Functions;
with Standard_Complex_Laur_SysFun;       use Standard_Complex_Laur_SysFun;
with Standard_Complex_Laur_Jacomats;     use Standard_Complex_Laur_Jacomats;
with Standard_Poly_Laur_Convertors;      use Standard_Poly_Laur_Convertors;
with Standard_Laur_Poly_Convertors;      use Standard_Laur_Poly_Convertors;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors; 
with Continuation_Parameters;
with Supports_of_Polynomial_Systems;
with Standard_Integer_Vertices;          use Standard_Integer_Vertices;
with Random_Coefficient_Systems;
with Integer_Lifting_Utilities;          use Integer_Lifting_Utilities;
with Floating_Integer_Convertors;
with Floating_Lifting_Functions;
with Floating_Lifting_Utilities;
with Induced_Permutations;
with Cayley_Trick;                       use Cayley_Trick;
with Standard_Integer_Triangulations;    use Standard_Integer_Triangulations;
with Standard_Dynamic_Triangulations;    use Standard_Dynamic_Triangulations;
with Triangulations_and_Subdivisions;    use Triangulations_and_Subdivisions;
with Flatten_Mixed_Subdivisions;         use Flatten_Mixed_Subdivisions;
with Mixed_Volume_Computation;           use Mixed_Volume_Computation;
with Exponent_Vectors;                   use Exponent_Vectors;
with Integer_Polyhedral_Continuation;    use Integer_Polyhedral_Continuation;
with Floating_Polyhedral_Continuation;   use Floating_Polyhedral_Continuation;
with Stable_Polyhedral_Continuation;     use Stable_Polyhedral_Continuation;
with Drivers_for_Static_Lifting;
with Drivers_for_MixedVol_Algorithm;     use Drivers_for_MixedVol_Algorithm;
with Multitasking_Polyhedral_Trackers;   use Multitasking_Polyhedral_Trackers;

package body Black_Mixed_Volume_Computations is

  procedure Black_Box_Mixed_Volume_Computation
               ( p : in out Poly_Sys; mix : out Link_to_Vector;
                 lifsup : out
                   Arrays_of_Integer_Vector_Lists.Link_to_Array_of_Lists;
                 mixsub : out Integer_Mixed_Subdivisions.Mixed_Subdivision;
                 mv : out natural32 ) is

    use Arrays_of_Integer_Vector_Lists;
    use Integer_Mixed_Subdivisions;

    n : constant integer32 := p'length;
    supports : constant Array_of_Lists(p'range)
             := Supports_of_Polynomial_Systems.Create(p);
    verpts : Array_of_Lists(p'range);
    tmix,perms : Standard_Integer_Vectors.Link_to_Vector;
    tmixsub : Mixed_Subdivision;

    procedure Collect_Flattening ( t : in Triangulation; l : List ) is

    -- DESCRIPTION :
    --   Updates the subdivision tmixsub with the flattened cells.
    --   The triangulation on entry contains the whole triangulation,
    --   not just the new cells.

      cells : Mixed_Subdivision;

    begin
      if Is_Null(tmixsub)
       then cells := Deep_Create(n,t);
       else cells := Non_Flat_Deep_Create(n,t);
            Construct(Head_Of(tmixsub),cells);
      end if;
      Flatten(cells);
      tmixsub := cells;
    end Collect_Flattening;
    procedure C_Dynamic_Lifting is
      new Standard_Dynamic_Triangulations.Dynamic_Lifting_with_Flat
            (Collect_Flattening);

  begin
    for i in supports'range loop
      verpts(i) := Vertex_Points(supports(i));
    end loop;
    Compute_Mixture(verpts,tmix,perms);
    p := Permute(p,perms);
    declare
      pts,lifted : Array_of_Lists(tmix'range);
      last : List;
      t : Triangulation;
      nt : natural32;
      lastcells : Mixed_Subdivision;
    begin
      if tmix'length = 1 then
        C_Dynamic_Lifting(verpts(1),false,false,0,lifted(1),last,t);
        if Is_Null(tmixsub) then
          tmixsub := Deep_Create(n,t);
        else
          lastcells := Non_Flat_Deep_Create(n,t);
          Construct(Head_Of(tmixsub),lastcells);
          tmixsub := lastcells;
        end if;
        Clear(t);
        Mixed_Volume(n,tmix.all,tmixsub,mv);
      elsif tmix'length <= n/2 then
        pts := Typed_Lists(tmix.all,verpts);
        Dynamic_Cayley(n,tmix.all,pts,false,false,0,lifted,tmixsub,nt);
        Mixed_Volume(n,tmix.all,tmixsub,mv);
      else
        Mixed_Volume(n,tmix.all,verpts,lifted,tmixsub,mv);
      end if;
      lifsup := new Array_of_Lists'(lifted);
    end;
    mix := tmix; mixsub := tmixsub;
  end Black_Box_Mixed_Volume_Computation;

  procedure Black_Box_Mixed_Volume_Computation
               ( p : in out Poly_Sys; mix,perm : out Link_to_Vector;
                 lifsup : out 
                   Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
                 mixsub : out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 mv : out natural32 ) is

   -- n : constant integer32 := p'last;
    sup : constant Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range)
        := Supports_of_Polynomial_Systems.Create(p);
    r : integer32;

  begin
    Drivers_for_MixedVol_Algorithm.Mixed_Volume_Computation
      (p'last,sup,0.0,r,mix,perm,mixsub,mv);
    declare
      fs : Arrays_of_Floating_Vector_Lists.Array_of_Lists(sup'range)
         := Floating_Integer_Convertors.Convert(sup);
      ls : constant Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
         := Floating_Lifting_Utilities.Lifted_Supports(mix'last,mixsub);
      ip : Standard_Integer_Vectors.Vector(fs'range);
    begin
      ip := Induced_Permutations.Permutation(fs,ls,mix.all);
      Induced_Permutations.Permute(ip,p);
      Arrays_of_Floating_Vector_Lists.Deep_Clear(fs);
      lifsup := new Arrays_of_Floating_Vector_Lists.Array_of_Lists'(ls);
    end;
  end Black_Box_Mixed_Volume_Computation;

  procedure Black_Box_Mixed_Volume_Computation
               ( p : in out Laur_Sys; mix,perm : out Link_to_Vector;
                 lifsup : out 
                   Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
                 mixsub : out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 mv : out natural32 ) is

   -- n : constant integer32 := p'last;
    sup : constant Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range)
        := Supports_of_Polynomial_Systems.Create(p);
    r : integer32;

  begin
    Drivers_for_MixedVol_Algorithm.Mixed_Volume_Computation
      (p'last,sup,0.0,r,mix,perm,mixsub,mv);
    declare
      fs : Arrays_of_Floating_Vector_Lists.Array_of_Lists(sup'range)
         := Floating_Integer_Convertors.Convert(sup);
      ls : constant Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
         := Floating_Lifting_Utilities.Lifted_Supports(mix'last,mixsub);
      ip : Standard_Integer_Vectors.Vector(fs'range);
    begin
      ip := Induced_Permutations.Permutation(fs,ls,mix.all);
      Induced_Permutations.Permute(ip,p);
      Arrays_of_Floating_Vector_Lists.Deep_Clear(fs);
      lifsup := new Arrays_of_Floating_Vector_Lists.Array_of_Lists'(ls);
    end;
  end Black_Box_Mixed_Volume_Computation;

  procedure Black_Box_Mixed_Volume_Computation
               ( p : in out Poly_Sys; mix,perm : out Link_to_Vector;
                 stlb : out double_float;
                 lifsup : out 
                   Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
                 mixsub : out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 orgmcc : out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 stbmcc : out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 mv,smv,tmv,orgcnt,stbcnt : out natural32 ) is

    n : constant integer32 := p'last;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range)
        := Supports_of_Polynomial_Systems.Create(p);
    r : integer32;

  begin
    stlb := Floating_Lifting_Functions.Lifting_Bound(p);
    Drivers_for_MixedVol_Algorithm.Mixed_Volume_Computation
      (p'last,sup,stlb,r,mix,perm,mixsub,mv);
    declare
      fs : Arrays_of_Floating_Vector_Lists.Array_of_Lists(sup'range)
         := Floating_Integer_Convertors.Convert(sup);
      ls : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
         := Floating_Lifting_Utilities.Lifted_Supports(mix'last,mixsub);
      ip : Standard_Integer_Vectors.Vector(fs'range);
    begin
      Induced_Permutations.Remove_Artificial_Origin(ls,stlb);
      ip := Induced_Permutations.Permutation(fs,ls,mix.all);
      Induced_Permutations.Permute(ip,p);
      Arrays_of_Floating_Vector_Lists.Deep_Clear(fs);
      lifsup := new Arrays_of_Floating_Vector_Lists.Array_of_Lists'(ls);
    end;
    Drivers_for_Static_Lifting.Floating_Volume_Computation
      (n,stlb,mix.all,mixsub,mv,smv,tmv);
    Floating_Mixed_Subdivisions.Split_Original_Cells
      (mixsub,stlb,orgmcc,stbmcc,orgcnt,stbcnt);
  end Black_Box_Mixed_Volume_Computation;

  procedure Black_Box_Polyhedral_Continuation
               ( p : in Poly_Sys; mix : in Vector;
                 lifsup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 mixsub : in Integer_Mixed_Subdivisions.Mixed_Subdivision;
                 q : in out Poly_Sys; qsols : in out Solution_List ) is

    n : constant integer32 := p'length;
    lq,llq : Laur_Sys(p'range);
    h : Eval_Coeff_Laur_Sys(q'range);
    c : Standard_Complex_VecVecs.VecVec(h'range);
    e : Exponent_Vectors_Array(h'range);
    j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
    m : Mult_Factors(j'range(1),j'range(2));

  begin
    q := Standard_Complex_Poly_Randomizers.Complex_Randomize1(p);
    lq := Polynomial_to_Laurent_System(q);
    llq := Perform_Lifting(n,mix,lifsup,lq);
    Clear(lq); Clear(q);
    lq := Eval(llq,Create(1.0),n+1);
    q := Laurent_to_Polynomial_System(lq);
    Continuation_Parameters.Tune(0);
   -- Mixed_Solve(llq,mix,mixsub,qsols);   too expensive !!!!
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
    Mixed_Solve(llq,lifsup,h,c,e,j,m,mix,mixsub,qsols);
    Set_Continuation_Parameter(qsols,Create(0.0));
    Clear(lq); Clear(llq);
    Clear(h); Clear(j); Clear(m);
    Standard_Complex_VecVecs.Clear(c);
  end Black_Box_Polyhedral_Continuation;

  procedure Black_Box_Polyhedral_Continuation
               ( p : in Laur_Sys; mix : in Vector;
                 lifsup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 mixsub : in Integer_Mixed_Subdivisions.Mixed_Subdivision;
                 q : in out Laur_Sys; qsols : in out Solution_List ) is

    n : constant integer32 := p'length;
    lq,llq : Laur_Sys(p'range);
    h : Eval_Coeff_Laur_Sys(q'range);
    c : Standard_Complex_VecVecs.VecVec(h'range);
    e : Exponent_Vectors_Array(h'range);
    j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
    m : Mult_Factors(j'range(1),j'range(2));

  begin
    lq := Standard_Complex_Laur_Randomizers.Complex_Randomize1(p);
    llq := Perform_Lifting(n,mix,lifsup,lq);
    Clear(lq); Clear(q);
    q := Eval(llq,Create(1.0),n+1);
    Continuation_Parameters.Tune(0);
   -- Mixed_Solve(llq,mix,mixsub,qsols);   too expensive !!!!
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
    Mixed_Solve(llq,lifsup,h,c,e,j,m,mix,mixsub,qsols);
    Set_Continuation_Parameter(qsols,Create(0.0));
    Clear(llq); Clear(h); Clear(j); Clear(m);
    Standard_Complex_VecVecs.Clear(c);
  end Black_Box_Polyhedral_Continuation;

  procedure Black_Box_Polyhedral_Continuation
               ( nt : in integer32;
                 p : in Laur_Sys; mix,perm : in Link_to_Vector;
                 lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 mcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 q : in out Laur_Sys; qsols : in out Solution_List ) is

    hq : Eval_Coeff_Laur_Sys(p'range);
    expvec : Exponent_Vectors_Array(p'range);
    coeffv : Standard_Complex_VecVecs.VecVec(p'range);
    jacmat : Eval_Coeff_Jaco_Mat(p'range,p'first..p'last+1);
    mulfac : Mult_Factors(jacmat'range(1),jacmat'range(2));

  begin
    q := Random_Coefficient_Systems.Create(natural32(p'last),mix.all,lifsup);
    hq := Create(q);
    expvec := Create(q);
    for i in q'range loop
      declare
        c : constant Standard_Complex_Vectors.Vector := Coeff(q(i));
      begin
        coeffv(i) := new Standard_Complex_Vectors.Vector(c'range);
        for j in c'range loop
          coeffv(i)(j) := c(j);
        end loop;
      end;
    end loop;
    Create(q,jacmat,mulfac);
    if nt = 0 then
      Mixed_Solve(q,lifsup,hq,coeffv,expvec,jacmat,mulfac,
                  mix.all,mcc,qsols);
    else
      Silent_Multitasking_Path_Tracker
        (q,nt,p'last,mix'last,mix.all,lifsup,mcc,
         hq,coeffv,expvec,jacmat,mulfac,qsols);
    end if;
    Set_Continuation_Parameter(qsols,Create(0.0));
    Clear(hq); Clear(jacmat); Clear(mulfac);
    Standard_Complex_VecVecs.Clear(coeffv);
  end Black_Box_Polyhedral_Continuation;

  procedure Black_Box_Polyhedral_Continuation
               ( nt : in integer32;
                 p : in Poly_Sys; mix,perm : in Link_to_Vector;
                 stlb : in double_float;
                 lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 orgmcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 stbmcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 q : in out Poly_Sys;
                 qsols,qsols0 : in out Solution_List ) is

    lq : Laur_Sys(p'range);
    hq : Eval_Coeff_Laur_Sys(p'range);
    expvec : Exponent_Vectors_Array(p'range);
    coeffv : Standard_Complex_VecVecs.VecVec(p'range);
    jacmat : Eval_Coeff_Jaco_Mat(p'range,p'first..p'last+1);
    mulfac : Mult_Factors(jacmat'range(1),jacmat'range(2));

  begin
    q := Random_Coefficient_Systems.Create(natural32(p'last),mix.all,lifsup);
    lq := Polynomial_to_Laurent_System(q);
    hq := Create(lq);
    expvec := Create(q);
    for i in q'range loop
      declare
        c : constant Standard_Complex_Vectors.Vector := Coeff(lq(i));
      begin
        coeffv(i) := new Standard_Complex_Vectors.Vector(c'range);
        for j in c'range loop
          coeffv(i)(j) := c(j);
        end loop;
      end;
    end loop;
    Create(lq,jacmat,mulfac);
    if nt = 0 then
      Mixed_Solve(lq,lifsup,hq,coeffv,expvec,jacmat,mulfac,
                  mix.all,orgmcc,qsols);
    else
      Silent_Multitasking_Path_Tracker
        (lq,nt,p'last,mix'last,mix.all,lifsup,orgmcc,
         hq,coeffv,expvec,jacmat,mulfac,qsols);
    end if;
    Set_Continuation_Parameter(qsols,Create(0.0));
    if not Floating_Mixed_Subdivisions.Is_Null(stbmcc) then
      Silent_Polyhedral_Continuation(lq,stlb,mix,lifsup,stbmcc,qsols0);
     -- put_line("looking at the stable mixed cells ...");
     -- Reporting_Polyhedral_Continuation
     --   (standard_output,lq,stlb,mix,lifsup,stbmcc,qsols0);
      Set_Continuation_Parameter(qsols0,Create(0.0));
     -- put("Length_Of(qsols0) = "); put(Length_Of(qsols0),1); new_line;
   -- else
   --   put_line("no stable mixed cells");
    end if;
    Clear(lq); Clear(hq); Clear(jacmat); Clear(mulfac);
    Standard_Complex_VecVecs.Clear(coeffv);
  end Black_Box_Polyhedral_Continuation;

end Black_Mixed_Volume_Computations;
