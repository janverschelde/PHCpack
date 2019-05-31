with text_io;                            use text_io;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors; 
with Supports_of_Polynomial_Systems;
with Standard_Integer32_Vertices;        use Standard_Integer32_Vertices;
with Floating_Integer_Convertors;
with Floating_Lifting_Functions;
with Floating_Lifting_Utilities;
with Induced_Permutations;
with Cayley_Trick;                       use Cayley_Trick;
with Standard_Integer32_Triangulations;  use Standard_Integer32_Triangulations;
with Standard_Dynamic32_Triangulations;
with Triangulations_and_Subdivisions;    use Triangulations_and_Subdivisions;
with Flatten_Mixed_Subdivisions;         use Flatten_Mixed_Subdivisions;
with Mixed_Volume_Computation;           use Mixed_Volume_Computation;
with Drivers_for_Static_Lifting;
with Drivers_for_MixedVol_Algorithm;

--with Standard_Integer_Numbers_io;
-- use Standard_Integer_Numbers_io;
--with Standard_Integer_Vectors_io;
-- use Standard_Integer_Vectors_io;

package body Black_Mixed_Volume_Computations is

  procedure Black_Box_Mixed_Volume_Computation
               ( p : in out Poly_Sys; mix : out Link_to_Vector;
                 lifsup : out
                   Arrays_of_Integer_Vector_Lists.Link_to_Array_of_Lists;
                 mixsub : out Integer_Mixed_Subdivisions.Mixed_Subdivision;
                 mv : out natural32; verbose : in integer32 := 0 ) is

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
      if Is_Null(tmixsub) then
        cells := Deep_Create(n,t);
      else
        cells := Non_Flat_Deep_Create(n,t);
        Construct(Head_Of(tmixsub),cells);
      end if;
      Flatten(cells);
      tmixsub := cells;
    end Collect_Flattening;
    procedure C_Dynamic_Lifting is
      new Standard_Dynamic32_Triangulations.Dynamic_Lifting_with_Flat
            (Collect_Flattening);

  begin
    if verbose > 0 then
      put("-> in black_mixed_volume_computation.");
      put_line("Black_Box_Mixed_Volume_Computation 1 ...");
    end if;
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

  procedure Make_Induced_Permutation
              ( sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                iprm : out Standard_Integer_Vectors.Link_to_Vector ) is

    fs : Arrays_of_Floating_Vector_Lists.Array_of_Lists(sup'range)
       := Floating_Integer_Convertors.Convert(sup);
    ls : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
       := Floating_Lifting_Utilities.Lifted_Supports(mix'last,mcc);
    ip : Standard_Integer_Vectors.Vector(fs'range);

  begin
    ip := Induced_Permutations.Permutation(fs,ls,mix);
    iprm := new Standard_Integer_Vectors.Vector'(ip);
    Arrays_of_Floating_Vector_Lists.Deep_Clear(fs);
    Arrays_of_Floating_Vector_Lists.Deep_Clear(ls);
  end Make_Induced_Permutation;

  procedure Make_Induced_Permutation
              ( sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                stlb : in double_float;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                iprm : out Standard_Integer_Vectors.Link_to_Vector ) is

    fs : Arrays_of_Floating_Vector_Lists.Array_of_Lists(sup'range)
       := Floating_Integer_Convertors.Convert(sup);
    ls : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
       := Floating_Lifting_Utilities.Lifted_Supports(mix'last,mcc);
    ip : Standard_Integer_Vectors.Vector(fs'range);

  begin
    Induced_Permutations.Remove_Artificial_Origin(ls,stlb);
    ip := Induced_Permutations.Permutation(fs,ls,mix);
    iprm := new Standard_Integer_Vectors.Vector'(ip);
    Arrays_of_Floating_Vector_Lists.Deep_Clear(fs);
    Arrays_of_Floating_Vector_Lists.Deep_Clear(ls);
  end Make_Induced_Permutation;

  procedure Make_Induced_Permutation
              ( p : in Poly_Sys;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                iprm : out Standard_Integer_Vectors.Link_to_Vector ) is

    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range)
        := Supports_of_Polynomial_Systems.Create(p);

  begin
    Make_Induced_Permutation(sup,mix,mcc,iprm);
    Arrays_of_Integer_Vector_Lists.Deep_Clear(sup);
  end Make_Induced_Permutation;

  procedure Make_Induced_Permutation
              ( sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                lifsup : out 
                  Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
                iprm : out Standard_Integer_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Auxiliary procedure to make the induced permutation
  --   as needed for semi-mixed polynomial systems.

  -- ON ENTRY :
  --   sup      the supports of the polynomials system;
  --   mix      type of mixture;
  --   mcc      a regular mixed-cell configuration.

  -- ON RETURN
  --   lifsup   the lifted supports;
  --   iprm     the induced permutation from the mixed-cell configuration.

    fs : Arrays_of_Floating_Vector_Lists.Array_of_Lists(sup'range)
       := Floating_Integer_Convertors.Convert(sup);
    ls : constant Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
       := Floating_Lifting_Utilities.Lifted_Supports(mix'last,mcc);
    ip : Standard_Integer_Vectors.Vector(fs'range);

  begin
    ip := Induced_Permutations.Permutation(fs,ls,mix);
    iprm := new Standard_Integer_Vectors.Vector'(ip);
    Arrays_of_Floating_Vector_Lists.Deep_Clear(fs);
    lifsup := new Arrays_of_Floating_Vector_Lists.Array_of_Lists'(ls);
  end Make_Induced_Permutation;

  procedure Make_Induced_Permutation
              ( sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                stlb : in double_float;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                lifsup : out 
                  Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
                iprm : out Standard_Integer_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Auxiliary procedure to make the induced permutation
  --   as needed for semi-mixed polynomial systems.

  -- ON ENTRY :
  --   sup      the supports of the polynomials system;
  --   stlb     stable lifting bound;
  --   mix      type of mixture;
  --   mcc      a regular mixed-cell configuration.

  -- ON RETURN
  --   lifsup   the lifted supports;
  --   iprm     the induced permutation from the mixed-cell configuration.

    fs : Arrays_of_Floating_Vector_Lists.Array_of_Lists(sup'range)
       := Floating_Integer_Convertors.Convert(sup);
    ls : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
       := Floating_Lifting_Utilities.Lifted_Supports(mix'last,mcc);
    ip : Standard_Integer_Vectors.Vector(fs'range);

  begin
    Induced_Permutations.Remove_Artificial_Origin(ls,stlb);
    ip := Induced_Permutations.Permutation(fs,ls,mix);
    iprm := new Standard_Integer_Vectors.Vector'(ip);
    Arrays_of_Floating_Vector_Lists.Deep_Clear(fs);
    lifsup := new Arrays_of_Floating_Vector_Lists.Array_of_Lists'(ls);
  end Make_Induced_Permutation;

  procedure Black_Box_Mixed_Volume_Computation
               ( p : in out Poly_Sys;
                 mix,perm,iprm : out Link_to_Vector;
                 lifsup : out 
                   Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
                 mixsub : out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 mv : out natural32; verbose : in integer32 := 0 ) is

   -- n : constant integer32 := p'last;
    sup : constant Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range)
        := Supports_of_Polynomial_Systems.Create(p);
    r : integer32;

  begin
    if verbose > 0 then
      put("-> in black_mixed_volume_computation.");
      put_line("Black_Box_Mixed_Volume_Computation 2 ...");
    end if;
    Drivers_for_MixedVol_Algorithm.Mixed_Volume_Computation
      (p'last,sup,0.0,r,mix,perm,mixsub,mv);
   -- declare
   --   fs : Arrays_of_Floating_Vector_Lists.Array_of_Lists(sup'range)
   --      := Floating_Integer_Convertors.Convert(sup);
   --   ls : constant Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
   --      := Floating_Lifting_Utilities.Lifted_Supports(mix'last,mixsub);
   --   ip : Standard_Integer_Vectors.Vector(fs'range);
   -- begin
   --   ip := Induced_Permutations.Permutation(fs,ls,mix.all);
   --   iprm := new Standard_Integer_Vectors.Vector'(ip);
     -- put("induced permutation : "); put(ip); new_line;
   --   Induced_Permutations.Permute(ip,p);
   --   Arrays_of_Floating_Vector_Lists.Deep_Clear(fs);
   --   lifsup := new Arrays_of_Floating_Vector_Lists.Array_of_Lists'(ls);
   -- end;
    Make_Induced_Permutation(sup,mix(1..r),mixsub,lifsup,iprm);
    Induced_Permutations.Permute(iprm.all,p);
  end Black_Box_Mixed_Volume_Computation;

  procedure Black_Box_Mixed_Volume_Computation
               ( p : in out Laur_Sys;
                 mix,perm,iprm : out Link_to_Vector;
                 lifsup : out 
                   Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
                 mixsub : out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 mv : out natural32; verbose : in integer32 := 0 ) is

   -- n : constant integer32 := p'last;
    sup : constant Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range)
        := Supports_of_Polynomial_Systems.Create(p);
    r : integer32;

  begin
    if verbose > 0 then
      put("-> in black_mixed_volume_computation.");
      put_line("Black_Box_Mixed_Volume_Computation 3 ...");
    end if;
    Drivers_for_MixedVol_Algorithm.Mixed_Volume_Computation
      (p'last,sup,0.0,r,mix,perm,mixsub,mv);
   -- declare
   --   fs : Arrays_of_Floating_Vector_Lists.Array_of_Lists(sup'range)
   --      := Floating_Integer_Convertors.Convert(sup);
   --   ls : constant Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
   --      := Floating_Lifting_Utilities.Lifted_Supports(mix'last,mixsub);
   --   ip : Standard_Integer_Vectors.Vector(fs'range);
   -- begin
   --   ip := Induced_Permutations.Permutation(fs,ls,mix.all);
   --   iprm := new Standard_Integer_Vectors.Vector'(ip);
     -- put("induced permutation : "); put(ip); new_line;
   --   Induced_Permutations.Permute(ip,p);
   --   Arrays_of_Floating_Vector_Lists.Deep_Clear(fs);
   --   lifsup := new Arrays_of_Floating_Vector_Lists.Array_of_Lists'(ls);
   -- end;
    Make_Induced_Permutation(sup,mix(1..r),mixsub,lifsup,iprm);
    Induced_Permutations.Permute(iprm.all,p);
  end Black_Box_Mixed_Volume_Computation;

  procedure Black_Box_Mixed_Volume_Computation
               ( p : in out Poly_Sys;
                 mix,perm,iprm : out Link_to_Vector;
                 stlb : out double_float;
                 lifsup : out 
                   Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
                 mixsub : out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 orgmcc : out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 stbmcc : out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 mv,smv,tmv,orgcnt,stbcnt : out natural32;
                 verbose : in integer32 := 0 ) is

    n : constant integer32 := p'last;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range)
        := Supports_of_Polynomial_Systems.Create(p);
    r : integer32;

  begin
    if verbose > 0 then
      put("-> in black_mixed_volume_computation.");
      put_line("Black_Box_Mixed_Volume_Computation 4 ...");
    end if;
    stlb := Floating_Lifting_Functions.Lifting_Bound(p);
    Drivers_for_MixedVol_Algorithm.Mixed_Volume_Computation
      (p'last,sup,stlb,r,mix,perm,mixsub,mv);
   -- declare
   --   fs : Arrays_of_Floating_Vector_Lists.Array_of_Lists(sup'range)
   --      := Floating_Integer_Convertors.Convert(sup);
   --   ls : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
   --      := Floating_Lifting_Utilities.Lifted_Supports(mix'last,mixsub);
   --   ip : Standard_Integer_Vectors.Vector(fs'range);
   -- begin
   --   Induced_Permutations.Remove_Artificial_Origin(ls,stlb);
   --   ip := Induced_Permutations.Permutation(fs,ls,mix.all);
   --   iprm := new Standard_Integer_Vectors.Vector'(ip);
     -- put("ip'first = "); put(ip'first,1); new_line;
     -- put("induced permutation : "); put(ip); new_line;
   --   Induced_Permutations.Permute(ip,p);
   --   Arrays_of_Floating_Vector_Lists.Deep_Clear(fs);
   --   lifsup := new Arrays_of_Floating_Vector_Lists.Array_of_Lists'(ls);
   -- end;
    Make_Induced_Permutation(sup,stlb,mix(1..r),mixsub,lifsup,iprm);
    Induced_Permutations.Permute(iprm.all,p);
    Drivers_for_Static_Lifting.Floating_Volume_Computation
      (n,stlb,mix(1..r),mixsub,mv,smv,tmv);
    Floating_Mixed_Subdivisions.Split_Original_Cells
      (mixsub,stlb,orgmcc,stbmcc,orgcnt,stbcnt);
  end Black_Box_Mixed_Volume_Computation;

end Black_Mixed_Volume_Computations;
