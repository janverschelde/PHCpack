with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Integer_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions;

package Black_Mixed_Volume_Computations is

-- DESCRIPTION :
--   Offers black-box routines to compute mixed volumes and to solve
--   a random coefficient start system.

  procedure Black_Box_Mixed_Volume_Computation
               ( p : in out Poly_Sys; mix : out Link_to_Vector;
                 lifsup : out 
                   Arrays_of_Integer_Vector_Lists.Link_to_Array_of_Lists;
                 mixsub : out Integer_Mixed_Subdivisions.Mixed_Subdivision;
                 mv : out natural32; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Selects the appropriate algorithm to compute the mixed volume,
  --   using integer-valued lifting functions.

  -- ON ENTRY :
  --   p         polynomial system.

  -- ON RETURN :
  --   p         permuted if type of mixture is not fully mixed;
  --   mix       type of mixture;
  --   lifsup    lifted supports of the system;
  --   mixsub    a regular mixed-cell configuration;
  --   mv        mixed volume.

  procedure Make_Induced_Permutation
              ( sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                iprm : out Standard_Integer_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Auxiliary procedure to make the induced permutation
  --   as needed for semi-mixed polynomial systems.

  -- ON ENTRY :
  --   sup      the supports of a polynomial system;
  --   mix      type of mixture;
  --   mcc      a regular mixed-cell configuration.

  -- ON RETURN
  --   iprm     the induced permutation from the mixed-cell configuration.

  procedure Make_Induced_Permutation
              ( sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                stlb : in double_float;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                iprm : out Standard_Integer_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Auxiliary procedure to make the induced permutation
  --   as needed for semi-mixed polynomial systems.

  -- ON ENTRY :
  --   sup      the supports of a polynomial system;
  --   stlb     the stable lifting bound;
  --   mix      type of mixture;
  --   mcc      a regular mixed-cell configuration.

  -- ON RETURN
  --   iprm     the induced permutation from the mixed-cell configuration.

  procedure Make_Induced_Permutation
              ( p : in Poly_Sys;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                iprm : out Standard_Integer_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Auxiliary procedure to make the induced permutation
  --   as needed for semi-mixed polynomial systems.

  -- ON ENTRY :
  --   p        a polynomial system;
  --   mix      type of mixture;
  --   mcc      a regular mixed-cell configuration.

  -- ON RETURN
  --   iprm     the induced permutation from the mixed-cell configuration.

  procedure Black_Box_Mixed_Volume_Computation
               ( p : in out Poly_Sys;
                 mix,perm,iprm : out Link_to_Vector;
                 lifsup : out 
                   Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
                 mixsub : out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 mv : out natural32; verbose : in integer32 := 0 );
  procedure Black_Box_Mixed_Volume_Computation
               ( p : in out Laur_Sys;
                 mix,perm,iprm : out Link_to_Vector;
                 lifsup : out 
                   Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
                 mixsub : out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 mv : out natural32; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Selects the appropriate algorithm to compute the mixed volume,
  --   using floating-valued random lifting functions.

  -- ON ENTRY :
  --   p         a (Laurent) polynomial system.

  -- ON RETURN :
  --   p         permuted if type of mixture is not fully mixed;
  --   mix       type of mixture;
  --   perm      permutation of equations in p for semimixed systems;
  --   iprm      induced permutation on the supports;
  --   lifsup    lifted supports of the system;
  --   mixsub    a regular mixed-cell configuration;
  --   mv        mixed volume.

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
                 verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Selects the appropriate algorithm to compute the stable mixed volume,
  --   using floating-valued random lifting functions.

  -- ON ENTRY :
  --   p         polynomial system;
  --   verbose   the verbose level.

  -- ON RETURN :
  --   p         permuted if type of mixture is not fully mixed;
  --   stlb      lifting used for the articial origin;
  --   mix       type of mixture;
  --   perm      permutation of equations in p for semimixed systems;
  --   iprm      induced permutation of the supports;
  --   lifsup    lifted supports of the system;
  --   mixsub    a regular mixed-cell configuration;
  --   orgmcc    original mixed cells, without artifical origin;
  --   stbmcc    extra stable mixed cells;
  --   mv        mixed volume;
  --   smv       stable mixed volume.
  --   tmv       total mixed volume;
  --   orgcnt    number of original mixed cells;
  --   stbcnt    number of extra stable mixed cells.

end Black_Mixed_Volume_Computations;
