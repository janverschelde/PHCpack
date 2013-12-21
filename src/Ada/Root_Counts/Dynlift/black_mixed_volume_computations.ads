with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
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
                 mv : out natural32 );

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

  procedure Black_Box_Mixed_Volume_Computation
               ( p : in out Poly_Sys; mix,perm : out Link_to_Vector;
                 lifsup : out 
                   Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
                 mixsub : out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 mv : out natural32 );
  procedure Black_Box_Mixed_Volume_Computation
               ( p : in out Laur_Sys; mix,perm : out Link_to_Vector;
                 lifsup : out 
                   Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
                 mixsub : out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 mv : out natural32 );

  -- DESCRIPTION :
  --   Selects the appropriate algorithm to compute the mixed volume,
  --   using floating-valued random lifting functions.

  -- ON ENTRY :
  --   p         a (Laurent) polynomial system.

  -- ON RETURN :
  --   p         permuted if type of mixture is not fully mixed;
  --   mix       type of mixture;
  --   perm      permutation of equations in p for semimixed systems;
  --   lifsup    lifted supports of the system;
  --   mixsub    a regular mixed-cell configuration;
  --   mv        mixed volume.

  procedure Black_Box_Mixed_Volume_Computation
               ( p : in out Poly_Sys; mix,perm : out Link_to_Vector;
                 stlb : out double_float;
                 lifsup : out 
                   Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
                 mixsub : out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 orgmcc : out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 stbmcc : out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 mv,smv,tmv,orgcnt,stbcnt : out natural32 );

  -- DESCRIPTION :
  --   Selects the appropriate algorithm to compute the stable mixed volume,
  --   using floating-valued random lifting functions.

  -- ON ENTRY :
  --   p         polynomial system.

  -- ON RETURN :
  --   p         permuted if type of mixture is not fully mixed;
  --   stlb      lifting used for the articial origin;
  --   mix       type of mixture;
  --   perm      permutation of equations in p for semimixed systems;
  --   lifsup    lifted supports of the system;
  --   mixsub    a regular mixed-cell configuration;
  --   orgmcc    original mixed cells, without artifical origin;
  --   stbmcc    extra stable mixed cells;
  --   mv        mixed volume;
  --   smv       stable mixed volume.
  --   tmv       total mixed volume;
  --   orgcnt    number of original mixed cells;
  --   stbcnt    number of extra stable mixed cells.

  procedure Black_Box_Polyhedral_Continuation
               ( p : in Poly_Sys; mix : in Vector;
                 lifsup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 mixsub : in Integer_Mixed_Subdivisions.Mixed_Subdivision;
                 q : in out Poly_Sys; qsols : in out Solution_List );
  procedure Black_Box_Polyhedral_Continuation
               ( p : in Laur_Sys; mix : in Vector;
                 lifsup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 mixsub : in Integer_Mixed_Subdivisions.Mixed_Subdivision;
                 q : in out Laur_Sys; qsols : in out Solution_List );

  -- DESCRIPTION :
  --   Creates and solves a random coefficient start system, based on 
  --   a regular mixed-cell configuration, induced by integer lifting.

  -- ON ENTRY :
  --   p         a (Laurent) polynomial system;
  --   mix       type of mixture;
  --   lifsup    lifted supports of the system;
  --   mixsub    regular mixed-cell configuration;
  --   mv        mixed volume.

  -- ON RETURN :
  --   q         random coefficient start system;
  --   qsols     solutions of q.

  procedure Black_Box_Polyhedral_Continuation
               ( nt : in integer32;
                 p : in Laur_Sys; mix,perm : in Link_to_Vector;
                 lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 mcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 q : in out Laur_Sys; qsols : in out Solution_List );

  procedure Black_Box_Polyhedral_Continuation
               ( nt : in integer32;
                 p : in Poly_Sys; mix,perm : in Link_to_Vector;
                 stlb : in double_float;
                 lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 orgmcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 stbmcc : in Floating_Mixed_Subdivisions.Mixed_Subdivision;
                 q : in out Poly_Sys;
                 qsols,qsols0 : in out Solution_List );

  -- DESCRIPTION :
  --   Creates and solves a random coefficient start system, based on 
  --   a regular mixed-cell configuration, induced by floating lifting.

  -- ON ENTRY :
  --   nt        number of tasks for the polyhedral continuation,
  --             if 0, then sequential execution;
  --   p         polynomial system;
  --   mix       type of mixture;
  --   perm      permutation of the equations in p;
  --   stlb      lifting used for the artificial origin;
  --   lifsup    lifted supports of the system;
  --   orgmcc    a regular mixed-cell configuration with only
  --             the original mixed cells, without artificial origin;
  --   stbmcc    a regular mixed-cell configuration, with only
  --             the extra stable mixed cells, with aritifical origin.

  -- ON RETURN :
  --   q         random coefficient start system;
  --   qsols     solutions of q, without zero components;
  --   qsols0    solutions of q with zero components.

end Black_Mixed_Volume_Computations;
