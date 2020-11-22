with text_io;                            use text_io;
with String_Splitters;                   use String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Natural_Numbers;           use Multprec_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Arrays_of_Floating_Vector_Lists;    use Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Partitions_of_Sets_Of_Unknowns;     use Partitions_of_Sets_of_Unknowns;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;

package Black_Box_Root_Counters is

-- DESCRIPTION :
--   This package provides an interface to the various root counting
--   capabilities in PHCpack, called by phc -b.
--   Three different levels of precision are supported:
--   standard double, double double, and quad double precision.

  procedure Compute_Total_Degree
               ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 tdeg : out natural64; mptdeg : out Natural_Number );
  procedure Compute_Total_Degree
               ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 tdeg : out natural64; mptdeg : out Natural_Number );
  procedure Compute_Total_Degree
               ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 tdeg : out natural64; mptdeg : out Natural_Number );

  -- DESCRIPTION :
  --   Computes the total degree of the system p.
  --   The output is in tdeg, except when the total degree is huge,
  --   larger than a 64-bit natural number.
  --   In that case, tdeg equals zero and mptdeg equals the total degree.

  function Set_Structure_Bound
               ( p : Standard_Complex_Poly_Systems.Poly_Sys ) return natural64;
  function Set_Structure_Bound
               ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys ) return natural64;
  function Set_Structure_Bound
               ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys ) return natural64;

  -- DESCRIPTION :
  --   Computes a set structure Bezout bound for the system p,
  --   via the construction of a random linear-product start system.

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
                 rocotime : out duration; verbose : in integer32 := 0 );
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
                 rocotime : out duration; verbose : in integer32 := 0 );
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
                 rocotime : out duration; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes four different root counts for the system p.
  --   If the flag "deg" is on, then the output parameter "mivo" is
  --   assigned to take the value of the total degree.

  -- ON ENTRY :
  --   p         a polynomial system.
  --   deg       if only degrees are used in the root count;
  --   verbose   the verbose level.

  -- ON RETURN :
  --   tode      total degree;
  --   mptode    multiprecision version of total degree (if overflow);
  --   mhbz      m-homogeneous Bezout number;
  --   setb      bound based on set structure;
  --   mivo      mixed volume;
  --   stmv      stable mixed volume;
  --   zz        partition used to compute mhbz;
  --   nz        number of sets in partition zz;
  --   lifsup    lifted supports of the system;
  --   stlb      lifting of the artificial origin;
  --   mix       type of mixture;
  --   perm      permutation of the equations in p;
  --   iprm      induced permutation of the blackbox mixed volume calculator;
  --   orgmcc    regular mixed-cell configuration to compute mivo;
  --   stbmcc    extra stable mixed cells that contribute to stmv;
  --   rocotime  is the time it took to compute the root count.

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
                 rocotime : out duration; verbose : in integer32 := 0 );
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
                 rocotime : out duration; verbose : in integer32 := 0 );
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
                 rocotime : out duration; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes four different root counts for the system p.
  --   If the flag "deg" is on, then the output parameter "mivo" is
  --   assigned to take the value of the total degree.

  -- ON ENTRY :
  --   file      output file;
  --   p         a polynomial system;
  --   deg       if only degrees are used in the root count;
  --   verbose   the verbose level.

  -- ON RETURN :
  --   tode      total degree;
  --   mptode    multiprecision version of total degree (if overflow)
  --   mhbz      m-homogeneous Bezout number;
  --   setb      bound based on set structure;
  --   mivo      mixed volume;
  --   stmv      stable mixed volume;
  --   zz        partition used to compute mhbz;
  --   nz        number of sets in partition zz;
  --   stlb      lifting for artificial origin;
  --   lifsup    lifted supports of the system;
  --   mix       type of mixture;
  --   perm      permutation of the equations in p;
  --   iprm      induced permutation of the blackbox mixed volume calculator;
  --   orgmcc    regular mixed-cell configuration to compute mivo;
  --   stbmcc    extra stable mixed cells that contribute to stmv;
  --   rocotime  is the time it took to compute the root count.

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
                 hocotime : out duration; verbose : in integer32 := 0 );
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
                 hocotime : out duration; verbose : in integer32 := 0 );
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
                 hocotime : out duration; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Constructs a start system for the minimal root count and least
  --   amount of work.

  -- ON ENTRY :
  --   p         polynomial system;
  --   d         total degree;
  --   bz        m-homogeneous Bezout number;
  --   bs        Bezout number based on set structure;
  --   mv        mixed volume;
  --   smv       stable mixed volume;
  --   z         partition that corresponds with bz;
  --   mix       type of mixture of the supports;
  --   stlb      lifting for the artificial origin;
  --   lifted    lifted supports;
  --   orgmcc    regular mixed-cell configuration to compute mv;
  --   stbmcc    extra stable mixed cells that contributed to smv;
  --   verbose   the verbose level.

  -- ON RETURN :
  --   roco      minimum(d,bz,bs,mv), provided mv /= 0;
  --   q         start system;
  --   qsols     solutions of q;
  --   qsols0    solutions of q with zero components;
  --   hocotime  is the time it took to construct the start system.

  procedure Construct_Start_System
               ( file : in file_type; nt : in integer32;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 d,bz,bs : in natural64; mv,smv : in natural32;
                 z : in Partition; mix : in Link_to_Vector;
                 stlb : in double_float; lifted : in Link_to_Array_of_Lists;
                 orgmcc,stbmcc : in Mixed_Subdivision; roco : out natural64;
                 q : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : in out Standard_Complex_Solutions.Solution_List;
                 hocotime : out duration; verbose : in integer32 := 0 );
  procedure Construct_Start_System
               ( file : in file_type; nt : in integer32;
                 p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 d,bz,bs : in natural64; mv,smv : in natural32;
                 z : in Partition; mix : in Link_to_Vector;
                 stlb : in double_float; lifted : in Link_to_Array_of_Lists;
                 orgmcc,stbmcc : in Mixed_Subdivision; roco : out natural64;
                 q : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : in out DoblDobl_Complex_Solutions.Solution_List;
                 hocotime : out duration; verbose : in integer32 := 0 );
  procedure Construct_Start_System
               ( file : in file_type; nt : in integer32;
                 p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 d,bz,bs : in natural64; mv,smv : in natural32;
                 z : in Partition; mix : in Link_to_Vector;
                 stlb : in double_float; lifted : in Link_to_Array_of_Lists;
                 orgmcc,stbmcc : in Mixed_Subdivision; roco : out natural64;
                 q : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : in out QuadDobl_Complex_Solutions.Solution_List;
                 hocotime : out duration; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Constructs a start system for the minimal root count and least
  --   amount of work, in double, double double, or quad double precision.

  -- ON ENTRY :
  --   file      output file;
  --   p         polynomial system;
  --   d         total degree;
  --   bz        m-homogeneous Bezout number;
  --   bs        Bezout number based on set structure;
  --   mv        mixed volume;
  --   smv       stable mixed volume;
  --   z         partition that corresponds with bz;
  --   mix       type of mixture of the supports;
  --   stlb      lifting of the artificial origin;
  --   orgmcc    regular mixed-cell configuration to compute mv;
  --   stbmcc    extra stable mixed cells that contribute to smv;
  --   verbose   the verbose level.

  -- ON RETURN :
  --   roco      minimum(d,bz,bs,mv), provided mv /= 0;
  --   q         start system;
  --   qsols     solutions of q;
  --   qsols0    solutions of q with zero components;
  --   hocotime  is the time it took to construct the start system.

  procedure Black_Box_Root_Counting 
               ( nt : in integer32; silent : in boolean;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );
  procedure Black_Box_Root_Counting 
               ( nt : in integer32; silent : in boolean;
                 p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out DoblDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );
  procedure Black_Box_Root_Counting 
               ( nt : in integer32; silent : in boolean;
                 p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out QuadDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Calculates four different root counts: total degree, m-homogeneous
  --   Bezout number, generalized Bezout number based on set structure,
  --   and mixed volume.  Heuristics are used for the Bezout numbers.
  --   Returns the start system with lowest root count and least amount
  --   of work, which means that linear-product start systems are preferred,
  --   when Bezout numbers equal the mixed volume.
  --   In case the mixed volume is zero because of a equation with only
  --   one monomial, the smallest Bezout bound is returned.

  -- ON ENTRY :
  --   nt        number of tasks for polyhedral homotopy continuation,
  --             if 0, then sequential execution;
  --   silent    if silent, then the computed root counts will not be shown
  --             on screen, otherwise, the root counter remains silent;
  --   p         a polynomial system;
  --   deg       if true, then only degree bounds are used,
  --             otherwise also the mixed volume is computed;
  --   verbose   the verbose level.

  -- ON RETURN :
  --   p         may have been permuted for semi-mixed inputs;
  --   rc        root count, Bezout number or stable mixed volume,
  --             rc = Length_Of(qsols) + Length_Of(qsols0);
  --   q         start system;
  --   qsols     solutions of q, without zero components;
  --   qsols0    solutions of q, with zero components;
  --   rocotime  is elapsed user cpu time for computation of the root counts;
  --   hocotime  is elapsed user cpu time for construction of start system.

  procedure Black_Box_Root_Counting 
               ( nt : in integer32;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 rocos : out Link_to_String;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );
  procedure Black_Box_Root_Counting 
               ( nt : in integer32;
                 p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 rocos : out Link_to_String;
                 q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out DoblDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );
  procedure Black_Box_Root_Counting 
               ( nt : in integer32;
                 p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 rocos : out Link_to_String;
                 q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out QuadDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Calculates four different root counts: total degree, m-homogeneous
  --   Bezout number, generalized Bezout number based on set structure,
  --   and mixed volume.  Heuristics are used for the Bezout numbers.
  --   Returns the start system with lowest root count and least amount
  --   of work, which means that linear-product start systems are preferred,
  --   when Bezout numbers equal the mixed volume.
  --   In case the mixed volume is zero because of a equation with only
  --   one monomial, the smallest Bezout bound is returned.

  -- ON ENTRY :
  --   nt        number of tasks for polyhedral homotopy continuation,
  --             if 0, then sequential execution;
  --   p         a polynomial system;
  --   deg       if true, then only degree bounds are used,
  --             otherwise also the mixed volume is computed;
  --   verbose   the verbose level.

  -- ON RETURN :
  --   p         may have been permuted for semi-mixed inputs;
  --   rc        root count, Bezout number or stable mixed volume,
  --             rc = Length_Of(qsols) + Length_Of(qsols0);
  --   rocos     string with the information about the root counts,
  --             with the same information as the above procedures
  --             when silent equals false;
  --   q         start system;
  --   qsols     solutions of q, without zero components;
  --   qsols0    solutions of q, with zero components;
  --   rocotime  is elapsed user cpu time for computation of the root counts;
  --   hocotime  is elapsed user cpu time for construction of start system.

  procedure Black_Box_Root_Counting 
               ( file : in file_type; nt : in integer32;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );
  procedure Black_Box_Root_Counting 
               ( file : in file_type; nt : in integer32;
                 p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out DoblDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );
  procedure Black_Box_Root_Counting 
               ( file : in file_type; nt : in integer32;
                 p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out QuadDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Calculates four different root counts: total degree, m-homogeneous
  --   Bezout number, generalized Bezout number based on set structure,
  --   and mixed volume.  Heuristics are used for the Bezout numbers.
  --   Returns the start system with lowest root count and least amount
  --   of work, which means that linear-product start systems are preferred,
  --   when Bezout numbers equal the mixed volume.
  --   In case the mixed volume is zero because of a equation with only
  --   one monomial, the smallest Bezout bound is returned.

  -- ON ENTRY :
  --   file      must be opened for output;
  --   nt        number of tasks for polyhedral homotopy continuation,
  --             if 0, then sequential execution;
  --   p         a polynomial system;
  --   deg       if true, then only degree bounds are used,
  --             otherwise also the mixed volume is computed;
  --   verbose   the verbose level.

  -- ON RETURN :
  --   p         may have been permuted for semi-mixed inputs;
  --   rc        root count, Bezout number or stable mixed volume,
  --             rc = Length_Of(qsols) + Length_Of(qsols0);
  --   q         start system;
  --   qsols     solutions of q, without zero components;
  --   qsols0    solutions of q, with zero components;
  --   rocotime  is elapsed user cpu time for computation of the root counts;
  --   hocotime  is elapsed user cpu time for construction of start system.

-- PIPELINED BLACK BOX ROOT COUNTING :

  procedure Pipelined_Stable_Continuation
              ( silent : in boolean; r : in integer32;
                mtype : in Standard_Integer_Vectors.Link_to_Vector;
                stlb : in double_float;
                lifsup: in Link_to_Array_of_Lists;
                mcc : in Mixed_Subdivision; tmv : in natural32;
                lq : in Standard_Complex_Laur_Systems.Laur_Sys;
                mv,smv : out natural32;
                qsols0 : out Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Extracts the stable mixed cells and applies polyhedral
  --   continuation to compute the solutions with zero coordinates,
  --   in standard double precision.

  -- ON ENTRY :
  --   silent   if not silent, then the mixed volume and the stable
  --            mixed volume are written to screen;
  --   r        number of distinct supports;
  --   mtype    type of mixture of the supports;
  --   stlb     bound on the artificial origin;
  --   lifsup   the lifted supports;
  --   mcc      regular mixed cell configuration for the supports;
  --   tmv      the total volume of cells without artificial origin in mcc;
  --   lq       a random coefficient system.

  -- ON RETURN :
  --   mv       the mixed volume of the supports;
  --   smv      the stable mixed volume;
  --   qsols0   solutions of lq with zero coordinates.

  procedure Pipelined_Stable_Continuation
              ( silent : in boolean; r : in integer32;
                mtype : in Standard_Integer_Vectors.Link_to_Vector;
                stlb : in double_float;
                lifsup: in Link_to_Array_of_Lists;
                mcc : in Mixed_Subdivision; tmv : in natural32;
                lq : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mv,smv : out natural32;
                qsols0 : out DoblDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Extracts the stable mixed cells and applies polyhedral
  --   continuation to compute the solutions with zero coordinates,
  --   in double double precision.

  -- ON ENTRY :
  --   silent   if not silent, then the mixed volume and the stable
  --            mixed volume are written to screen;
  --   r        number of distinct supports;
  --   mtype    type of mixture of the supports;
  --   stlb     bound on the artificial origin;
  --   lifsup   the lifted supports;
  --   mcc      regular mixed cell configuration for the supports;
  --   tmv      the total volume of cells without artificial origin in mcc;
  --   lq       a random coefficient system.

  -- ON RETURN :
  --   mv       the mixed volume of the supports;
  --   smv      the stable mixed volume;
  --   qsols0   solutions of lq with zero coordinates.

  procedure Pipelined_Stable_Continuation
              ( silent : in boolean; r : in integer32;
                mtype : in Standard_Integer_Vectors.Link_to_Vector;
                stlb : in double_float;
                lifsup: in Link_to_Array_of_Lists;
                mcc : in Mixed_Subdivision; tmv : in natural32;
                lq : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mv,smv : out natural32;
                qsols0 : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Extracts the stable mixed cells and applies polyhedral
  --   continuation to compute the solutions with zero coordinates,
  --   in double double precision.

  -- ON ENTRY :
  --   silent   if not silent, then the mixed volume and the stable
  --            mixed volume are written to screen;
  --   r        number of distinct supports;
  --   mtype    type of mixture of the supports;
  --   stlb     bound on the artificial origin;
  --   lifsup   the lifted supports;
  --   mcc      regular mixed cell configuration for the supports;
  --   tmv      the total volume of cells without artificial origin in mcc;
  --   lq       a random coefficient system.

  -- ON RETURN :
  --   mv       the mixed volume of the supports;
  --   smv      the stable mixed volume;
  --   qsols0   solutions of lq with zero coordinates.

  procedure Pipelined_Root_Counting 
               ( nt : in integer32; silent : in boolean;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );
  procedure Pipelined_Root_Counting 
               ( nt : in integer32; silent : in boolean;
                 p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out DoblDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );
  procedure Pipelined_Root_Counting 
               ( nt : in integer32; silent : in boolean;
                 p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out QuadDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Calculates four different root counts: total degree, m-homogeneous
  --   Bezout number, generalized Bezout number based on set structure,
  --   and mixed volume.  Heuristics are used for the Bezout numbers.
  --   For nt >= 2, the random coefficient system is solved with pipelined
  --   polyhedral homotopies and then that is the returned start system,
  --   unless deg is true or the mixed volume equals zero.
  --   In case the mixed volume is zero because of a equation with only
  --   one monomial, the smallest Bezout bound is returned.

  -- ON ENTRY :
  --   nt        number of tasks for polyhedral homotopy continuation,
  --             if 0, then sequential execution;
  --   silent    if silent, then the computed root counts will not be shown
  --             on screen, otherwise, the root counter remains silent;
  --   p         a polynomial system;
  --   deg       if true, then only degree bounds are used,
  --             otherwise also the mixed volume is computed;
  --   verbose   the verbose level.

  -- ON RETURN :
  --   p         may have been permuted for semi-mixed inputs;
  --   rc        root count, Bezout number or stable mixed volume,
  --             rc = Length_Of(qsols) + Length_Of(qsols0);
  --   q         start system;
  --   qsols     solutions of q, without zero components;
  --   qsols0    solutions of q, with zero components;
  --   rocotime  is elapsed user cpu time for computation of the root counts;
  --   hocotime  is elapsed user cpu time for construction of start system.

  procedure Pipelined_Root_Counting 
               ( nt : in integer32;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 rocos : out Link_to_String;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );
  procedure Pipelined_Root_Counting 
               ( nt : in integer32;
                 p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 rocos : out Link_to_String;
                 q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out DoblDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );
  procedure Pipelined_Root_Counting 
               ( nt : in integer32;
                 p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 rocos : out Link_to_String;
                 q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out QuadDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Calculates four different root counts: total degree, m-homogeneous
  --   Bezout number, generalized Bezout number based on set structure,
  --   and mixed volume.  Heuristics are used for the Bezout numbers.
  --   For nt >= 2, the random coefficient system is solved with pipelined
  --   polyhedral homotopies and then that is the returned start system,
  --   unless deg is true or the mixed volume equals zero.
  --   In case the mixed volume is zero because of a equation with only
  --   one monomial, the smallest Bezout bound is returned.

  -- ON ENTRY :
  --   nt        number of tasks for polyhedral homotopy continuation,
  --             if 0, then sequential execution;
  --   p         a polynomial system;
  --   deg       if true, then only degree bounds are used,
  --             otherwise also the mixed volume is computed;
  --   verbose   the verbose level.

  -- ON RETURN :
  --   p         may have been permuted for semi-mixed inputs;
  --   rc        root count, Bezout number or stable mixed volume,
  --             rc = Length_Of(qsols) + Length_Of(qsols0);
  --   rocos     string with the information about the root counts,
  --             with the same information as the above procedures
  --             when silent equals false;
  --   q         start system;
  --   qsols     solutions of q, without zero components;
  --   qsols0    solutions of q, with zero components;
  --   rocotime  is elapsed user cpu time for computation of the root counts;
  --   hocotime  is elapsed user cpu time for construction of start system.

  procedure Pipelined_Root_Counting 
               ( file : in file_type; nt : in integer32;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );
  procedure Pipelined_Root_Counting 
               ( file : in file_type; nt : in integer32;
                 p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out DoblDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );
  procedure Pipelined_Root_Counting 
               ( file : in file_type; nt : in integer32;
                 p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 deg : in boolean; rc : out natural32;
                 q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 qsols,qsols0 : out QuadDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Calculates four different root counts: total degree, m-homogeneous
  --   Bezout number, generalized Bezout number based on set structure,
  --   and mixed volume.  Heuristics are used for the Bezout numbers.
  --   For nt >= 2, the random coefficient system is solved with pipelined
  --   polyhedral homotopies and then that is the returned start system,
  --   unless deg is true or the mixed volume equals zero.
  --   In case the mixed volume is zero because of a equation with only
  --   one monomial, the smallest Bezout bound is returned.

  -- ON ENTRY :
  --   file      must be opened for output;
  --   nt        number of tasks for polyhedral homotopy continuation,
  --             if 0, then sequential execution;
  --   p         a polynomial system;
  --   deg       if true, then only degree bounds are used,
  --             otherwise also the mixed volume is computed;
  --   verbose   the verbose level.

  -- ON RETURN :
  --   p         may have been permuted for semi-mixed inputs;
  --   rc        root count, Bezout number or stable mixed volume,
  --             rc = Length_Of(qsols) + Length_Of(qsols0);
  --   q         start system;
  --   qsols     solutions of q, without zero components;
  --   qsols0    solutions of q, with zero components;
  --   rocotime  is elapsed user cpu time for computation of the root counts;
  --   hocotime  is elapsed user cpu time for construction of start system.

-- FOR LAURENT SYSTEMS :

  procedure Black_Box_Root_Counting 
               ( nt : in integer32; silent : in boolean;
                 p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out Standard_Complex_Laur_Systems.Laur_Sys;
                 qsols : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );
  procedure Black_Box_Root_Counting 
               ( nt : in integer32; silent : in boolean;
                 p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out DoblDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );
  procedure Black_Box_Root_Counting 
               ( nt : in integer32; silent : in boolean;
                 p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out QuadDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Wrapper to a mixed-volume computation for a Laurent polynomial system.

  -- ON ENTRY :
  --   nt        number of tasks for polyhedral homotopy continuation,
  --             if 0, then sequential execution;
  --   silent    if silent, then the mixed volume will not be written
  --             to screen, otherwise, if silent, then the user will
  --             be shown the mixed volume on screen;
  --   p         a Laurent polynomial system.

  -- ON RETURN :
  --   p         may have been permuted for semi-mixed inputs;
  --   rc        the root count is the mixed volume, equals Lenght_Of(qsols);
  --   q         start system;
  --   qsols     solutions of q, without zero components;
  --   rocotime  is elapsed user cpu time for computation of the root counts;
  --   hocotime  is elapsed user cpu time for construction of start system.

  procedure Black_Box_Root_Counting 
               ( nt : in integer32;
                 p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32; rocos : out Link_to_String;
                 q : out Standard_Complex_Laur_Systems.Laur_Sys;
                 qsols : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );
  procedure Black_Box_Root_Counting 
               ( nt : in integer32;
                 p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32; rocos : out Link_to_String;
                 q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out DoblDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );
  procedure Black_Box_Root_Counting 
               ( nt : in integer32;
                 p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32; rocos : out Link_to_String;
                 q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out QuadDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Wrapper to a mixed-volume computation for a Laurent polynomial system.

  -- ON ENTRY :
  --   nt        number of tasks for polyhedral homotopy continuation,
  --             if 0, then sequential execution;
  --   p         a Laurent polynomial system.

  -- ON RETURN :
  --   p         may have been permuted for semi-mixed inputs;
  --   rc        the root count is the mixed volume, equals Length_Of(qsols);
  --   rocos     information about the root count, identical what the
  --             above procedures write when silent is false;
  --   q         start system;
  --   qsols     solutions of q, without zero components;
  --   rocotime  is elapsed user cpu time for computation of the root counts;
  --   hocotime  is elapsed user cpu time for construction of start system.

  procedure Black_Box_Root_Counting 
               ( file : in file_type; nt : in integer32;
                 p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out Standard_Complex_Laur_Systems.Laur_Sys;
                 qsols : out Standard_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );
  procedure Black_Box_Root_Counting 
               ( file : in file_type; nt : in integer32;
                 p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out DoblDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );
  procedure Black_Box_Root_Counting 
               ( file : in file_type; nt : in integer32;
                 p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out QuadDobl_Complex_Solutions.Solution_List;
                 rocotime,hocotime : out duration;
                 verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Wrapper to a mixed-volume computation for a Laurent polynomial system,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   file      must be opened for output;
  --   nt        number of tasks for polyhedral homotopy continuation,
  --             if 0, then sequential execution;
  --   p         a Laurent polynomial system.

  -- ON RETURN :
  --   p         may have been permuted for semi-mixed inputs;
  --   rc        the root count is the mixed volume, equals Lenght_Of(qsols);
  --   q         start system;
  --   qsols     solutions of q, without zero components;
  --   rocotime  is elapsed user cpu time for computation of the root counts;
  --   hocotime  is elapsed user cpu time for construction of start system.

-- PIPELINING FOR LAURENT SYSTEMS :

  procedure Pipelined_Root_Counting 
               ( nt : in integer32; silent : in boolean;
                 p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out Standard_Complex_Laur_Systems.Laur_Sys;
                 qsols : out Standard_Complex_Solutions.Solution_List;
                 elaptime : out duration; verbose : in integer32 := 0 );
  procedure Pipelined_Root_Counting 
               ( nt : in integer32; silent : in boolean;
                 p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out DoblDobl_Complex_Solutions.Solution_List;
                 elaptime : out duration; verbose : in integer32 := 0 );
  procedure Pipelined_Root_Counting 
               ( nt : in integer32; silent : in boolean;
                 p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out QuadDobl_Complex_Solutions.Solution_List;
                 elaptime : out duration; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Wrapper to the pipelined polyhedral homotopies to solve
  --   a random coefficient system for a given Laurent polynomial system,
  --   in double, double double, or quad double precision.
  --   No output is written to screen.

  -- ON ENTRY :
  --   nt        number of tasks for polyhedral homotopy continuation,
  --             nt must be at least 2;
  --   silent    if silent, then the mixed volume will not be written
  --             to screen, otherwise, if silent, then the user will
  --             be shown the mixed volume on screen;
  --   p         a Laurent polynomial system;
  --   verbose   the verbose level.

  -- ON RETURN :
  --   p         may have been permuted for semi-mixed inputs;
  --   rc        the root count is the mixed volume, equals Length_Of(qsols);
  --   q         start system;
  --   qsols     solutions of q, without zero components;
  --   elaptime  is elapsed user cpu time for the computation.

  procedure Pipelined_Root_Counting 
               ( nt : in integer32;
                 p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32; rocos : out Link_to_String;
                 q : out Standard_Complex_Laur_Systems.Laur_Sys;
                 qsols : out Standard_Complex_Solutions.Solution_List;
                 elaptime : out duration; verbose : in integer32 := 0 );
  procedure Pipelined_Root_Counting 
               ( nt : in integer32;
                 p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32; rocos : out Link_to_String;
                 q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out DoblDobl_Complex_Solutions.Solution_List;
                 elaptime : out duration; verbose : in integer32 := 0 );
  procedure Pipelined_Root_Counting 
               ( nt : in integer32;
                 p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32; rocos : out Link_to_String;
                 q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out QuadDobl_Complex_Solutions.Solution_List;
                 elaptime : out duration; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Wrapper to the pipelined polyhedral homotopies to solve
  --   a random coefficient system for a given Laurent polynomial system,
  --   in double, double double, or quad double precision.
  --   No output is written to screen.

  -- ON ENTRY :
  --   nt        number of tasks for polyhedral homotopy continuation,
  --             nt must be at least 2;
  --   p         a Laurent polynomial system;
  --   verbose   the verbose level.

  -- ON RETURN :
  --   p         may have been permuted for semi-mixed inputs;
  --   rc        the root count is the mixed volume, equals Length_Of(qsols);
  --   rocos     information about the root count, identical what the
  --             above procedures write when silent is false;
  --   q         start system;
  --   qsols     solutions of q, without zero components;
  --   elaptime  is elapsed user cpu time for the computation.

  procedure Pipelined_Root_Counting 
               ( file : in file_type; nt : in integer32;
                 p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out Standard_Complex_Laur_Systems.Laur_Sys;
                 qsols : out Standard_Complex_Solutions.Solution_List;
                 elaptime : out duration; verbose : in integer32 := 0 );
  procedure Pipelined_Root_Counting 
               ( file : in file_type; nt : in integer32;
                 p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out DoblDobl_Complex_Solutions.Solution_List;
                 elaptime : out duration; verbose : in integer32 := 0 );
  procedure Pipelined_Root_Counting 
               ( file : in file_type; nt : in integer32;
                 p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 rc : out natural32;
                 q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 qsols : out QuadDobl_Complex_Solutions.Solution_List;
                 elaptime : out duration; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Wrapper to the pipelined polyhedral homotopies to solve
  --   a random coefficient system for a given Laurent polynomial system,
  --   in double, double double, or quad double precision.
  --   Output is written to file.

  -- ON ENTRY :
  --   file      must be opened for output;
  --   nt        number of tasks for polyhedral homotopy continuation,
  --             nt must be at least 2;
  --   p         a Laurent polynomial system.

  -- ON RETURN :
  --   p         may have been permuted for semi-mixed inputs;
  --   rc        the root count is the mixed volume, equals Lenght_Of(qsols);
  --   q         start system;
  --   qsols     solutions of q, without zero components;
  --   elaptime  is elapsed user cpu time for the computation.


  procedure Main ( nt : in natural32; infilename,outfilename : in string;
                   verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --    Defines phc -b -r, the main blackbox root counter procedure.

  -- ON ENTRY :
  --   nt             the number of tasks, if 0 then no multitasking,
  --                  otherwise nt tasks will be used to track the paths;
  --   infilename     the name of the input file;
  --   outfilename    the name of the output file;
  --   verbose        the verbose level.

end Black_Box_Root_Counters;
