with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Arrays_of_Floating_Vector_Lists;   use Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Floating_Mixed_Subdivisions;       use Floating_Mixed_Subdivisions;

package Pipelined_Polyhedral_Trackers is

-- DESCRIPTION :
--   This package offers multitasked path trackers to solve a random
--   coefficient system with a polyhedral homotopy.
--   The path tracking is interlaced with the production of the 
--   mixed cells, computed by the MixedVol Algorithm.

-- AFTER PREPROCESSING AND LIFTING, ON LAURENT SYSTEMS, SILENT :

  procedure Silent_Multitasking_Tracker
              ( nt,nbequ,r : in integer32; stlb : in double_float;
                mtype,perm,idx : in Standard_Integer_Vectors.Link_to_Vector;
                vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                lft : in Standard_Floating_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                sols : out Standard_Complex_Solutions.Solution_List );
  procedure Silent_Multitasking_Tracker
              ( nt,nbequ,r : in integer32; stlb : in double_float;
                mtype,perm,idx : in Standard_Integer_Vectors.Link_to_Vector;
                vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                lft : in Standard_Floating_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Silent_Multitasking_Tracker
              ( nt,nbequ,r : in integer32; stlb : in double_float;
                mtype,perm,idx : in Standard_Integer_Vectors.Link_to_Vector;
                vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                lft : in Standard_Floating_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Uses nt tasks to solve a random coefficient system q,
  --   in standard double, double double or quad double precision.
  --   The silent version produce no intermediate output.
  --   This procedure is called after the preprocessing and the lifting
  --   done by mv_upto_pre4mv and mv_lift for the MixedVol Algorithm.

  -- ON ENTRY :
  --   nt       number of tasks;
  --   nbequ    number of equations and variables;
  --   r        number of distinct supports;
  --   stlb     lifting of the artificial origin, 0.0 if no stable mv;
  --   mtype    the type of mixture of the supports;
  --   perm     permutation of the supports;
  --   idx      index to the vertex set;
  --   vtx      vertex points;
  --   lft      lifting values for the vertex points.
 
  -- ON RETURN :
  --   lif      configuration of lifted supports;
  --   mcc      a regular mixed-cell configuration;
  --   mv       the mixed volume of the polytopes spanned by the supports;
  --   q        a random coefficient system;
  --   sols     all solutions of the system q.

-- AFTER PREPROCESSING AND LIFTING, ON LAURENT SYSTEMS, REPORTING :

  procedure Reporting_Multitasking_Tracker
              ( file : in file_type;
                nt,nbequ,r : in integer32; stlb : in double_float;
                mtype,perm,idx : in Standard_Integer_Vectors.Link_to_Vector;
                vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                lft : in Standard_Floating_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                sols : out Standard_Complex_Solutions.Solution_List );
  procedure Reporting_Multitasking_Tracker
              ( file : in file_type;
                nt,nbequ,r : in integer32; stlb : in double_float;
                mtype,perm,idx : in Standard_Integer_Vectors.Link_to_Vector;
                vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                lft : in Standard_Floating_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Reporting_Multitasking_Tracker
              ( file : in file_type;
                nt,nbequ,r : in integer32; stlb : in double_float;
                mtype,perm,idx : in Standard_Integer_Vectors.Link_to_Vector;
                vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                lft : in Standard_Floating_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Uses nt tasks to solve a random coefficient system q,
  --   in standard double, double double or quad double precision.
  --   The reporting version allows to monitor the computations.
  --   This procedure is called after the preprocessing and the lifting
  --   done by mv_upto_pre4mv and mv_lift for the MixedVol Algorithm.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   nt       number of tasks;
  --   nbequ    number of equations and variables;
  --   r        number of distinct supports;
  --   stlb     lifting of the artificial origin, 0.0 if no stable mv;
  --   mtype    the type of mixture of the supports;
  --   perm     permutation of the supports;
  --   idx      index to the vertex set;
  --   vtx      vertex points;
  --   lft      lifting values for the vertex points.
 
  -- ON RETURN :
  --   lif      configuration of lifted supports;
  --   mcc      a regular mixed-cell configuration;
  --   mv       the mixed volume of the polytopes spanned by the supports;
  --   q        a random coefficient system;
  --   sols     all solutions of the system q.

-- DOES PREPROCESSING AND LIFTING, ON LAURENT SYSTEMS, SILENT :

  procedure Silent_Multitasking_Tracker
              ( nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                sols : out Standard_Complex_Solutions.Solution_List );
  procedure Silent_Multitasking_Tracker
              ( nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Silent_Multitasking_Tracker
              ( nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Uses nt tasks to solve a random coefficient system q,
  --   in standard double, double double, or quad double precision.
  --   The silent version produces no intermediate output.
  --   This procedure calls the preprocessing and the lifting done
  --   respectively by the mv_upto_pre4mv and mv_lift procedures in
  --   the MixedVol Algorithm, before calling the other silent tracker.

  -- ON ENTRY :
  --   nt       number of tasks;
  --   nbequ    number of equations and variables;
  --   nbpts    the total number of points in the supports;
  --   ind      ind(k) marks the beginning of the k-th support;
  --   cnt      cnt(k) counts the number of points in the k-th support;
  --   support  vector range 1..nbequ*nbpts with the coordinates of
  --            all points in the supports.
 
  -- ON RETURN :
  --   r        number of distinct supports;
  --   mtype    the type of mixture of the supports;
  --   perm     permutation of the supports;
  --   lif      configuration of lifted supports;
  --   mcc      a regular mixed-cell configuration;
  --   mv       the mixed volume of the polytopes spanned by the supports;
  --   q        a random coefficient system;
  --   sols     all solutions of the system q.

-- DOES PREPROCESSING AND LIFTING, ON LAURENT SYSTEMS, REPORTING :

  procedure Reporting_Multitasking_Tracker
              ( file : in file_type; nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                sols : out Standard_Complex_Solutions.Solution_List );
  procedure Reporting_Multitasking_Tracker
              ( file : in file_type; nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Reporting_Multitasking_Tracker
              ( file : in file_type; nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Uses nt tasks to solve a random coefficient system q,
  --   in standard double, double double, or quad double precision.
  --   The reporting version allows monitoring the computations.
  --   This procedure calls the preprocessing and the lifting done
  --   respectively by the mv_upto_pre4mv and mv_lift procedures in the
  --   MixedVol Algorithm, before calling the other reporting tracker.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   nt       number of tasks;
  --   nbequ    number of equations and variables;
  --   nbpts    the total number of points in the supports;
  --   ind      ind(k) marks the beginning of the k-th support;
  --   cnt      cnt(k) counts the number of points in the k-th support;
  --   support  vector range 1..nbequ*nbpts with the coordinates of
  --            all points in the supports.
 
  -- ON RETURN :
  --   r        number of distinct supports;
  --   mtype    the type of mixture of the supports;
  --   perm     permutation of the supports;
  --   lif      configuration of lifted supports;
  --   mcc      a regular mixed-cell configuration;
  --   mv       the mixed volume of the polytopes spanned by the supports;
  --   q        a random coefficient system;
  --   sols     all solutions of the system q.

-- DOES PREPROCESSING AND LIFTING, STABLE MIXED VOLUMES, SILENT :

  procedure Silent_Multitasking_Tracker
              ( nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                stable : in boolean; stlb : in double_float;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                sols : out Standard_Complex_Solutions.Solution_List );
  procedure Silent_Multitasking_Tracker
              ( nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                stable : in boolean; stlb : in double_float;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Silent_Multitasking_Tracker
              ( nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                stable : in boolean; stlb : in double_float;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Uses nt tasks to solve a random coefficient system q,
  --   in standard double, double double, or quad double precision.
  --   The silent version produces no intermediate output.
  --   This procedure calls the preprocessing and the lifting done
  --   respectively by the mv_upto_pre4mv and mv_lift procedures in
  --   the MixedVol Algorithm, before calling the other silent tracker.

  -- ON ENTRY :
  --   nt       number of tasks;
  --   nbequ    number of equations and variables;
  --   nbpts    the total number of points in the supports;
  --   ind      ind(k) marks the beginning of the k-th support;
  --   cnt      cnt(k) counts the number of points in the k-th support;
  --   support  vector range 1..nbequ*nbpts with the coordinates of
  --            all points in the supports;
  --   stable   indicates whether stable mixed volumes are requested;
  --   stlb     lifting bound to use for stable mixed volumes,
  --            equals 0.0 if no stable mixed volumes are requested.
 
  -- ON RETURN :
  --   r        number of distinct supports;
  --   mtype    the type of mixture of the supports;
  --   perm     permutation of the supports;
  --   lif      configuration of lifted supports;
  --   mcc      a regular mixed-cell configuration;
  --   mv       the mixed volume of the polytopes spanned by the supports;
  --   q        a random coefficient system;
  --   sols     all solutions of the system q.

-- DOES PREPROCESSING AND LIFTING, STABLE MIXED VOLUMES, REPORTING :

  procedure Reporting_Multitasking_Tracker
              ( file : in file_type; nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                stable : in boolean; stlb : in double_float;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                sols : out Standard_Complex_Solutions.Solution_List );
  procedure Reporting_Multitasking_Tracker
              ( file : in file_type; nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                stable : in boolean; stlb : in double_float;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Reporting_Multitasking_Tracker
              ( file : in file_type; nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                stable : in boolean; stlb : in double_float;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Uses nt tasks to solve a random coefficient system q,
  --   in standard double, double double, or quad double precision.
  --   The reporting version allows monitoring the computations.
  --   This procedure calls the preprocessing and the lifting done
  --   respectively by the mv_upto_pre4mv and mv_lift procedures in the
  --   MixedVol Algorithm, before calling the other reporting tracker.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   nt       number of tasks;
  --   nbequ    number of equations and variables;
  --   nbpts    the total number of points in the supports;
  --   ind      ind(k) marks the beginning of the k-th support;
  --   cnt      cnt(k) counts the number of points in the k-th support;
  --   support  vector range 1..nbequ*nbpts with the coordinates of
  --            all points in the supports;
  --   stable   indicates whether stable mixed volumes are requested;
  --   stlb     lifting bound to use for stable mixed volumes,
  --            equals 0.0 if no stable mixed volumes are requested.
 
  -- ON RETURN :
  --   r        number of distinct supports;
  --   mtype    the type of mixture of the supports;
  --   perm     permutation of the supports;
  --   lif      configuration of lifted supports;
  --   mcc      a regular mixed-cell configuration;
  --   mv       the mixed volume of the polytopes spanned by the supports;
  --   q        a random coefficient system;
  --   sols     all solutions of the system q.

end Pipelined_Polyhedral_Trackers;
