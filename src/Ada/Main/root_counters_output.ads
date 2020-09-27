with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Natural_Numbers;           use Multprec_Natural_Numbers;
with Partitions_of_Sets_Of_Unknowns;     use Partitions_of_Sets_of_Unknowns;

package Root_Counters_Output is

-- DESCRIPTION :
--   The procedures in the package define the output of the various
--   root counts computed by the black box root counters.

  procedure Write_Total_Degree
               ( file : in file_type;
                 d : in natural64; mp_d : in Natural_Number );
 
  -- DESCRIPTION :
  --   Writes the total degree to file.

  -- ON ENTRY :
  --   file      to write root counts on (could be standard_output);
  --   d         total degree;
  --   mp_d      multiprecision version of total degree (if overflow).

  procedure Write_Mixed_Volumes
               ( file : in file_type; mv,smv : in natural32 );

  -- DESCRIPTION :
  --   Writes the mixed volumes to file.

  -- ON ENTRY :
  --   file      to write root counts to;
  --   mv        mixed volume;
  --   smv       stable mixed volume.

  procedure Write_Root_Counts
               ( file : in file_type; no_mv : in boolean;
                 d : in natural64; mp_d : in Natural_Number;
                 m : in natural32; bz,bs : in natural64;
                 mv,smv : in natural32; z : in Partition );

  -- DESCRIPTION :
  --   Writes root counts and set structures to file.

  -- ON ENTRY :
  --   file      to write root counts to;
  --   no_mv     if no mixed volume was computed;
  --   d         total degree;
  --   mp_d      multiprecision version of total degree (if overflow);
  --   m         number of sets in the m-homogeneous Bezout number;
  --   bz        m-homogeneous Bezout number;
  --   bs        set structure Bezout bound;
  --   mv        mixed volume;
  --   smv       stable mixed volume;
  --   z         partition of variables for m-homogeneous Bezout number.

  function Root_Counts_to_String
               ( no_mv : in boolean;
                 d : in natural64; mp_d : in Natural_Number;
                 m : in natural32; bz,bs : in natural64;
                 mv,smv : in natural32; z : in Partition ) return string;

  -- DESCRIPTION :
  --   Writes root counts and set structures to file to string.

  -- ON ENTRY :
  --   no_mv     if no mixed volume was computed;
  --   d         total degree;
  --   mp_d      multiprecision version of total degree (if overflow);
  --   m         number of sets in the m-homogeneous Bezout number;
  --   bz        m-homogeneous Bezout number;
  --   bs        set structure Bezout bound;
  --   mv        mixed volume;
  --   smv       stable mixed volume;
  --   z         partition of variables for m-homogeneous Bezout number.

  function Mixed_Volumes_to_String
             ( d : in natural64; mv,smv : in natural32 ) return string;

  -- DESCRIPTION :
  --   Writes total degree d, mixed volume mv, and stable mixed volume smv
  --   to a string, which is returned.

  function Mixed_Volumes_to_String
             ( mv,smv : in natural32 ) return string;

  -- DESCRIPTION :
  --   Writes the mixed volume mv and stable mixed volume smv
  --   to a string, which is returned.

end Root_Counters_Output;
