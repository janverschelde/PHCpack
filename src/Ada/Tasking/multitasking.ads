with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package Multitasking is

-- DESCRIPTION :
--   Provides an encapsulation to multitasking, along with the definition
--   of a array of booleans to store the status of running tasks,
--   needed for synchronization.

  type boolean_array is array ( integer32 range <> ) of boolean;

  function all_true ( n : integer32; b : boolean_array ) return boolean;
  
  -- DESCRIPTION  :
  --   Returns true if all entries from 1 to n in b are true.

  function all_false ( n : integer32; b : boolean_array ) return boolean;

  -- DESCRIPTION :
  --   Returns true if all entries from 1 to n in b are false;

  function to_string ( n : natural32 ) return string;
  function to_string ( n : integer32 ) return string;

  -- DESCRIPTION :
  --   Returns the string representation of a natural number.

  procedure show ( n : in integer32; b : in boolean_array; s : in string );

  -- DESCRIPTION :
  --   Shows the value for state s for all entries from 1 to n in b.
  --   The string s indicates the meaning of the flags in b.

-- WORKERS executing one job :

  generic

     with procedure Job ( i,n : in integer32 );

     -- DESCRIPTION :
     --   Job to be done by the i-th task, i out of n.

  procedure Silent_Workers ( n : in integer32 );

  -- DESCRIPTION :
  --   Launches as many workers as n.
  --   Every worker calls the procedure Job with its identification number.
  --   This version remains silent.

  generic

     with procedure Job ( i,n : in integer32 );

     -- DESCRIPTION :
     --   Job to be done by the i-th task, i out of n.

  procedure Reporting_Workers ( n : in integer32 );

  -- DESCRIPTION :
  --   Launches as many workers as n.
  --   Every worker calls the procedure Job with its identification number.
  --   This version writes extra information to screen,
  --   allowing to monitor the progress of the jobs.

-- WORKERS executing jobs in loop :

  generic

     with procedure Job ( i,n : in integer32; continue : out boolean );

     -- DESCRIPTION :
     --   Job to be done by the i-th task, i out of n.
     --   Proceeds to the next run of Job if continue is set to true,
     --   otherwise stops if continue is set to false.

  procedure Silent_Looping_Workers ( n : in integer32 );

  -- DESCRIPTION :
  --   Launches as many workers as n.
  --   Every worker calls the procedure Job with its identification number.
  --   This version remains silent.
  --   As long as the continue in the job routine remains true,
  --   the Job is called again after synchronization with all other tasks.

  generic

     with procedure Job ( i,n : in integer32; continue : out boolean );

     -- DESCRIPTION :
     --   Job to be done by the i-th task, i out of n.
     --   Proceeds to the next run of Job is continue is set to true,
     --   otherwise stops if continue is set to false.

  procedure Reporting_Looping_Workers ( n : in integer32 );

  -- DESCRIPTION :
  --   Launches as many workers as n.
  --   Every worker calls the procedure Job with its identification number.
  --   This version writes extra information to screen,
  --   allowing to monitor the progress of the jobs.
  --   As long as the continue in Job remains true, the Job procedure
  --   is called again for each task, after synchronization.

end Multitasking;
