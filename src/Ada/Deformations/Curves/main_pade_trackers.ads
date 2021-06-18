with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package Main_Pade_Trackers is

-- DESCRIPTION :
--   The main procedures for phc -u are defined by this package.

  procedure Run_Regular_Newton_Puiseux ( valprc : in character );

  -- DESCRIPTION :
  --   Runs the Newton-Puiseux algorithm for the generic case,
  --   in double, double double, or quad double precision,
  --   when valprc is respectively '0', '1', or '4'.

  procedure Run_Power_Series_Newton
              ( infilename,outfilename : in string;
                valprc : in character; vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes power series expansions with Newton's method,
  --   in double, double double, or quad double precision,
  --   when valprc is respectively '0', '1', or '4'.
  --   The value of the verbose level is vrb.

  procedure Run_Path_Trackers
              ( valprc : in character; vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Tracks paths in double, double double, or quad double precision,
  --   when valprc is respectively '0', '1', or '4'.
  --   The value of the verbose level is vrb.

  procedure Run_Path_Convolution_Trackers
              ( nbtasks : in natural32; valprc : in character;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Tracks paths with convolution circuits 
  --   and multitasking if nbtasks > 0,
  --   If nbtasks equals zero, then the number of tasks is prompted for.
  --   in double, double double, or quad double precision,
  --   when valprc is respectively '0', '1', or '4'.
  --   The value of the verbose level is vrb.

  function Prompt_for_Precision_Level return character;

  -- DESCRIPTION :
  --   Prompts for the precision and returns '1', '2', '3', '4', '5', '6',
  --   or '7', respectively for double, double double, triple double,
  --   quad double, penta double, octo double, or deca double precision.

  procedure Run_Regular_Newton_Puiseux;

  -- DESCRIPTION :
  --   Prompts for the precision and then runs the Newton-Puiseux
  --   algorithm for the (sufficiently) generic case.
  --   The verbose level is given in vrb.

  procedure Run_Power_Series_Newton
              ( infilename,outfilename : in string; vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts for the precision and then runs Newton's method
  --   to compute power series expansions.
  --   The verbose level is given in vrb.

  procedure Run_Path_Trackers ( vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts for the precision and then tracks paths.
  --   The verbose level is given in vrb.

  procedure Run_Path_Convolution_Trackers
              ( nbt : in natural32; vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts for the precision and then tracks paths.
  --   The number of tasks is given in nbt.
  --   The verbose level is given in vrb.

  procedure Run_Newton_Fabry
              ( nbtasks : in natural32; prc : in character;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   If prc = '0', then prompts for the precision to compute
  --   the convergence radius estimate of a series solution.
  --   The number of tasks is in nbtasks and the verbose level in vrblvl.

  function Prompt_for_Method return character;

  -- DESCRIPTION :
  --   Displays the menu with methods and returns the character
  --   that corresponds with the selected method.

  procedure Nonzero_Precision_Main
              ( infilename,outfilename : in string;
                nbtasks : in natural32; valprc : in character;
                vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   This main procedure is called when the precision is determined,
  --   either interactively or via the value of precision in valprc.
  --   The value for the precision in valprc is not '0', it is
  --   '1' for double, '2' for double double, '3' for quad double,
  --   '4' for quad double, '5' for penta double, '6' for octo double,
  --   or '7' for deca double,
  --   The procedure can be called immediately if the precision is
  --   set at the command line.

  procedure Zero_Precision_Main
              ( infilename,outfilename : in string;
                nbtasks : in natural32; vrb : in integer32 := 0 );

  -- DESCRIPTION :
  --   This main procedure is called when the precision is undetermined,
  --   The value for the precision will be set after the prompt
  --   for the kind of power series method.

  procedure Main ( infilename,outfilename : in string;
                   nbtasks : in natural32; precision : in character;
                   verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Defines what phc -u does.

  -- ON ENTRY :
  --   nbtasks      number of tasks, 0 for no multitasking;
  --   precision    indicates the precision, is one of the following:
  --                '0' : the user will be prompted for the precision,
  --                '1' : standard double precision,
  --                '2' : double double precision,
  --                '3' : triple double precision,
  --                '4' : quad double precision;
  --                '5' : penta double precision;
  --                '6' : octo double precision;
  --                '7' : deca double precision;
  --   infilename   name of the input file, if "", then the user will
  --                be prompted to provide the name of the input file;
  --   outfilename  name of the output file, if "", then the user will
  --                be prompted to provide the name of the output file;
  --   verbose      the verbose level.

end Main_Pade_Trackers;
