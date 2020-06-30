with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

procedure mainseries ( nt : in natural32; precision : in character;
                       infilename,outfilename : in string;
                       verbose : in integer32 := 0 );

-- DESCRIPTION :
--   This procedure is called by dispatch.adb to define the phc -u,
--   to compute power series solutions of polynomial systems.
--
-- ON ENTRY :
--   nt           number of tasks, 0 for no multitasking;
--   precision    indicates the precision, is one of the following:
--                '0' : the user will be prompted for the precision,
--                '1' : standard double precision,
--                '2' : double double precision,
--                '4' : quad double precision;
--   infilename   name of the input file, if "", then the user will
--                be prompted to provide the name of the input file;
--   outfilename  name of the output file, if "", then the user will
--                be prompted to provide the name of the output file;
--   verbose      the verbose level.
