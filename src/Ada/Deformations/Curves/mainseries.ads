procedure mainseries ( precision : in character;
                       infilename,outfilename : in string );

-- DESCRIPTION :
--   This procedure is called by dispatch.adb to define the phc -u,
--   to compute power series solutions of polynomial systems.
--
-- ON ENTRY :
--   precision    indicates the precision, is one of the following:
--                '0' : the user will be prompted for the precision,
--                '1' : standard double precision,
--                '2' : double double precision,
--                '4' : quad double precision;
--   infilename   name of the input file, if "", then the user will
--                be prompted to provide the name of the input file;
--   outfilename  name of the output file, if "", then the user will
--                be prompted to provide the name of the output file.
