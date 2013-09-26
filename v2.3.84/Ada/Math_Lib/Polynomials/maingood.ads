procedure maingood ( infilename,outfilename : in string );

-- DESCRIPTION :
--   This procedure checks whether the polynomial system given on input is
--   in a valid input format, by parsing it as a Laurent polynomial system.
--   This procedure gets called with the option '-g'.
--   Since phc -b adds to the original input file, the phc -g input copy
--   will create a copy of the system in the file input to the file copy,
--   along with some parsing diagnostics and a time stamp.
--   For systems that are parsed well, the output file starts with
--   the string version submitted to phc (not the parsed polynomials).
--   This copying and time stamp of phc -g is a good feature
--   to keep a catalog of all systems submitted to PHCpack.
