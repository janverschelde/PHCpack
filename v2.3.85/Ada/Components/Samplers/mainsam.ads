procedure mainsam ( witset,logfile : in string );

-- DESCRIPTION :
--   Main driver to the computation of a new witness set,
--   or sampling new points, as called by dispatcher.

-- ON ENTRY :
--   witset       witness set for an algebraic set;
--   logfile      file name to write diagnostics on.

-- ON RETURN :
--   The user is prompted for a an output fill which will contains
--   new samples of an algebraic set, written as witness set.
