procedure mainwit ( witset_one,witset_two,logfile : in string );

-- DESCRIPTION :
--   Main driver to the witness set intersection as called by dispatcher.

-- REQUIRED : the algebraic sets live in the same ambient space,
--   if the file names are empty strings, then the user will be
--   prompted to provide the file names.

-- ON ENTRY :
--   witset_one   witness set for the first algebraic set;
--   witset_two   witness set for the second algebraic set;
--   logfile      file name to write diagnostics on.

-- ON RETURN :
--   logfile_w0, logfile_w1, etc. will contain the super witness sets
--   for components of dimension 0, 1, etc. of the intersection of the
--   two witness sets.
