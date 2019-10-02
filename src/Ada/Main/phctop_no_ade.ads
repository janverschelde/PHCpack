procedure phctop_no_ade;

-- DESCRIPTION :
--   This procedure processes the command line arguments
--   and calls the proper procedures corresponding to the options.
--   This is the top level procedure for the executable phc.
--   The "no_ade" means that the QD library will not be needed during
--   the linking, as the Algorithmic Differentiation for the Path library
--   will not be linked in.
