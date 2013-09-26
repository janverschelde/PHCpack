procedure Dispatch;

-- DESCRIPTION :
--   This procedure scans the arguments of the command line and calls the
--   appropriate drivers.
--   The arguments may be
--     * the names of the input and output file, given in this order;
--     * four different options, preceeded by a hyphen (-).
--   The first type of option is one of the following:
--     s : scal => scaling of the polynomial system
--     d : redu => reduce w.r.t. total degree d
--     p : poco => polynomial continuation
--     r : roco => root counting methods
--     m : mvc  => mixed volume computation
--     v : vali => validation of the solutions
--   The second option is the `-b', to switch to batch processing or
--   black box computation.  This option makes sense in combination with
--   the `-p', `-m' and `-v'.
--   Calling `phc -b' is equivalent to `phc -b -m' followed by `phc -b -p'.
