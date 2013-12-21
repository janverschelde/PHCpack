procedure mainvali ( infilename,outfilename : in string );

-- DESCRIPTION :
--   This is the routine for validating the solutions of a polynomial system,
--   as called by the central dispatcher.
--   The arguments are the names of the input and output file respectively.

-- FUNCTIONALITY :
--   There are three types of functionality offered:
--   1. Weeding out the solution sets:
--     check residuals to separate solutions from path failures;
--     compute local condition numbers, which solution are regular/singular;
--     see whether solutions are clustered, eventually multiple solutions;
--     compute generators of the solution set in case of symmetry.
--   2. Computation of winding numbers:
--     this requires a homotopy with solution paths at t < 1.
--   3. Validation of polyhedral end game:
--     computation of frequency table of the path directions;
--     this requires an output file of poco where the path directions
--     have been computed.
