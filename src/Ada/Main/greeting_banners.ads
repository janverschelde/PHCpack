package Greeting_Banners is

-- DESCRIPTION :
--   Defines banners to greet the user of phc.
--   Exports the version string.

  welcome : constant string :=
    "Welcome to PHC (Polynomial Homotopy Continuation) v2.4.90 20 Mar 2024";

  compban : constant string :=
    "a numerical irreducible decomposition for solution sets";

  enumban : constant string :=
    "SAGBI/Pieri and Littlewood-Richardson homotopies";

  facban : constant string :=
    "Factor a pure dimensional solution set into irreducibles";

  goodban : constant string :=
    "Checking whether a given input system has the right syntax.";

  symbban : constant string :=
    "Writing the contents of the symbol table after reading an input system.";

  hypban : constant string :=
    "Witness Set for Hypersurface cutting with Random Line.";

  mvcban : constant string :=
    "Mixed Volumes, MixedVol, DEMiCs, polyhedral homotopies";

  pocoban : constant string :=
    "Polynomial Continuation by a homotopy in 1 parameter";

  reduban : constant string :=
    "Linear and nonlinear Reduction w.r.t the total degree of the system.";

  rocoban : constant string :=
    "Root counting and construction of start system";

  samban : constant string :=
    "Sample points of a pure dimensional solution set, given a witness set.";

  scalban : constant string :=
    "Equation/variable Scaling on polynomial system and solution list.";

  slvban : constant string :=
    "Running the equation-by-equation solver, no tasking (yet).";

  trackban : constant string :=
    "Tracking Solution Paths with incremental read/write of solutions.";

  adepban : constant string :=
    "Tracking Solution Paths with Algorithmic Differentiation Methods.";

  seriesban : constant string :=
    "Power Series, Pade Approximants";

  veriban : constant string :=
    "Verification, refinement and purification of computed solution lists.";

  witban : constant string :=
    "Witness Set Intersection using Diagonal Homotopies.";

  function Version return string;

  -- DESCRIPTION :
  --   Returns a string with the version number and release date.

  procedure How_to_Cite;

  -- DESCRIPTION :
  --   Writes how to cite PHCpack in publications.

  procedure show_help;

  -- DESCRIPTION :
  --   Writes essential information about the use of phc to screen.

  procedure help4setseed;

  -- DESCRIPTION :
  --   Writes help on setting the seed in the random number generators.

  procedure help4eqnbyeqn;

  -- DESCRIPTION :
  --   Writes information about the equation-by-equation solver.

  procedure help4blackbox;

  -- DESCRIPTION :
  --   Writes information about the blackbox solver.

  procedure help4compsolve;

  -- DESCRIPTION :
  --   Writes information about computing a numerical irreducible
  --   decomposition in blackbox mode.

  procedure help4components;

  -- DESCRIPTION :
  --   Writes information about the numerical irreducible decomposition.

  procedure help4reduction;

  -- DESCRIPTION :
  --   Writes information about degree reduction.

  procedure help4enumeration;

  -- DESCRIPTION :
  --   Writes information about numerical Schubert calculus.

  procedure help4factor;

  -- DESCRIPTION :
  --   Writes help on the factorization into irreducible components.

  procedure help4feedback;

  -- DESCRIPTION :
  --   Writes help on the computation of feedback laws.

  procedure help4goodformat;

  -- DESCRIPTION :
  --   Writes help on the checking of input formats.

  procedure help4help;

  -- DESCRIPTION :
  --   Writes help on the help system.

  procedure help4hypersurface;

  -- DESCRIPTION :
  --   Writes help for a witness set for a hypersurface.

  procedure help4symbols;

  -- DESCRIPTION :
  --   Writes help on the writing of the symbol table.

  procedure help4continuation;

  -- DESCRIPTION :
  --   Writes help on the polynomial continuation.

  procedure help4jumpstart;

  -- DESCRIPTION :
  --   Writes help on the jumpstarting path tracking.

  procedure help4adepath;

  -- DESCRIPTION :
  --   Writes help for the path trackers with algorithmic differentiation.

  procedure help4mixvol;

  -- DESCRIPTION :
  --   Writes help on the mixed volume computation and polyhedral homotopies.

  procedure help4rootcounts;

  -- DESCRIPTION :
  --   Writes help on the root counting and start system construction.

  procedure help4scaling;

  -- DESCRIPTION :
  --   Writes help on equation and variable scaling.

  procedure help4tasking;

  -- DESCRIPTION :
  --   Writes help on multitasking.

  procedure help4series;

  -- DESCRIPTION :
  --   Writes help on Newton's method for power series solutions.

  procedure help4verification;

  -- DESCRIPTION :
  --   Writes help on the verification of lists of solutions.

  procedure help4verbose;

  -- DESCRIPTION :
  --   Writes help on the verbose option.

  procedure help4witsetinsect;

  -- DESCRIPTION :
  --   Writes help on witness set intersection with diagonal homotopies.

  procedure help4pythondict;

  -- DESCRIPTION :
  --   Writes help on converting solutions to Python dictionaries.

  procedure help4sampler;

  -- DESCRIPTION :
  --   Writes help on sampling points on a positive dimensional solution set.

  procedure help4mapleform;

  -- DESCRIPTION :
  --   Writes help on converting solutions to Maple format.

  procedure help4getstart;

  -- DESCRIPTION :
  --   Writes help on extracting the start system and start solutions
  --   from the output of the blackbox solver.

end Greeting_Banners;
