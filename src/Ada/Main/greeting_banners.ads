package Greeting_Banners is

-- DESCRIPTION :
--   Defines banners to greet the user of phc.
--   Exports the version string.

  welcome : constant string :=
    "Welcome to PHC (Polynomial Homotopy Continuation) v2.4.19  3 Jun 2016";

  compban : constant string :=
    "a numerical irreducible decomposition for solution sets";

  enumban : constant string :=
    "SAGBI/Pieri homotopies to solve a problem in enumerative geometry.";

  facban : constant string :=
    "Factor a pure dimensional solution set into irreducible components.";

  goodban : constant string :=
    "Checking whether a given input system has the right syntax.";

  symbban : constant string :=
    "Writing the contents of the symbol table after reading an input system.";

  hypban : constant string :=
    "Witness Set for Hypersurface cutting with Random Line.";

  mvcban : constant string :=
    "Mixed-Volume Computation, MixedVol, polyhedral homotopies";

  pocoban : constant string :=
    "Polynomial Continuation by a homotopy in 1 parameter";

  reduban : constant string :=
    "Linear and nonlinear Reduction w.r.t the total degree of the system.";

  rocoban : constant string :=
    "Root counting and Construction of product and polyhedral start systems.";

  samban : constant string :=
    "Sample points of a pure dimensional solution set, given a witness set.";

  scalban : constant string :=
    "Equation/variable Scaling on polynomial system and solution list.";

  slvban : constant string :=
    "Running the equation-by-equation solver, no tasking (yet).";

  trackban : constant string :=
    "Tracking Solution Paths with incremental read/write of solutions.";

  valiban : constant string :=
    "Verification, refinement and purification of computed solution lists.";

  witban : constant string :=
    "Witness Set Intersection using Diagonal Homotopies.";

  function Version return string;

  -- DESCRIPTION :
  --   Returns a string with the version number and release date.

  procedure show_help;

  -- DESCRIPTION :
  --   Writes essential information about the use of phc to screen.

  procedure help4eqnbyeqn;

  -- DESCRIPTION :
  --   Writes information about the equation-by-equation solver.

  procedure help4blackbox;

  -- DESCRIPTION :
  --   Writes information about the blackbox solver.

  procedure help4components;

  -- DESCRIPTION :
  --   Writes information about the numerical irreducible decomposition.

  procedure help4reduction;

  -- DESCRIPTION :
  --   Writes information about degree reduction.

  procedure help4enumeration;

  -- DESCRIPTION :
  --   Writes information about numerical Schubert calculus.

end Greeting_Banners;
