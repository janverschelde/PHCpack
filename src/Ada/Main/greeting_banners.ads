package Greeting_Banners is

-- DESCRIPTION :
--   Defines banners to greet the user of phc.
--   Exports the version string.

  welcome : constant string :=
    "Welcome to PHC (Polynomial Homotopy Continuation) v2.4.12 19 Feb 2016";

  compban : constant string :=
    "A numerical irreducible decomposition for the solution components.";

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

end Greeting_Banners;
