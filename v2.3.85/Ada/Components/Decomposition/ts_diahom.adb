with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Drivers_to_Intrinsic_Solvers;      use Drivers_to_Intrinsic_Solvers;
with Extrinsic_Diagonal_Solvers;        use Extrinsic_Diagonal_Solvers;

procedure ts_diahom is

-- DESCRIPTION :
--   Test drivers for extrinsic and intrinsic versions of the diagonal
--   homotopies to intersect pure dimensional solution sets.

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Diagonal homotopies to compute witness points on components");
    new_line;
    put_line("MENU to test the diagonal homotopies :");
    put_line("  Bring polynomial into complete intersection format :");
    put_line("    1. Randomize to have as many equations as co-dimension.");
    put_line("  Manipulate diagonal homotopies extrinsically : ");
    put_line("    2. Build a diagonal cascade from two embedded systems.");
    put_line("    3. Collapse a diagonal system eliminating the diagonal.");
    put_line("  Intrinsic version provides incremental solver :");
    put_line("    4. Intersect two hypersurfaces, given by two polynomials.");
    put_line("    5. Compute witness points for given polynomial system.");
    put("Type 1, 2, 3, 4, or 5 to make your choice : ");
    Ask_Alternative(ans,"12345");
    case ans is
      when '1' => Randomize_System;
      when '2' => Build_Diagonal_Cascade;
      when '3' => Collapse_Diagonal_System;
      when '4' => Intersection_of_Two_Hypersurfaces;
      when '5' => Driver_to_Hypersurface_Solvers;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_diahom;
