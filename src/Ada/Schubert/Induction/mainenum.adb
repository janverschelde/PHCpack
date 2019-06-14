with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Drivers_for_SAGBI_Homotopies;       use Drivers_for_SAGBI_Homotopies;
with Drivers_for_Pieri_Homotopies;       use Drivers_for_Pieri_Homotopies;
with Drivers_for_Quantum_Pieri;          use Drivers_for_Quantum_Pieri;
with Drivers_for_Schubert_Induction;     use Drivers_for_Schubert_Induction;

procedure mainenum ( verbose : in integer32 := 0 ) is

  m,p,q,n : natural32 := 0;
  ans : character;

begin
  if verbose > 0 then
    put("At verbose level "); put(verbose,1);
    put_line(", in mainenum ...");
  end if;
  new_line;
  put_line("MENU of Homotopies to solve Enumerative Geometry Problems");
  put_line("  1. SAGBI for intersection hypersurface conditions;");
  put_line("  2. Pieri for hypersurface and general co-dimensions;");
  put_line("  3. Pieri to compute maps of degree q that produce p-planes;");
  put_line("  4. Count solutions to a general Schubert intersection problem;");
  put_line("  5. Compute solutions to Schubert intersection conditions."); 
  put("Type 1, 2, 3, 4, or 5 to select : "); Ask_Alternative(ans,"12345");
  new_line;
  if ans = '4' or ans = '5' then
    put("Give n (dimension of ambient space) : "); get(n);
  else
    put("Give p, dimension of the solution planes : "); get(p);
    put("Give m, the co-dimension so that n = m+p : "); get(m);
    if ans = '3' then
      put("Give q, the degree of the maps : "); get(q);
    end if;
  end if;
  skip_line;
  new_line;
  case ans is
    when '1' => Driver_for_SAGBI_Homotopies(m+p,p);
    when '2' => Driver_for_Pieri_Homotopies(m+p,p);
    when '3' => Driver_for_Quantum_Pieri(m+p,p,q);
    when '4' => Resolve_Intersection_Condition(n);
    when '5' => Solve_Schubert_Problems(integer32(n));
    when others => put_line("Option not recognized.  Please try again...");
  end case;
end mainenum;
