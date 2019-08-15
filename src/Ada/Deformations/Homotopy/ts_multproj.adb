with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Characters_and_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Symbol_Table;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;
with Partitions_of_Sets_of_Unknowns_io;  use Partitions_of_Sets_of_Unknowns_io;
with Multi_Projective_Transformations;   use Multi_Projective_Transformations;

procedure ts_multproj is

-- DESCRIPTION :
--   Development of multiprojective coordinate transformations.

  procedure Standard_Test
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                m : in natural32 ) is

  -- DESCRIPTION :
  --   Runs a test on multi-projective transformation on a system
  --   with coefficients in double precision.

    mpart : constant Standard_Natural_Vectors.Vector := iget(m);
    nbr : constant natural32 := Symbol_Table.Number;
    spart : constant Partition := Make_Partition(nbr,m,mpart);
    deg : Standard_Integer_Vectors.Vector(1..integer32(m));
    mhp : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    lhp : constant Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(m))
        := Standard_Random_Linear_Polynomials(nbr,m,spart);
    shp : constant Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(m))
        := Standard_Start_Linear_Polynomials(nbr,m);

  begin
    put("The partition :"); put(mpart); new_line;
    put("Its symbolic form : "); put(spart); new_line;
    for i in p'range loop
      deg := Multiset_Degrees(p(i),m,spart);
      put("degrees of polynomial "); put(i,1);
      put(" :"); put(deg); new_line;
    end loop;
    Symbol_Table.Enlarge(m);
    for i in 1..m loop
      declare
        sv : constant string := "Z" & Characters_and_Numbers.nConvert(i);
      begin
        Symbol_Table.Add_String(sv);
      end;
    end loop;
    mhp := Make_Homogeneous(p,m,spart);
    put("The "); put(m,1); put_line("-homogeneous form :"); put(mhp);
    put_line("The linear equations :"); put_line(lhp);
    put_line("The start equations :"); put(shp);
  end Standard_Test;

  procedure DoblDobl_Test
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                m : in natural32 ) is

  -- DESCRIPTION :
  --   Runs a test on multi-projective transformation on a system
  --   with coefficients in double double precision.

    mpart : constant Standard_Natural_Vectors.Vector := iget(m);
    nbr : constant natural32 := Symbol_Table.Number;
    spart : constant Partition := Make_Partition(nbr,m,mpart);
    deg : Standard_Integer_Vectors.Vector(1..integer32(m));
    mhp : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    lhp : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m))
        := DoblDobl_Random_Linear_Polynomials(nbr,m,spart);
    shp : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m))
        := DoblDobl_Start_Linear_Polynomials(nbr,m);

  begin
    put("The partition :"); put(mpart); new_line;
    put("Its symbolic form : "); put(spart); new_line;
    for i in p'range loop
      deg := Multiset_Degrees(p(i),m,spart);
      put("degrees of polynomial "); put(i,1);
      put(" :"); put(deg); new_line;
    end loop;
    Symbol_Table.Enlarge(m);
    for i in 1..m loop
      declare
        sv : constant string := "Z" & Characters_and_Numbers.nConvert(i);
      begin
        Symbol_Table.Add_String(sv);
      end;
    end loop;
    mhp := Make_Homogeneous(p,m,spart);
    put("The "); put(m,1); put_line("-homogeneous form :"); put(mhp);
    put_line("The linear equations :"); put_line(lhp);
    put_line("The start equations :"); put(shp);
  end DoblDobl_Test;

  procedure QuadDobl_Test
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                m : in natural32 ) is

  -- DESCRIPTION :
  --   Runs a test on multi-projective transformation on a system
  --   with coefficients in quad double precision.

    mpart : constant Standard_Natural_Vectors.Vector := iget(m);
    nbr : constant natural32 := Symbol_Table.Number;
    spart : constant Partition := Make_Partition(nbr,m,mpart);
    deg : Standard_Integer_Vectors.Vector(1..integer32(m));
    mhp : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    lhp : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m))
        := QuadDobl_Random_Linear_Polynomials(nbr,m,spart);
    shp : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m))
        := QuadDobl_Start_Linear_Polynomials(nbr,m);

  begin
    put("The partition :"); put(mpart); new_line;
    put("Its symbolic form : "); put(spart); new_line;
    for i in p'range loop
      deg := Multiset_Degrees(p(i),m,spart);
      put("degrees of polynomial "); put(i,1);
      put(" :"); put(deg); new_line;
    end loop;
    Symbol_Table.Enlarge(m);
    for i in 1..m loop
      declare
        sv : constant string := "Z" & Characters_and_Numbers.nConvert(i);
      begin
        Symbol_Table.Add_String(sv);
      end;
    end loop;
    mhp := Make_Homogeneous(p,m,spart);
    put("The "); put(m,1); put_line("-homogeneous form :"); put(mhp);
    put_line("The linear equations :"); put_line(lhp);
    put_line("The start equations :"); put(shp);
  end QuadDobl_Test;

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system,
  --   interactively defines a partition for the set of unknowns,
  --   and then runs a test in standard double precision.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    m : natural32 := 0;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    new_line;
    put("-> read "); put(lp'last,1); put_line(" polynomials.");
    new_line;
    put("Give the number of sets in the partition : "); get(m);
    Standard_Test(lp.all,m);
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system,
  --   interactively defines a partition for the set of unknowns,
  --   and then runs a test in double double precision.

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    m : natural32 := 0;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    new_line;
    put("-> read "); put(lp'last,1); put_line(" polynomials.");
    new_line;
    put("Give the number of sets in the partition : "); get(m);
    DoblDobl_Test(lp.all,m);
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system,
  --   interactively defines a partition for the set of unknowns,
  --   and then runs a test in quad double precision.

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    m : natural32 := 0;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    new_line;
    put("-> read "); put(lp'last,1); put_line(" polynomials.");
    new_line;
    put("Give the number of sets in the partition : "); get(m);
    QuadDobl_Test(lp.all,m);
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the level of precision.

    ans : character;
  
  begin
    new_line;
    put_line("MENU for the working precision to test multihomogenizations :");
    put_line("  0. standard double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Main;
      when '1' => DoblDobl_Main;
      when '2' => QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_multproj;
