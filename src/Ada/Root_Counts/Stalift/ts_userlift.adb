with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Mixed_Volume_Computation;
with Main_Lifting_Functions;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;
with Integer_Mixed_Subdivisions_io;
with Drivers_for_Static_Lifting;

procedure ts_userlift is

-- DESCRIPTION :
--   Standalone test program for the mixed volume computation
--   with a user defined integer lifting function.

  procedure Lift_and_Prune
              ( p : in Poly_Sys;
                mix : out Standard_Integer_Vectors.Link_to_Vector;
                mcc : out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Extracts the supports and computes the mixture.

    sup : Array_of_Lists(p'range) := Create(p);
    perms : Standard_Integer_Vectors.Link_to_Vector;
    n : constant integer32 := p'last;
    r : integer32;

    use Drivers_for_Static_Lifting;

  begin
    Mixed_Volume_Computation.Compute_Mixture(sup,mix,perms);
    r := mix'last;
    put("Number of distinct supports : "); put(r,1); new_line;
    declare
      mixsup : constant Array_of_Lists(mix'range)
             := Mixed_Volume_Computation.Typed_Lists(mix.all,sup);
      lifsup : Array_of_Lists(mix'range);
    begin
      for k in mixsup'range loop
        lifsup(k) := Main_Lifting_Functions.Read_Integer_Lifting(mixsup(k));
      end loop;
      put_line("The lifted supports : "); put(lifsup);
      Integer_Create_Mixed_Cells(n,mix.all,lifsup,mcc);
    end;
  end Lift_and_Prune;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system,
  --   extracts the supports, computes the type of mixture,
  --   and asks the user for a lifting for each point.

    lp : Link_to_Poly_Sys;
    mcc : Mixed_Subdivision;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    n,mv : natural32;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    Lift_and_Prune(lp.all,mix,mcc);
    put_line("The mixed cell configuration :");
    n := natural32(lp'last);
    Integer_Mixed_Subdivisions_io.put(n,mix.all,mcc,mv);
    put("The mixed volume : "); put(mv,1); new_line;
  end Main;

begin
  Main;
end ts_userlift;
