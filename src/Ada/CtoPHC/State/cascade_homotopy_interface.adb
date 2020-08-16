with text_io;                           use text_io;
with Interfaces.C;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Witness_Sets;
with Standard_PolySys_Container;
with Standard_Solutions_Container;
with PHCpack_Operations;

package body Cascade_Homotopy_Interface is

  function Cascade_Homotopy_Standard_Polynomial
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in cascade_homotopy_interface.");
      put_line("Cascade_Homotopy_Standard_Polynomial ...");
    end if;
    PHCpack_Operations.Standard_Cascade_Homotopy;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cascade_homotopy_interface.");
        put_line("Cascade_Homotopy_Standard_Polynomial.");
      end if;
      return 164;
  end Cascade_Homotopy_Standard_Polynomial;

  function Cascade_Homotopy_DoblDobl_Polynomial
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in cascade_homotopy_interface.");
      put_line("Cascade_Homotopy_DoblDobl_Polynomial ...");
    end if;
    PHCpack_Operations.DoblDobl_Cascade_Homotopy;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cascade_homotopy_interface.");
        put_line("Cascade_Homotopy_DoblDobl_Polynomial.");
      end if;
      return 178;
  end Cascade_Homotopy_DoblDobl_Polynomial;

  function Cascade_Homotopy_QuadDobl_Polynomial
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in cascade_homotopy_interface.");
      put_line("Cascade_Homotopy_QuadDobl_Polynomial ...");
    end if;
    PHCpack_Operations.QuadDobl_Cascade_Homotopy;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cascade_homotopy_interface.");
        put_line("Cascade_Homotopy_QuadDobl_Polynomial.");
      end if;
      return 188;
  end Cascade_Homotopy_QuadDobl_Polynomial;

  function Cascade_Homotopy_Standard_Laurent
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in cascade_homotopy_interface.");
      put_line("Cascade_Homotopy_Standard_Laurent ...");
    end if;
    PHCpack_Operations.Standard_Cascade_Laurent_Homotopy;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cascade_homotopy_interface.");
        put_line("Cascade_Homotopy_Standard_Laurent.");
      end if;
      return 789;
  end Cascade_Homotopy_Standard_Laurent;

  function Cascade_Homotopy_DoblDobl_Laurent
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in cascade_homotopy_interface.");
      put_line("Cascade_Homotopy_DoblDobl_Laurent ...");
    end if;
    PHCpack_Operations.DoblDobl_Cascade_Laurent_Homotopy;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cascade_homotopy_interface.");
        put_line("Cascade_Homotopy_DoblDobl_Laurent.");
      end if;
      return 790;
  end Cascade_Homotopy_DoblDobl_Laurent;

  function Cascade_Homotopy_QuadDobl_Laurent
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in cascade_homotopy_interface.");
      put_line("Cascade_Homotopy_QuadDobl_Laurent ...");
    end if;
    PHCpack_Operations.QuadDobl_Cascade_Laurent_Homotopy;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cascade_homotopy_interface.");
        put_line("Cascade_Homotopy_QuadDobl_Laurent.");
      end if;
      return 791;
  end Cascade_Homotopy_QuadDobl_Laurent;

  function Cascade_Homotopy_Cut_Slack
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant natural := natural(v_a(v_a'first));
    lp : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
    sols : Standard_Complex_Solutions.Solution_List
         := Standard_Solutions_Container.Retrieve;
    np : Standard_Complex_Poly_Systems.Poly_Sys(lp'first..lp'last-1);

  begin
    if vrblvl > 0 then
      put("-> in cascade_homotopy_interface.");
      put_line("Cascade_Homotopy_Cut_Slack ...");
    end if;
   -- put("#equations in systems container : "); put(lp'last,1); new_line;
   -- put("#solutions in solutions container : ");
   -- put(Length_Of(sols),1); new_line;
    if k > 0 then
      Witness_Sets.Remove_Component(sols);
      np := Witness_Sets.Remove_Embedding1(lp.all,1);
      Standard_PolySys_Container.Clear;
      Standard_PolySys_Container.Initialize(np);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in cascade_homotopy_interface.");
        put_line("Cascade_Homotopy_Cut_Slack.");
      end if;
      return 171;
  end Cascade_Homotopy_Cut_Slack;

end Cascade_Homotopy_Interface;
