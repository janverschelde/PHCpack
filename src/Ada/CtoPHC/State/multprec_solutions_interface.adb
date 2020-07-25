with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Multprec_Complex_Solutions;
with Multprec_Complex_Solutions_io;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with PHCpack_Operations;
with Multprec_Solutions_Container;

package body Multprec_Solutions_Interface is

  function Multprec_Solutions_Read 
             ( vrblvl : integer32 := 0 ) return integer32 is

    sols : Multprec_Complex_Solutions.Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in Multprec_Solutions_interface.");
      put_line("Multprec_Solutions_Read ...");
    end if;
    Multprec_Complex_Solutions_io.Read(sols);
    Multprec_Solutions_Container.Initialize(sols);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in Multprec_Solutions_interface.");
        put_line("Multprec_Solutions_Read.");
      end if;
      return 450;
  end Multprec_Solutions_Read;

  function Multprec_Solutions_Write
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Multprec_Complex_Solutions;
    use Multprec_Complex_Solutions_io;

    sols : constant Solution_List := Multprec_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in Multprec_Solutions_interface.");
      put_line("Multprec_Solutions_Write ...");
    end if;
    if not Is_Null(sols) then
      if PHCpack_Operations.Is_File_Defined then
        put(PHCpack_Operations.output_file,
            Length_Of(sols),natural32(Head_Of(sols).n),sols);
      else
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in Multprec_Solutions_interface.");
        put_line("Multprec_Solutions_Write.");
      end if;
      return 451;
  end Multprec_Solutions_Write;

  function Multprec_Solutions_Size 
             ( b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in Multprec_Solutions_interface.");
      put_line("Multprec_Solutions_Size ...");
    end if;
    Assign(integer32(Multprec_Solutions_Container.Length),b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in Multprec_Solutions_interface.");
        put_line("Multprec_Solutions_Size.");
      end if;
      return 452;
  end Multprec_Solutions_Size;

  function Multprec_Solutions_Dimension
             ( b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in Multprec_Solutions_interface.");
      put_line("Multprec_Solutions_Dimension ...");
    end if;
    Assign(integer32(Multprec_Solutions_Container.Dimension),b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in Multprec_Solutions_interface.");
        put_line("Multprec_Solutions_Dimension.");
      end if;
      return 453;
  end Multprec_Solutions_Dimension;

  function Multprec_Solutions_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in Multprec_Solutions_interface.");
      put_line("Multprec_Solutions_Clear ...");
    end if;
    Multprec_Solutions_Container.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in Multprec_Solutions_interface.");
        put_line("Multprec_Solutions_Clear.");
      end if;
      return 457;
  end Multprec_Solutions_Clear;

end Multprec_Solutions_Interface;
