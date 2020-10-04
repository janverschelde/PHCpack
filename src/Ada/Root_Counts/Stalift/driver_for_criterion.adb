with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Communications_with_User;           use Communications_with_User;
with Main_Vertex_Points;                 use Main_Vertex_Points;
with Drivers_for_Mixed_Contributions;    use Drivers_for_Mixed_Contributions;

procedure Driver_for_Criterion
             ( file : in file_type; points : in out Array_of_Lists ) is

  function Menu_for_Criterion return character is

  -- DESCRIPTION :
  --   Shows the menu for computing the set of essential points.

    ans : character;

  begin
    new_line;
    put_line("MENU for removing non-contributing points :");
    put_line("  0. elimination of non-vertex points");
    put_line("  1. no computation of vertex points.");
    put_line("  2. apply simple criterion once");
    put_line("  3. exhaustive sweep through supports");
    put("Make your choice : "); Ask_Alternative(ans,"0123");
    return ans;
  end Menu_for_Criterion;

  procedure Dispatch_Criterion ( choice : character ) is

  -- DESCRIPTION :
  --   Dispatches the selected choice.

    nred : natural32 := 0;

  begin
    if choice /= '1' then
      Vertex_Points(file,points);
      case choice is
        when '2' => Once_Simple_Sweep(file,points,nred);
        when '3' => Full_Simple_Sweep(file,points,nred);
        when others => null;
      end case;
    end if;
  end Dispatch_Criterion;

  procedure Driver is

    choice : constant character := Menu_for_Criterion;

  begin
    Dispatch_Criterion(choice);
  end Driver;

begin
  Driver;
end Driver_for_Criterion;
