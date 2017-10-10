with Quad_Double_Numbers;                 use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;            use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Solutions;          use QuadDobl_Complex_Solutions;

package body QuadDobl_Solution_Manipulators is

  procedure Remove_Imaginary_Part ( t : in out Complex_Number ) is

    zero : constant quad_double := create(0.0);

  begin
    t := Create(REAL_PART(t),zero);
  end Remove_Imaginary_Part;

  procedure Remove_Imaginary_Target ( s : in out Link_to_Solution ) is
  begin
    Remove_Imaginary_Part(s.t);
  end Remove_Imaginary_Target;

  procedure Remove_Imaginary_Target ( s : in out Solution_List ) is

    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Remove_Imaginary_Target(ls);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Remove_Imaginary_Target;

end QuadDobl_Solution_Manipulators;
