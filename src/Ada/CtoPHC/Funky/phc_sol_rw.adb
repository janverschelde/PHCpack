with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions; 
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io; 
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Integer_Arrays_io;               use C_Integer_Arrays_io;
with Coefficient_Solution_Vectors;      use Coefficient_Solution_Vectors;

function phc_sol_rw ( rw,size_p : integer32; p : C_DblArrs.Pointer )
                    return C_DblArrs.Pointer is

  res : C_DblArrs.Pointer;

  procedure Coefficient_Vector_Representation ( s : in Solution_List ) is

    n : constant integer32 := Head_Of(s).n;
    m : constant C_Integer_Array := Multiplicities(s);
    c : constant C_Double_Array := Coefficients(s);
    cct : C_Double_Array(0..Interfaces.C.size_T(size_p-1))
        := C_DblArrs.Value(p,Interfaces.C.ptrdiff_T(size_p));

  begin
    cct := Concat(natural32(n),m,c);
    res := cct(0)'unchecked_access;
  end Coefficient_Vector_Representation;

  procedure Write_Solution_List is

    cct : constant C_Double_Array
        := C_DblArrs.Value(p,Interfaces.C.ptrdiff_T(size_p));
    size : constant integer32 := integer32(cct(cct'first));
    x : C_Double_Array(0..Interfaces.C.size_t(size))
      := cct(0..Interfaces.C.size_t(size));

  begin
    if size > size_p then
      put_line("WARNING: MORE ELEMENTS THAN SIZE OF BUFFER !!!");
      put(size,1); put(" > "); put(size_p,1); new_line;
    end if;
    declare
      n : constant natural32 := Dimension(x);
      m : constant C_Integer_Array := Multiplicities(cct);
      c : constant C_Double_Array := Coefficients(x);
      s : Solution_List := Create(n,m,c);
    begin
      if not Is_Null(s)
       then put(standard_output,Length_Of(s),natural32(Head_Of(s).n),s);
      end if;
    end;
  end Write_Solution_List;

begin
  if rw = 0 then
    new_line;
    declare
      sols : Solution_List;
    begin
      Read(sols);
      Coefficient_Vector_Representation(sols);
    end;
  else
    Write_Solution_List;
  end if;
  return res;
end phc_sol_rw;
