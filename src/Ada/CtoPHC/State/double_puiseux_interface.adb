with Interfaces.C;
with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;

package body Double_Puiseux_Interface is

  function Linear_Solver
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in double_puiseux_interface.Linear_Solver ...");
    end if;
    declare
      v_a : constant C_Integer_Array
          := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
      nbr : constant integer32 := integer32(v_a(v_a'first));
    begin
      if vrblvl > 0
       then put("  nbr : "); put(nbr,1); put_line(" ...");
      end if;
    end;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in double_puiseux_interface.");
        put_line("Linear_Solver.");
      end if;
      return -1;
  end Linear_Solver;

  function Newton_Steps
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in double_puiseux_interface.Newton_Steps ...");
    end if;
    declare
      v_a : constant C_Integer_Array
          := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
      nbr : constant integer32 := integer32(v_a(v_a'first));
    begin
      if vrblvl > 0
       then put("  nbr : "); put(nbr,1); put_line(" ...");
      end if;
    end;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in double_puiseux_interface.");
        put_line("Newton_Steps.");
      end if;
      return -1;
  end Newton_Steps;

end Double_Puiseux_Interface;
