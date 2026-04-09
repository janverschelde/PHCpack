with Interfaces.C;
with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_VecVecs;
with Standard_Complex_VecVecs;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Standard_LaurSys_Container;
with Double_VecVecs_Container;
with DCMPLX_VecVecs_Container;
with Real_Powered_Homotopy_IO;

package body Double_Puiseux_Interface is

  procedure Show_Data ( vrblvl : in integer32 := 0) is

  -- DESCRIPTION :
  --   Retrieves the Laurent polynomial system and the corresponding
  --   coefficients of the series coefficients and the arrays of vectors.
  --   Writes the series system to screen as a test.

    p : constant Standard_Complex_Laur_Systems.Link_to_Laur_Sys
      := Standard_LaurSys_Container.Retrieve;
    pwr : constant Standard_Floating_VecVecs.Link_to_Array_of_VecVecs
        := Double_VecVecs_Container.Get(vrblvl-1);
    cff : constant Standard_Complex_VecVecs.Link_to_Array_of_VecVecs
        := DCMPLX_VecVecs_Container.Get(vrblvl-1);

    use Standard_Complex_Laur_Systems;

  begin
    if vrblvl > 0
     then put_line("-> in double_puiseux_interface.Show_Data ...");
    end if;
    if p = null
     then put_line("The Laurent system is not defined!");
     else put_line("The Laurent system :"); put(p.all);
    end if;
    Real_Powered_Homotopy_IO.put_line
      (p'last,p'last,pwr(pwr'first)'last,p.all,cff.all,pwr.all,
       vrblvl=>vrblvl-1);
  end Show_Data;

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
      if vrblvl > 0 then
        put("  nbr : "); put(nbr,1); put_line(" ...");
        Show_Data(vrblvl);
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
