with text_io;                           use text_io;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Multprec_Complex_Solutions;
with Standard_PolySys_Container;
with DoblDobl_PolySys_Container;
with QuadDobl_PolySys_Container;
with Multprec_PolySys_Container;
with Standard_Solutions_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_Solutions_Container;
with Multprec_Solutions_Container;
with PHCpack_Operations;

package body Job_Containers is

  function Standard_Copy_Target_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Standard_Copy_Target_System_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Target_System(lp);
    if lp = null then
      return 1;
    else
      Standard_PolySys_Container.Initialize(lp.all);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.Standard_Copy_Target_System_to_Container.");
      end if;
      return 1;
  end Standard_Copy_Target_System_to_Container;

  function DoblDobl_Copy_Target_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("DoblDobl_Copy_Target_System_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Target_System(lp);
    if lp = null then
      return 251;
    else
      DoblDobl_PolySys_Container.Initialize(lp.all);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.DoblDobl_Copy_Target_System_to_Container.");
      end if;
      return 251;
  end DoblDobl_Copy_Target_System_to_Container;

  function QuadDobl_Copy_Target_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("QuadDobl_Copy_Target_System_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Target_System(lp);
    if lp = null then
      return 261;
    else
      QuadDobl_PolySys_Container.Initialize(lp.all);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.DoblDobl_Copy_Target_System_to_Container.");
      end if;
      return 261;
  end QuadDobl_Copy_Target_System_to_Container;

  function Multprec_Copy_Target_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Multprec_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Multprec_Copy_Target_System_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Target_System(lp);
    if lp = null then
      return 281;
    else
      Multprec_PolySys_Container.Initialize(lp.all);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.Multprec_Copy_Target_System_to_Container.");
      end if;
      return 281;
  end Multprec_Copy_Target_System_to_Container;

  function Standard_Copy_Container_System_to_Target
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Standard_Copy_Container_System_to_Target.");
    end if;
    if lp = null then
      return 2;
    else
      PHCpack_Operations.Store_Target_System(lp.all);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.Standard_Copy_Container_System_to_Target.");
      end if;
      return 2;
  end Standard_Copy_Container_System_to_Target;

  function DoblDobl_Copy_Container_System_to_Target
             ( vrblvl : integer32 := 0 ) return integer32 is 

    use DoblDobl_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("DoblDobl_Copy_Container_System_to_Target.");
    end if;
    if lp = null then
      return 252;
    else
      PHCpack_Operations.Store_Target_System(lp.all);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.DoblDobl_Copy_Container_System_to_Target.");
      end if;
      return 252;
  end DoblDobl_Copy_Container_System_to_Target;

  function QuadDobl_Copy_Container_System_to_Target
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("QuadDobl_Copy_Container_System_to_Target.");
    end if;
    if lp = null then
      return 262;
    else
      PHCpack_Operations.Store_Target_System(lp.all);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.QuadDobl_Copy_Container_System_to_Target.");
      end if;
      return 262;
  end QuadDobl_Copy_Container_System_to_Target;

  function Multprec_Copy_Container_System_to_Target
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Multprec_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := Multprec_PolySys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Multprec_Copy_Container_System_to_Target.");
    end if;
    if lp = null then
      return 282;
    else
      PHCpack_Operations.Store_Target_System(lp.all);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.Multprec_Copy_Container_System_to_Target.");
      end if;
      return 282;
  end Multprec_Copy_Container_System_to_Target;

  function Standard_Copy_Start_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Standard_Copy_Start_System_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Start_System(lp);
    if lp = null then
      return 3;
    else
      Standard_PolySys_Container.Initialize(lp.all);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.Standard_Copy_Start_System_to_Container.");
      end if;
      return 3;
  end Standard_Copy_Start_System_to_Container;

  function DoblDobl_Copy_Start_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is


    use DoblDobl_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("DoblDobl_Copy_Start_System_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Start_System(lp);
    if lp = null then
      return 253;
    else
      DoblDobl_PolySys_Container.Initialize(lp.all);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.DoblDobl_Copy_Start_System_to_Container.");
      end if;
      return 253;
  end DoblDobl_Copy_Start_System_to_Container;

  function QuadDobl_Copy_Start_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("QuadDobl_Copy_Start_System_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Start_System(lp);
    if lp = null then
      return 263;
    else
      QuadDobl_PolySys_Container.Initialize(lp.all);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.QuadDobl_Copy_Start_System_to_Container.");
      end if;
      return 263;
  end QuadDobl_Copy_Start_System_to_Container;

  function Multprec_Copy_Start_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Multprec_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Multprec_Copy_Start_System_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Start_System(lp);
    if lp = null then
      return 283;
    else
      Multprec_PolySys_Container.Initialize(lp.all);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.Multprec_Copy_Start_System_to_Container.");
      end if;
      return 283;
  end Multprec_Copy_Start_System_to_Container;

  function Standard_Copy_Container_System_to_Start
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Standard_Copy_Container_System_to_Start.");
    end if;
    if lp = null then
      return 4;
    else
      PHCpack_Operations.Store_Start_System(lp.all);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.Standard_Copy_Container_System_to_Start.");
      end if;
      return 4;
  end Standard_Copy_Container_System_to_Start;

  function DoblDobl_Copy_Container_System_to_Start
             ( vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("DoblDobl_Copy_Container_System_to_Start.");
    end if;
    if lp = null then
      return 254;
    else
      PHCpack_Operations.Store_Start_System(lp.all);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.DoblDobl_Copy_Container_System_to_Start.");
      end if;
      return 254;
  end DoblDobl_Copy_Container_System_to_Start;

  function QuadDobl_Copy_Container_System_to_Start
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("QuadDobl_Copy_Container_System_to_Start.");
    end if;
    if lp = null then
      return 264;
    else
      PHCpack_Operations.Store_Start_System(lp.all);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.QuadDobl_Copy_Container_System_to_Start.");
      end if;
      return 264;
  end QuadDobl_Copy_Container_System_to_Start;

  function Multprec_Copy_Container_System_to_Start
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Multprec_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := Multprec_PolySys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Multprec_Copy_Container_System_to_Start.");
    end if;
    if lp = null then
      return 284;
    else
      PHCpack_Operations.Store_Start_System(lp.all);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.Multprec_Copy_Container_System_to_Start.");
      end if;
      return 284;
  end Multprec_Copy_Container_System_to_Start;

  function Standard_Copy_Target_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions;
    sols : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Standard_Copy_Target_Solutions_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Target_Solutions(sols);
    if Is_Null(sols) then 
      return 5;
    else
      Standard_Solutions_Container.Initialize(sols);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("Standard_Copy_Target_Solutions_to_Container.");
      end if;
      put_line("Exception 5 at copy of target solutions to container.");
      return 5;
  end Standard_Copy_Target_Solutions_to_Container;

  function DoblDobl_Copy_Target_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Solutions;
    sols : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("DoblDobl_Copy_Target_Solutions_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Target_Solutions(sols);
    if Is_Null(sols) then 
      return 255;
    else
      DoblDobl_Solutions_Container.Initialize(sols);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("DoblDobl_Copy_Target_Solutions_to_Container.");
      end if;
      return 255;
  end DoblDobl_Copy_Target_Solutions_to_Container;

  function QuadDobl_Copy_Target_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Solutions;
    sols : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("QuadDobl_Copy_Target_Solutions_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Target_Solutions(sols);
    if Is_Null(sols) then 
      return 265;
    else
      QuadDobl_Solutions_Container.Initialize(sols);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("QuadDobl_Copy_Target_Solutions_to_Container.");
      end if;
      return 265;
  end QuadDobl_Copy_Target_Solutions_to_Container;

  function Multprec_Copy_Target_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Multprec_Complex_Solutions;
    sols : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Multprec_Copy_Target_Solutions_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Target_Solutions(sols);
    if Is_Null(sols) then 
      return 285;
    else
      Multprec_Solutions_Container.Initialize(sols);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("Multprec_Copy_Target_Solutions_to_Container.");
      end if;
      return 285;
  end Multprec_Copy_Target_Solutions_to_Container;

  function Standard_Copy_Container_Solutions_to_Target
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions;  
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Standard_Copy_Container_Solutions_to_Target.");
    end if;
    if Is_Null(sols) then
      return 6;
    else
      PHCpack_Operations.Store_Target_Solutions(sols);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("Standard_Copy_Container_Solutions_to_Target.");
      end if;
      return 6;
  end Standard_Copy_Container_Solutions_to_Target;

  function DoblDobl_Copy_Container_Solutions_to_Target
             ( vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Solutions;  
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("DoblDobl_Copy_Container_Solutions_to_Target.");
    end if;
    if Is_Null(sols) then
      return 256;
    else
      PHCpack_Operations.Store_Target_Solutions(sols);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("DoblDobl_Copy_Container_Solutions_to_Target.");
      end if;
      return 256;
  end DoblDobl_Copy_Container_Solutions_to_Target;

  function QuadDobl_Copy_Container_Solutions_to_Target
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Solutions;  
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("QuadDobl_Copy_Container_Solutions_to_Target.");
    end if;
    if Is_Null(sols) then
      return 266;
    else
      PHCpack_Operations.Store_Target_Solutions(sols);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("QuadDobl_Copy_Container_Solutions_to_Target.");
      end if;
      return 266;
  end QuadDobl_Copy_Container_Solutions_to_Target;

  function Multprec_Copy_Container_Solutions_to_Target
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Multprec_Complex_Solutions;  
    sols : constant Solution_List := Multprec_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Multprec_Copy_Container_Solutions_to_Target.");
    end if;
    if Is_Null(sols) then
      return 286;
    else
      PHCpack_Operations.Store_Target_Solutions(sols);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("Multprec_Copy_Container_Solutions_to_Target.");
      end if;
      return 286;
  end Multprec_Copy_Container_Solutions_to_Target;

  function Standard_Copy_Start_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions;
    sols : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Standard_Copy_Start_Solutions_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    if Is_Null(sols) then
      return 7;
    else
      Standard_Solutions_Container.Initialize(sols);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("Standard_Copy_Start_Solutions_to_Container.");
      end if;
      return 7;
  end Standard_Copy_Start_Solutions_to_Container;

  function DoblDobl_Copy_Start_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Solutions;
    sols : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("DoblDobl_Copy_Start_Solutions_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    if Is_Null(sols) then
      return 257;
    else
      DoblDobl_Solutions_Container.Initialize(sols);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("DoblDobl_Copy_Start_Solutions_to_Container.");
      end if;
      return 257;
  end DoblDobl_Copy_Start_Solutions_to_Container;

  function QuadDobl_Copy_Start_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Solutions;
    sols : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("QuadDobl_Copy_Start_Solutions_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    if Is_Null(sols) then
      return 267;
    else
      QuadDobl_Solutions_Container.Initialize(sols);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("QuadDobl_Copy_Start_Solutions_to_Container.");
      end if;
      return 267;
  end QuadDobl_Copy_Start_Solutions_to_Container;

  function Multprec_Copy_Start_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Multprec_Complex_Solutions;
    sols : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Multprec_Copy_Start_Solutions_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    if Is_Null(sols) then
      return 287;
    else
      Multprec_Solutions_Container.Initialize(sols);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("Multprec_Copy_Start_Solutions_to_Container.");
      end if;
      return 287;
  end Multprec_Copy_Start_Solutions_to_Container;

  function Standard_Copy_Container_Solutions_to_Start
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Standard_Copy_Container_Solutions_to_Start.");
    end if;
    if Is_Null(sols) then
      return 8;
    else
      PHCpack_Operations.Store_Start_Solutions(sols);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("Standard_Copy_Container_Solutions_to_Start.");
      end if;
      return 8;
  end Standard_Copy_Container_Solutions_to_Start;

  function DoblDobl_Copy_Container_Solutions_to_Start
             ( vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Solutions;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("DoblDobl_Copy_Container_Solutions_to_Start.");
    end if;
    if Is_Null(sols) then
      return 258;
    else
      PHCpack_Operations.Store_Start_Solutions(sols);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("DoblDobl_Copy_Container_Solutions_to_Start.");
      end if;
      return 258;
  end DoblDobl_Copy_Container_Solutions_to_Start;

  function QuadDobl_Copy_Container_Solutions_to_Start
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Solutions;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("QuadDobl_Copy_Container_Solutions_to_Start.");
    end if;
    if Is_Null(sols) then
      return 268;
    else
      PHCpack_Operations.Store_Start_Solutions(sols);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("QuadDobl_Copy_Container_Solutions_to_Start.");
      end if;
      return 268;
  end QuadDobl_Copy_Container_Solutions_to_Start;

  function Multprec_Copy_Container_Solutions_to_Start
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Multprec_Complex_Solutions;
    sols : constant Solution_List := Multprec_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Multprec_Copy_Container_Solutions_to_Start.");
    end if;
    if Is_Null(sols) then
      return 288;
    else
      PHCpack_Operations.Store_Start_Solutions(sols);
      return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("Multprec_Copy_Container_Solutions_to_Start.");
      end if;
      return 288;
  end Multprec_Copy_Container_Solutions_to_Start;

end Job_Containers;
