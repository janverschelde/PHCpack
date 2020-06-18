with text_io;                           use text_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Multprec_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Multprec_Complex_Solutions;
with Standard_PolySys_Container;
with Standard_LaurSys_Container;
with DoblDobl_PolySys_Container;
with DoblDobl_LaurSys_Container;
with QuadDobl_PolySys_Container;
with QuadDobl_LaurSys_Container;
with Multprec_PolySys_Container;
with Standard_Solutions_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_Solutions_Container;
with Multprec_Solutions_Container;
with PHCpack_Operations;

package body Job_Containers is

-- COPYING POLYNOMIAL SYSTEMS :

  function Standard_Target_Poly_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Standard_Target_Poly_System_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Target_System(lp);
    if lp = null
     then return 1;
     else Standard_PolySys_Container.Initialize(lp.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.Standard_Target_Poly_System_to_Container.");
      end if;
      return 1;
  end Standard_Target_Poly_System_to_Container;

  function DoblDobl_Target_Poly_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("DoblDobl_Target_Poly_System_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Target_System(lp);
    if lp = null
     then return 251;
     else DoblDobl_PolySys_Container.Initialize(lp.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.DoblDobl_Target_Poly_System_to_Container.");
      end if;
      return 251;
  end DoblDobl_Target_Poly_System_to_Container;

  function QuadDobl_Target_Poly_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("QuadDobl_Target_Poly_System_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Target_System(lp);
    if lp = null
     then return 261;
     else QuadDobl_PolySys_Container.Initialize(lp.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.DoblDobl_Target_Poly_System_to_Container.");
      end if;
      return 261;
  end QuadDobl_Target_Poly_System_to_Container;

  function Multprec_Target_Poly_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Multprec_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Multprec_Target_Poly_System_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Target_System(lp);
    if lp = null
     then return 281;
     else Multprec_PolySys_Container.Initialize(lp.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.Multprec_Target_Poly_System_to_Container.");
      end if;
      return 281;
  end Multprec_Target_Poly_System_to_Container;

  function Standard_Container_Poly_System_to_Target
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Standard_Container_Poly_System_to_Target.");
    end if;
    if lp = null
     then return 2;
     else PHCpack_Operations.Store_Target_System(lp.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.Standard_Container_Poly_System_to_Target.");
      end if;
      return 2;
  end Standard_Container_Poly_System_to_Target;

  function DoblDobl_Container_Poly_System_to_Target
             ( vrblvl : integer32 := 0 ) return integer32 is 

    use DoblDobl_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("DoblDobl_Container_Poly_System_to_Target.");
    end if;
    if lp = null
     then return 252;
     else PHCpack_Operations.Store_Target_System(lp.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.DoblDobl_Container_Poly_System_to_Target.");
      end if;
      return 252;
  end DoblDobl_Container_Poly_System_to_Target;

  function QuadDobl_Container_Poly_System_to_Target
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("QuadDobl_Container_Poly_System_to_Target.");
    end if;
    if lp = null
     then return 262;
     else PHCpack_Operations.Store_Target_System(lp.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.QuadDobl_Container_Poly_System_to_Target.");
      end if;
      return 262;
  end QuadDobl_Container_Poly_System_to_Target;

  function Multprec_Container_Poly_System_to_Target
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Multprec_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := Multprec_PolySys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Multprec_Container_Poly_System_to_Target.");
    end if;
    if lp = null
     then return 282;
     else PHCpack_Operations.Store_Target_System(lp.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.Multprec_Container_Poly_System_to_Target.");
      end if;
      return 282;
  end Multprec_Container_Poly_System_to_Target;

  function Standard_Start_Poly_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Standard_Start_Poly_System_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Start_System(lp);
    if lp = null
     then return 3;
     else Standard_PolySys_Container.Initialize(lp.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.Standard_Start_Poly_System_to_Container.");
      end if;
      return 3;
  end Standard_Start_Poly_System_to_Container;

  function DoblDobl_Start_Poly_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("DoblDobl_Start_Poly_System_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Start_System(lp);
    if lp = null
     then return 253;
     else DoblDobl_PolySys_Container.Initialize(lp.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.DoblDobl_Start_Poly_System_to_Container.");
      end if;
      return 253;
  end DoblDobl_Start_Poly_System_to_Container;

  function QuadDobl_Start_Poly_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("QuadDobl_Start_Poly_System_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Start_System(lp);
    if lp = null
     then return 263;
     else QuadDobl_PolySys_Container.Initialize(lp.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.QuadDobl_Start_Poly_System_to_Container.");
      end if;
      return 263;
  end QuadDobl_Start_Poly_System_to_Container;

  function Multprec_Start_Poly_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Multprec_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Multprec_Start_Poly_System_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Start_System(lp);
    if lp = null
     then return 283;
     else Multprec_PolySys_Container.Initialize(lp.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.Multprec_Start_Poly_System_to_Container.");
      end if;
      return 283;
  end Multprec_Start_Poly_System_to_Container;

  function Standard_Container_Poly_System_to_Start
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Standard_Container_Poly_System_to_Start.");
    end if;
    if lp = null
     then return 4;
     else PHCpack_Operations.Store_Start_System(lp.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.Standard_Container_Poly_System_to_Start.");
      end if;
      return 4;
  end Standard_Container_Poly_System_to_Start;

  function DoblDobl_Container_Poly_System_to_Start
             ( vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("DoblDobl_Container_Poly_System_to_Start.");
    end if;
    if lp = null
     then return 254;
     else PHCpack_Operations.Store_Start_System(lp.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.DoblDobl_Container_Poly_System_to_Start.");
      end if;
      return 254;
  end DoblDobl_Container_Poly_System_to_Start;

  function QuadDobl_Container_Poly_System_to_Start
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("QuadDobl_Container_Poly_System_to_Start.");
    end if;
    if lp = null
     then return 264;
     else PHCpack_Operations.Store_Start_System(lp.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.QuadDobl_Container_Poly_System_to_Start.");
      end if;
      return 264;
  end QuadDobl_Container_Poly_System_to_Start;

  function Multprec_Container_Poly_System_to_Start
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Multprec_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := Multprec_PolySys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Multprec_Container_Poly_System_to_Start.");
    end if;
    if lp = null
     then return 284;
     else PHCpack_Operations.Store_Start_System(lp.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.Multprec_Container_Poly_System_to_Start.");
      end if;
      return 284;
  end Multprec_Container_Poly_System_to_Start;

-- COPYING LAURENT POLYNOMIAL SYSTEMS :

  function Standard_Target_Laur_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Laur_Systems;
    q : Link_to_Laur_Sys;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Standard_Target_Laur_System_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Target_System(q);
    if q = null
     then return 786;
     else Standard_LaurSys_Container.Initialize(q.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.Standard_Target_Laur_System_to_Container.");
      end if;
      return 786;
  end Standard_Target_Laur_System_to_Container;

  function DoblDobl_Target_Laur_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Laur_Systems;
    q : Link_to_Laur_Sys;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("DoblDobl_Target_Laur_System_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Target_System(q);
    if q = null
     then return 787;
     else DoblDobl_LaurSys_Container.Initialize(q.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.DoblDobl_Target_Laur_System_to_Container.");
      end if;
      return 787;
  end DoblDobl_Target_Laur_System_to_Container;

  function QuadDobl_Target_Laur_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Laur_Systems;
    q : Link_to_Laur_Sys;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("QuadDobl_Target_Laur_System_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Target_System(q);
    if q = null
     then return 788;
     else QuadDobl_LaurSys_Container.Initialize(q.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.QuadDobl_Target_Laur_System_to_Container.");
      end if;
      return 788;
  end QuadDobl_Target_Laur_System_to_Container;

  function Standard_Container_Laur_System_to_Target
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Laur_Systems;
    q : constant Link_to_Laur_Sys := Standard_LaurSys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Standard_Container_Laur_System_to_Target.");
    end if;
    if q = null
     then return 780;
     else PHCpack_Operations.Store_Target_System(q.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.Standard_Container_Laur_System_to_Target.");
      end if;
      return 780;
  end Standard_Container_Laur_System_to_Target;

  function DoblDobl_Container_Laur_System_to_Target
             ( vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Laur_Systems;
    q : constant Link_to_Laur_Sys := DoblDobl_LaurSys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("DoblDobl_Container_Laur_System_to_Target.");
    end if;
    if q = null
     then return 781;
     else PHCpack_Operations.Store_Target_System(q.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.DoblDobl_Container_Laur_System_to_Target.");
      end if;
      return 781;
  end DoblDobl_Container_Laur_System_to_Target;

  function QuadDobl_Container_Laur_System_to_Target
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Laur_Systems;
    q : constant Link_to_Laur_Sys := QuadDobl_LaurSys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("QuadDobl_Container_Laur_System_to_Target.");
    end if;
    if q = null
     then return 782;
     else PHCpack_Operations.Store_Target_System(q.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.QuadDobl_Container_Laur_System_to_Target.");
      end if;
      return 782;
  end QuadDobl_Container_Laur_System_to_Target;

  function Standard_Start_Laur_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Laur_Systems;
    q : Link_to_Laur_Sys;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Standard_Start_Laur_System_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Start_System(q);
    if q = null
     then return 783;
     else Standard_LaurSys_Container.Initialize(q.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.Standard_Start_Laur_System_to_Container.");
      end if;
      return 783;
  end Standard_Start_Laur_System_to_Container;

  function DoblDobl_Start_Laur_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Laur_Systems;
    q : Link_to_Laur_Sys;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("DoblDobl_Start_Laur_System_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Start_System(q);
    if q = null
     then return 784;
     else DoblDobl_LaurSys_Container.Initialize(q.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.DoblDobl_Start_Laur_System_to_Container.");
      end if;
      return 784;
  end DoblDobl_Start_Laur_System_to_Container;

  function QuadDobl_Start_Laur_System_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Laur_Systems;
    q : Link_to_Laur_Sys;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("QuadDobl_Start_Laur_System_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Start_System(q);
    if q = null
     then return 785;
     else QuadDobl_LaurSys_Container.Initialize(q.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.QuadDobl_Start_Laur_System_to_Container.");
      end if;
      return 785;
  end QuadDobl_Start_Laur_System_to_Container;

  function Standard_Container_Laur_System_to_Start
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Laur_Systems;
    q : constant Link_to_Laur_Sys := Standard_LaurSys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Standard_Container_Laur_System_to_Start.");
    end if;
    if q = null
     then return 777;
     else PHCpack_Operations.Store_Start_System(q.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.Standard_Container_Laur_System_to_Start.");
      end if;
      return 777;
  end Standard_Container_Laur_System_to_Start;

  function DoblDobl_Container_Laur_System_to_Start
             ( vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Laur_Systems;
    q : constant Link_to_Laur_Sys := DoblDobl_LaurSys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("DoblDobl_Container_Laur_System_to_Start.");
    end if;
    if q = null
     then return 778;
     else PHCpack_Operations.Store_Start_System(q.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.DoblDobl_Container_Laur_System_to_Start.");
      end if;
      return 778;
  end DoblDobl_Container_Laur_System_to_Start;

  function QuadDobl_Container_Laur_System_to_Start
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Laur_Systems;
    q : constant Link_to_Laur_Sys := QuadDobl_LaurSys_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("QuadDobl_Container_Laur_System_to_Start.");
    end if;
    if q = null
     then return 779;
     else PHCpack_Operations.Store_Start_System(q.all); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at ");
        put_line("job_containers.QuadDobl_Container_Laur_System_to_Start.");
      end if;
      return 779;
  end QuadDobl_Container_Laur_System_to_Start;

-- COPYING SOLUTIONS :

  function Standard_Target_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions;
    sols : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Standard_Target_Solutions_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Target_Solutions(sols);
    if Is_Null(sols)
     then return 5;
     else Standard_Solutions_Container.Initialize(sols); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("Standard_Target_Solutions_to_Container.");
      end if;
      put_line("Exception 5 at copy of target solutions to container.");
      return 5;
  end Standard_Target_Solutions_to_Container;

  function DoblDobl_Target_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Solutions;
    sols : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("DoblDobl_Target_Solutions_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Target_Solutions(sols);
    if Is_Null(sols)
     then return 255;
     else DoblDobl_Solutions_Container.Initialize(sols); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("DoblDobl_Target_Solutions_to_Container.");
      end if;
      return 255;
  end DoblDobl_Target_Solutions_to_Container;

  function QuadDobl_Target_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Solutions;
    sols : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("QuadDobl_Target_Solutions_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Target_Solutions(sols);
    if Is_Null(sols)
     then return 265;
     else QuadDobl_Solutions_Container.Initialize(sols); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("QuadDobl_Target_Solutions_to_Container.");
      end if;
      return 265;
  end QuadDobl_Target_Solutions_to_Container;

  function Multprec_Target_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Multprec_Complex_Solutions;
    sols : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Multprec_Target_Solutions_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Target_Solutions(sols);
    if Is_Null(sols)
     then return 285;
     else Multprec_Solutions_Container.Initialize(sols); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("Multprec_Target_Solutions_to_Container.");
      end if;
      return 285;
  end Multprec_Target_Solutions_to_Container;

  function Standard_Container_Solutions_to_Target
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions;  
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Standard_Container_Solutions_to_Target.");
    end if;
    if Is_Null(sols)
     then return 6;
     else PHCpack_Operations.Store_Target_Solutions(sols); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("Standard_Container_Solutions_to_Target.");
      end if;
      return 6;
  end Standard_Container_Solutions_to_Target;

  function DoblDobl_Container_Solutions_to_Target
             ( vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Solutions;  
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("DoblDobl_Container_Solutions_to_Target.");
    end if;
    if Is_Null(sols)
     then return 256;
     else PHCpack_Operations.Store_Target_Solutions(sols); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("DoblDobl_Container_Solutions_to_Target.");
      end if;
      return 256;
  end DoblDobl_Container_Solutions_to_Target;

  function QuadDobl_Container_Solutions_to_Target
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Solutions;  
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("QuadDobl_Container_Solutions_to_Target.");
    end if;
    if Is_Null(sols)
     then return 266;
     else PHCpack_Operations.Store_Target_Solutions(sols); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("QuadDobl_Container_Solutions_to_Target.");
      end if;
      return 266;
  end QuadDobl_Container_Solutions_to_Target;

  function Multprec_Container_Solutions_to_Target
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Multprec_Complex_Solutions;  
    sols : constant Solution_List := Multprec_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Multprec_Container_Solutions_to_Target.");
    end if;
    if Is_Null(sols)
     then return 286;
     else PHCpack_Operations.Store_Target_Solutions(sols); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("Multprec_Container_Solutions_to_Target.");
      end if;
      return 286;
  end Multprec_Container_Solutions_to_Target;

  function Standard_Start_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions;
    sols : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Standard_Start_Solutions_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    if Is_Null(sols)
     then return 7;
     else Standard_Solutions_Container.Initialize(sols); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("Standard_Start_Solutions_to_Container.");
      end if;
      return 7;
  end Standard_Start_Solutions_to_Container;

  function DoblDobl_Start_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Solutions;
    sols : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("DoblDobl_Start_Solutions_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    if Is_Null(sols)
     then return 257;
     else DoblDobl_Solutions_Container.Initialize(sols); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("DoblDobl_Start_Solutions_to_Container.");
      end if;
      return 257;
  end DoblDobl_Start_Solutions_to_Container;

  function QuadDobl_Start_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Solutions;
    sols : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("QuadDobl_Start_Solutions_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    if Is_Null(sols)
     then return 267;
     else QuadDobl_Solutions_Container.Initialize(sols); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("QuadDobl_Start_Solutions_to_Container.");
      end if;
      return 267;
  end QuadDobl_Start_Solutions_to_Container;

  function Multprec_Start_Solutions_to_Container
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Multprec_Complex_Solutions;
    sols : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Multprec_Start_Solutions_to_Container.");
    end if;
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    if Is_Null(sols)
     then return 287;
     else Multprec_Solutions_Container.Initialize(sols); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("Multprec_Start_Solutions_to_Container.");
      end if;
      return 287;
  end Multprec_Start_Solutions_to_Container;

  function Standard_Container_Solutions_to_Start
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Solutions;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Standard_Container_Solutions_to_Start.");
    end if;
    if Is_Null(sols)
     then return 8;
     else PHCpack_Operations.Store_Start_Solutions(sols); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("Standard_Container_Solutions_to_Start.");
      end if;
      return 8;
  end Standard_Container_Solutions_to_Start;

  function DoblDobl_Container_Solutions_to_Start
             ( vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Solutions;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("DoblDobl_Container_Solutions_to_Start.");
    end if;
    if Is_Null(sols)
     then return 258;
     else PHCpack_Operations.Store_Start_Solutions(sols); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("DoblDobl_Container_Solutions_to_Start.");
      end if;
      return 258;
  end DoblDobl_Container_Solutions_to_Start;

  function QuadDobl_Container_Solutions_to_Start
             ( vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Solutions;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("QuadDobl_Container_Solutions_to_Start.");
    end if;
    if Is_Null(sols)
     then return 268;
     else PHCpack_Operations.Store_Start_Solutions(sols); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("QuadDobl_Container_Solutions_to_Start.");
      end if;
      return 268;
  end QuadDobl_Container_Solutions_to_Start;

  function Multprec_Container_Solutions_to_Start
             ( vrblvl : integer32 := 0 ) return integer32 is

    use Multprec_Complex_Solutions;
    sols : constant Solution_List := Multprec_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in job_containers.");
      put_line("Multprec_Container_Solutions_to_Start.");
    end if;
    if Is_Null(sols)
     then return 288;
     else PHCpack_Operations.Store_Start_Solutions(sols); return 0;
    end if;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised at job_containers.");
        put_line("Multprec_Container_Solutions_to_Start.");
      end if;
      return 288;
  end Multprec_Container_Solutions_to_Start;

end Job_Containers;
