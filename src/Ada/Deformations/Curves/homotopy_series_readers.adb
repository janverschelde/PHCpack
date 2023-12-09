with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Random_Numbers;
with DoblDobl_Random_Numbers;
with TripDobl_Random_Numbers;
with QuadDobl_Random_Numbers;
with PentDobl_Random_Numbers;
with OctoDobl_Random_Numbers;
with DecaDobl_Random_Numbers;
with HexaDobl_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Homotopy;
with Standard_Coefficient_Homotopy;
with DoblDobl_Homotopy;
with DoblDobl_Coefficient_Homotopy;
with TripDobl_Homotopy;
with TripDobl_Coefficient_Homotopy;
with QuadDobl_Homotopy;
with QuadDobl_Coefficient_Homotopy;
with PentDobl_Homotopy;
with PentDobl_Coefficient_Homotopy;
with OctoDobl_Homotopy;
with OctoDobl_Coefficient_Homotopy;
with DecaDobl_Homotopy;
with DecaDobl_Coefficient_Homotopy;
with HexaDobl_Homotopy;
with HexaDobl_Coefficient_Homotopy;
with Artificial_Parameter_Homotopy_io;
with Projective_Transformations;         use Projective_Transformations;
with Multi_Projective_Transformations;   use Multi_Projective_Transformations;
with Homogenization;                     use Homogenization;
with Standard_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_Systems;
with PentDobl_CSeries_Poly_Systems;
with OctoDobl_CSeries_Poly_Systems;
with DecaDobl_CSeries_Poly_Systems;
with HexaDobl_CSeries_Poly_Systems;
with Standard_Parameter_Systems;
with DoblDobl_Parameter_Systems;
with TripDobl_Parameter_Systems;
with QuadDobl_Parameter_Systems;
with PentDobl_Parameter_Systems;
with OctoDobl_Parameter_Systems;
with DecaDobl_Parameter_Systems;
with HexaDobl_Parameter_Systems;
with Series_and_Homotopies;
with Series_and_Predictors;
with Jacobian_Rabinowitsch_Trick;        use Jacobian_Rabinowitsch_Trick;

package body Homotopy_Series_Readers is

  procedure Standard_Projective_Transformation
              ( target,start
                 : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    ptar : Standard_Complex_Poly_Systems.Poly_Sys(target'first..target'last+1);
    pstr : Standard_Complex_Poly_Systems.Poly_Sys(start'first..start'last+1);

  begin
    Projective_Transformation(target.all);
    Projective_Transformation(start.all);
    ptar := Add_Random_Hyperplanes(target.all,1,false);
    pstr := Add_Standard_Hyperplanes(start.all,1);
    Standard_Complex_Poly_Systems.Clear(target);
    target := new Standard_Complex_Poly_Systems.Poly_Sys'(ptar);
    Standard_Complex_Poly_Systems.Clear(start);
    start := new Standard_Complex_Poly_Systems.Poly_Sys'(pstr);
  end Standard_Projective_Transformation;

  procedure DoblDobl_Projective_Transformation
              ( target,start
                 : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    ptar : DoblDobl_Complex_Poly_Systems.Poly_Sys(target'first..target'last+1);
    pstr : DoblDobl_Complex_Poly_Systems.Poly_Sys(start'first..start'last+1);

  begin
    Projective_Transformation(target.all);
    Projective_Transformation(start.all);
    ptar := Add_Random_Hyperplanes(target.all,1,false);
    pstr := Add_Standard_Hyperplanes(start.all,1);
    DoblDobl_Complex_Poly_Systems.Clear(target);
    target := new DoblDobl_Complex_Poly_Systems.Poly_Sys'(ptar);
    DoblDobl_Complex_Poly_Systems.Clear(start);
    start := new DoblDobl_Complex_Poly_Systems.Poly_Sys'(pstr);
  end DoblDobl_Projective_Transformation;

  procedure TripDobl_Projective_Transformation
              ( target,start
                 : in out TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    ptar : TripDobl_Complex_Poly_Systems.Poly_Sys(target'first..target'last+1);
    pstr : TripDobl_Complex_Poly_Systems.Poly_Sys(start'first..start'last+1);

  begin
    Projective_Transformation(target.all);
    Projective_Transformation(start.all);
    ptar := Add_Random_Hyperplanes(target.all,1,false);
    pstr := Add_Standard_Hyperplanes(start.all,1);
    TripDobl_Complex_Poly_Systems.Clear(target);
    target := new TripDobl_Complex_Poly_Systems.Poly_Sys'(ptar);
    TripDobl_Complex_Poly_Systems.Clear(start);
    start := new TripDobl_Complex_Poly_Systems.Poly_Sys'(pstr);
  end TripDobl_Projective_Transformation;

  procedure QuadDobl_Projective_Transformation
              ( target,start
                 : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    ptar : QuadDobl_Complex_Poly_Systems.Poly_Sys(target'first..target'last+1);
    pstr : QuadDobl_Complex_Poly_Systems.Poly_Sys(start'first..start'last+1);

  begin
    Projective_Transformation(target.all);
    Projective_Transformation(start.all);
    ptar := Add_Random_Hyperplanes(target.all,1,false);
    pstr := Add_Standard_Hyperplanes(start.all,1);
    QuadDobl_Complex_Poly_Systems.Clear(target);
    target := new QuadDobl_Complex_Poly_Systems.Poly_Sys'(ptar);
    QuadDobl_Complex_Poly_Systems.Clear(start);
    start := new QuadDobl_Complex_Poly_Systems.Poly_Sys'(pstr);
  end QuadDobl_Projective_Transformation;

  procedure PentDobl_Projective_Transformation
              ( target,start
                 : in out PentDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    ptar : PentDobl_Complex_Poly_Systems.Poly_Sys(target'first..target'last+1);
    pstr : PentDobl_Complex_Poly_Systems.Poly_Sys(start'first..start'last+1);

  begin
    Projective_Transformation(target.all);
    Projective_Transformation(start.all);
    ptar := Add_Random_Hyperplanes(target.all,1,false);
    pstr := Add_Standard_Hyperplanes(start.all,1);
    PentDobl_Complex_Poly_Systems.Clear(target);
    target := new PentDobl_Complex_Poly_Systems.Poly_Sys'(ptar);
    PentDobl_Complex_Poly_Systems.Clear(start);
    start := new PentDobl_Complex_Poly_Systems.Poly_Sys'(pstr);
  end PentDobl_Projective_Transformation;

  procedure OctoDobl_Projective_Transformation
              ( target,start
                 : in out OctoDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    ptar : OctoDobl_Complex_Poly_Systems.Poly_Sys(target'first..target'last+1);
    pstr : OctoDobl_Complex_Poly_Systems.Poly_Sys(start'first..start'last+1);

  begin
    Projective_Transformation(target.all);
    Projective_Transformation(start.all);
    ptar := Add_Random_Hyperplanes(target.all,1,false);
    pstr := Add_Standard_Hyperplanes(start.all,1);
    OctoDobl_Complex_Poly_Systems.Clear(target);
    target := new OctoDobl_Complex_Poly_Systems.Poly_Sys'(ptar);
    OctoDobl_Complex_Poly_Systems.Clear(start);
    start := new OctoDobl_Complex_Poly_Systems.Poly_Sys'(pstr);
  end OctoDobl_Projective_Transformation;

  procedure DecaDobl_Projective_Transformation
              ( target,start
                 : in out DecaDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    ptar : DecaDobl_Complex_Poly_Systems.Poly_Sys(target'first..target'last+1);
    pstr : DecaDobl_Complex_Poly_Systems.Poly_Sys(start'first..start'last+1);

  begin
    Projective_Transformation(target.all);
    Projective_Transformation(start.all);
    ptar := Add_Random_Hyperplanes(target.all,1,false);
    pstr := Add_Standard_Hyperplanes(start.all,1);
    DecaDobl_Complex_Poly_Systems.Clear(target);
    target := new DecaDobl_Complex_Poly_Systems.Poly_Sys'(ptar);
    DecaDobl_Complex_Poly_Systems.Clear(start);
    start := new DecaDobl_Complex_Poly_Systems.Poly_Sys'(pstr);
  end DecaDobl_Projective_Transformation;

  procedure HexaDobl_Projective_Transformation
              ( target,start
                 : in out HexaDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    ptar : HexaDobl_Complex_Poly_Systems.Poly_Sys(target'first..target'last+1);
    pstr : HexaDobl_Complex_Poly_Systems.Poly_Sys(start'first..start'last+1);

  begin
    Projective_Transformation(target.all);
    Projective_Transformation(start.all);
    ptar := Add_Random_Hyperplanes(target.all,1,false);
    pstr := Add_Standard_Hyperplanes(start.all,1);
    HexaDobl_Complex_Poly_Systems.Clear(target);
    target := new HexaDobl_Complex_Poly_Systems.Poly_Sys'(ptar);
    HexaDobl_Complex_Poly_Systems.Clear(start);
    start := new HexaDobl_Complex_Poly_Systems.Poly_Sys'(pstr);
  end HexaDobl_Projective_Transformation;

  procedure Standard_Projective_Transformation
              ( target : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                start : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List ) is
  begin
    Standard_Projective_Transformation(target,start);
    Projective_Transformation(sols);
  end Standard_Projective_Transformation;

  procedure DoblDobl_Projective_Transformation
              ( target : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                start : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is
  begin
    DoblDobl_Projective_Transformation(target,start);
    Projective_Transformation(sols);
  end DoblDobl_Projective_Transformation;

  procedure TripDobl_Projective_Transformation
              ( target : in out TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                start : in out TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out TripDobl_Complex_Solutions.Solution_List ) is
  begin
    TripDobl_Projective_Transformation(target,start);
    Projective_Transformation(sols);
  end TripDobl_Projective_Transformation;

  procedure QuadDobl_Projective_Transformation
              ( target : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                start : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is
  begin
    QuadDobl_Projective_Transformation(target,start);
    Projective_Transformation(sols);
  end QuadDobl_Projective_Transformation;

  procedure PentDobl_Projective_Transformation
              ( target : in out PentDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                start : in out PentDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out PentDobl_Complex_Solutions.Solution_List ) is
  begin
    PentDobl_Projective_Transformation(target,start);
    Projective_Transformation(sols);
  end PentDobl_Projective_Transformation;

  procedure OctoDobl_Projective_Transformation
              ( target : in out OctoDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                start : in out OctoDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out OctoDobl_Complex_Solutions.Solution_List ) is
  begin
    OctoDobl_Projective_Transformation(target,start);
    Projective_Transformation(sols);
  end OctoDobl_Projective_Transformation;

  procedure DecaDobl_Projective_Transformation
              ( target : in out DecaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                start : in out DecaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out DecaDobl_Complex_Solutions.Solution_List ) is
  begin
    DecaDobl_Projective_Transformation(target,start);
    Projective_Transformation(sols);
  end DecaDobl_Projective_Transformation;

  procedure HexaDobl_Projective_Transformation
              ( target : in out HexaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                start : in out HexaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out HexaDobl_Complex_Solutions.Solution_List ) is
  begin
    HexaDobl_Projective_Transformation(target,start);
    Projective_Transformation(sols);
  end HexaDobl_Projective_Transformation;

  procedure Standard_Multi_Projective_Transformation
              ( target,start
                  : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                m : in natural32; z : in Partition ) is

    p : constant Standard_Complex_Poly_Systems.Poly_Sys
      := Multi_Projective_Transformation(target.all,m,z);
    q : constant Standard_Complex_Poly_Systems.Poly_Sys
      := Multi_Projective_Transformation(start.all,m,z,true);

  begin
    Standard_Complex_Poly_Systems.Clear(target);
    target := new Standard_Complex_Poly_Systems.Poly_Sys'(p);
    Standard_Complex_Poly_Systems.Clear(start);
    start := new Standard_Complex_Poly_Systems.Poly_Sys'(q);
  end Standard_Multi_Projective_Transformation;

  procedure DoblDobl_Multi_Projective_Transformation
              ( target,start
                  : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                m : in natural32; z : in Partition ) is

    p : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
      := Multi_Projective_Transformation(target.all,m,z);
    q : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
      := Multi_Projective_Transformation(start.all,m,z,true);

  begin
    DoblDobl_Complex_Poly_Systems.Clear(target);
    target := new DoblDobl_Complex_Poly_Systems.Poly_Sys'(p);
    DoblDobl_Complex_Poly_Systems.Clear(start);
    start := new DoblDobl_Complex_Poly_Systems.Poly_Sys'(q);
  end DoblDobl_Multi_Projective_Transformation;

  procedure TripDobl_Multi_Projective_Transformation
              ( target,start
                  : in out TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                m : in natural32; z : in Partition ) is

    p : constant TripDobl_Complex_Poly_Systems.Poly_Sys
      := Multi_Projective_Transformation(target.all,m,z);
    q : constant TripDobl_Complex_Poly_Systems.Poly_Sys
      := Multi_Projective_Transformation(start.all,m,z,true);

  begin
    TripDobl_Complex_Poly_Systems.Clear(target);
    target := new TripDobl_Complex_Poly_Systems.Poly_Sys'(p);
    TripDobl_Complex_Poly_Systems.Clear(start);
    start := new TripDobl_Complex_Poly_Systems.Poly_Sys'(q);
  end TripDobl_Multi_Projective_Transformation;

  procedure QuadDobl_Multi_Projective_Transformation
              ( target,start
                  : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                m : in natural32; z : in Partition ) is

    p : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
      := Multi_Projective_Transformation(target.all,m,z);
    q : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
      := Multi_Projective_Transformation(start.all,m,z,true);

  begin
    QuadDobl_Complex_Poly_Systems.Clear(target);
    target := new QuadDobl_Complex_Poly_Systems.Poly_Sys'(p);
    QuadDobl_Complex_Poly_Systems.Clear(start);
    start := new QuadDobl_Complex_Poly_Systems.Poly_Sys'(q);
  end QuadDobl_Multi_Projective_Transformation;

  procedure PentDobl_Multi_Projective_Transformation
              ( target,start
                  : in out PentDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                m : in natural32; z : in Partition ) is

    p : constant PentDobl_Complex_Poly_Systems.Poly_Sys
      := Multi_Projective_Transformation(target.all,m,z);
    q : constant PentDobl_Complex_Poly_Systems.Poly_Sys
      := Multi_Projective_Transformation(start.all,m,z,true);

  begin
    PentDobl_Complex_Poly_Systems.Clear(target);
    target := new PentDobl_Complex_Poly_Systems.Poly_Sys'(p);
    PentDobl_Complex_Poly_Systems.Clear(start);
    start := new PentDobl_Complex_Poly_Systems.Poly_Sys'(q);
  end PentDobl_Multi_Projective_Transformation;

  procedure OctoDobl_Multi_Projective_Transformation
              ( target,start
                  : in out OctoDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                m : in natural32; z : in Partition ) is

    p : constant OctoDobl_Complex_Poly_Systems.Poly_Sys
      := Multi_Projective_Transformation(target.all,m,z);
    q : constant OctoDobl_Complex_Poly_Systems.Poly_Sys
      := Multi_Projective_Transformation(start.all,m,z,true);

  begin
    OctoDobl_Complex_Poly_Systems.Clear(target);
    target := new OctoDobl_Complex_Poly_Systems.Poly_Sys'(p);
    OctoDobl_Complex_Poly_Systems.Clear(start);
    start := new OctoDobl_Complex_Poly_Systems.Poly_Sys'(q);
  end OctoDobl_Multi_Projective_Transformation;

  procedure DecaDobl_Multi_Projective_Transformation
              ( target,start
                  : in out DecaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                m : in natural32; z : in Partition ) is

    p : constant DecaDobl_Complex_Poly_Systems.Poly_Sys
      := Multi_Projective_Transformation(target.all,m,z);
    q : constant DecaDobl_Complex_Poly_Systems.Poly_Sys
      := Multi_Projective_Transformation(start.all,m,z,true);

  begin
    DecaDobl_Complex_Poly_Systems.Clear(target);
    target := new DecaDobl_Complex_Poly_Systems.Poly_Sys'(p);
    DecaDobl_Complex_Poly_Systems.Clear(start);
    start := new DecaDobl_Complex_Poly_Systems.Poly_Sys'(q);
  end DecaDobl_Multi_Projective_Transformation;

  procedure HexaDobl_Multi_Projective_Transformation
              ( target,start
                  : in out HexaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                m : in natural32; z : in Partition ) is

    p : constant HexaDobl_Complex_Poly_Systems.Poly_Sys
      := Multi_Projective_Transformation(target.all,m,z);
    q : constant HexaDobl_Complex_Poly_Systems.Poly_Sys
      := Multi_Projective_Transformation(start.all,m,z,true);

  begin
    HexaDobl_Complex_Poly_Systems.Clear(target);
    target := new HexaDobl_Complex_Poly_Systems.Poly_Sys'(p);
    HexaDobl_Complex_Poly_Systems.Clear(start);
    start := new HexaDobl_Complex_Poly_Systems.Poly_Sys'(q);
  end HexaDobl_Multi_Projective_Transformation;

  procedure Standard_Multi_Projective_Transformation
              ( target : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                start : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                m : in natural32; z : in Partition ) is
  begin
    Standard_Multi_Projective_Transformation(target,start,m,z);
    Add_Ones(sols,m);
  end Standard_Multi_Projective_Transformation;

  procedure DoblDobl_Multi_Projective_Transformation
              ( target : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                start : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                m : in natural32; z : in Partition ) is
  begin
    DoblDobl_Multi_Projective_Transformation(target,start,m,z);
    Add_Ones(sols,m);
  end DoblDobl_Multi_Projective_Transformation;

  procedure TripDobl_Multi_Projective_Transformation
              ( target : in out TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                start : in out TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out TripDobl_Complex_Solutions.Solution_List;
                m : in natural32; z : in Partition ) is
  begin
    TripDobl_Multi_Projective_Transformation(target,start,m,z);
    Add_Ones(sols,m);
  end TripDobl_Multi_Projective_Transformation;

  procedure QuadDobl_Multi_Projective_Transformation
              ( target : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                start : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                m : in natural32; z : in Partition ) is
  begin
    QuadDobl_Multi_Projective_Transformation(target,start,m,z);
    Add_Ones(sols,m);
  end QuadDobl_Multi_Projective_Transformation;

  procedure PentDobl_Multi_Projective_Transformation
              ( target : in out PentDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                start : in out PentDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out PentDobl_Complex_Solutions.Solution_List;
                m : in natural32; z : in Partition ) is
  begin
    PentDobl_Multi_Projective_Transformation(target,start,m,z);
    Add_Ones(sols,m);
  end PentDobl_Multi_Projective_Transformation;

  procedure OctoDobl_Multi_Projective_Transformation
              ( target : in out OctoDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                start : in out OctoDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out OctoDobl_Complex_Solutions.Solution_List;
                m : in natural32; z : in Partition ) is
  begin
    OctoDobl_Multi_Projective_Transformation(target,start,m,z);
    Add_Ones(sols,m);
  end OctoDobl_Multi_Projective_Transformation;

  procedure DecaDobl_Multi_Projective_Transformation
              ( target : in out DecaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                start : in out DecaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out DecaDobl_Complex_Solutions.Solution_List;
                m : in natural32; z : in Partition ) is
  begin
    DecaDobl_Multi_Projective_Transformation(target,start,m,z);
    Add_Ones(sols,m);
  end DecaDobl_Multi_Projective_Transformation;

  procedure HexaDobl_Multi_Projective_Transformation
              ( target : in out HexaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                start : in out HexaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out HexaDobl_Complex_Solutions.Solution_List;
                m : in natural32; z : in Partition ) is
  begin
    HexaDobl_Multi_Projective_Transformation(target,start,m,z);
    Add_Ones(sols,m);
  end HexaDobl_Multi_Projective_Transformation;

  procedure Standard_Reader
              ( nbequ : out integer32;
                sols : out Standard_Complex_Solutions.Solution_List;
               -- tpow : in natural32;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                homcrd,rabin : in boolean := false ) is

    target,start : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    ans : character;

  begin
    Artificial_Parameter_Homotopy_io.get(target,start,sols);
    if not rabin then
      if homcrd then
        Standard_Projective_Transformation(target,start,sols);
      end if;
      nbequ := target'last;
      if not homcrd then
        Standard_Homotopy.Create(target.all,start.all,1,gamma); -- tpow,gamma);
      else -- set the value for tpow to one if homcrd
        Standard_Homotopy.Create(target.all,start.all,1,gamma);
        Standard_Coefficient_Homotopy.Create(start.all,target.all,1,gamma);
      end if;
    else
      new_line;
      put("Apply Rabinowitsch trick to put singularities at infinity ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans /= 'y' then
        nbequ := target'last;
        Standard_Homotopy.Create(target.all,start.all,1,gamma); -- tpow,gamma);
      else
        declare
          jrbtarget : constant Standard_Complex_Poly_Systems.Poly_Sys
                    := Jacobian_Rabinowitsch(target.all);
          jrbstart : constant Standard_Complex_Poly_Systems.Poly_Sys
                   := Jacobian_Rabinowitsch(start.all);
          jrbsols : constant Standard_Complex_Solutions.Solution_List
                  := Jacobian_Rabinowitsch(sols);
        begin
          nbequ := jrbtarget'last;
          Standard_Homotopy.Create(jrbtarget,jrbstart,1,gamma); -- tpow,gamma);
          Standard_Complex_Solutions.Deep_Clear(sols);
          sols := jrbsols;
        end;
      end if;
    end if;
  end Standard_Reader;

  procedure DoblDobl_Reader
              ( nbequ : out integer32;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
               -- tpow : in natural32;
                gamma : in DoblDobl_Complex_Numbers.Complex_Number;
                homcrd,rabin : in boolean := false ) is

    target,start : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    ans : character;

  begin
    Artificial_Parameter_Homotopy_io.get(target,start,sols);
    if not rabin then
      if homcrd
       then DoblDobl_Projective_Transformation(target,start,sols);
      end if;
      nbequ := target'last;
      if not homcrd then
        DoblDobl_Homotopy.Create(target.all,start.all,1,gamma); -- tpow,gamma);
      else -- set tpow to one for homogeneous coordinates
        DoblDobl_Homotopy.Create(target.all,start.all,1,gamma);
        DoblDobl_Coefficient_Homotopy.Create(start.all,target.all,1,gamma);
      end if;
    else
      new_line;
      put("Apply Rabinowitsch trick to put singularities at infinity ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans /= 'y' then
        nbequ := target'last;
        DoblDobl_Homotopy.Create(target.all,start.all,1,gamma); -- tpow,gamma);
      else
        declare
          jrbtarget : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
                    := Jacobian_Rabinowitsch(target.all);
          jrbstart : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
                   := Jacobian_Rabinowitsch(start.all);
          jrbsols : constant DoblDobl_Complex_Solutions.Solution_List
                  := Jacobian_Rabinowitsch(sols);
        begin
          nbequ := jrbtarget'last;
          DoblDobl_Homotopy.Create(jrbtarget,jrbstart,1,gamma); -- tpow,gamma);
          DoblDobl_Complex_Solutions.Deep_Clear(sols);
          sols := jrbsols;
        end;
      end if;
    end if;
  end DoblDobl_Reader;

  procedure TripDobl_Reader
              ( nbequ : out integer32;
                sols : out TripDobl_Complex_Solutions.Solution_List;
               -- tpow : in natural32;
                gamma : in TripDobl_Complex_Numbers.Complex_Number;
                homcrd,rabin : in boolean := false ) is

    target,start : TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    Artificial_Parameter_Homotopy_io.get(target,start,sols);
    if not rabin then
      if homcrd
       then TripDobl_Projective_Transformation(target,start,sols);
      end if;
      nbequ := target'last;
      if not homcrd then
        TripDobl_Homotopy.Create(target.all,start.all,1,gamma); -- tpow,gamma);
      else -- set tpow to one for homogeneous coordinates
        TripDobl_Homotopy.Create(target.all,start.all,1,gamma);
        TripDobl_Coefficient_Homotopy.Create(start.all,target.all,1,gamma);
      end if;
    else
      nbequ := target'last;
      TripDobl_Homotopy.Create(target.all,start.all,1,gamma); -- tpow,gamma);
    end if;
  end TripDobl_Reader;

  procedure QuadDobl_Reader
              ( nbequ : out integer32;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
               -- tpow : in natural32;
                gamma : in QuadDobl_Complex_Numbers.Complex_Number;
                homcrd,rabin : in boolean := false ) is

    target,start : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    ans : character;

  begin
    Artificial_Parameter_Homotopy_io.get(target,start,sols);
    if not rabin then
      if homcrd
       then QuadDobl_Projective_Transformation(target,start,sols);
      end if;
      nbequ := target'last;
      if not homcrd then
        QuadDobl_Homotopy.Create(target.all,start.all,1,gamma); -- tpow,gamma);
      else -- set tpow to one for homogeneous coordinates
        QuadDobl_Homotopy.Create(target.all,start.all,1,gamma);
        QuadDobl_Coefficient_Homotopy.Create(start.all,target.all,1,gamma);
      end if;
    else
      new_line;
      put("Apply Rabinowitsch trick to put singularities at infinity ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans /= 'y' then
        nbequ := target'last;
        QuadDobl_Homotopy.Create(target.all,start.all,1,gamma); -- tpow,gamma);
      else
        declare
          jrbtarget : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
                    := Jacobian_Rabinowitsch(target.all);
          jrbstart : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
                   := Jacobian_Rabinowitsch(start.all);
          jrbsols : constant QuadDobl_Complex_Solutions.Solution_List
                  := Jacobian_Rabinowitsch(sols);
        begin
          nbequ := jrbtarget'last;
          QuadDobl_Homotopy.Create(jrbtarget,jrbstart,1,gamma); -- tpow,gamma);
          QuadDobl_Complex_Solutions.Deep_Clear(sols);
          sols := jrbsols;
        end;
      end if;
    end if;
  end QuadDobl_Reader;

  procedure PentDobl_Reader
              ( nbequ : out integer32;
                sols : out PentDobl_Complex_Solutions.Solution_List;
               -- tpow : in natural32;
                gamma : in PentDobl_Complex_Numbers.Complex_Number;
                homcrd,rabin : in boolean := false ) is

    target,start : PentDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    Artificial_Parameter_Homotopy_io.get(target,start,sols);
    if not rabin then
      if homcrd
       then PentDobl_Projective_Transformation(target,start,sols);
      end if;
      nbequ := target'last;
      if not homcrd then
        PentDobl_Homotopy.Create(target.all,start.all,1,gamma); -- tpow,gamma);
      else -- set tpow to one for homogeneous coordinates
        PentDobl_Homotopy.Create(target.all,start.all,1,gamma);
        PentDobl_Coefficient_Homotopy.Create(start.all,target.all,1,gamma);
      end if;
    else
      nbequ := target'last;
      PentDobl_Homotopy.Create(target.all,start.all,1,gamma); -- tpow,gamma);
    end if;
  end PentDobl_Reader;

  procedure OctoDobl_Reader
              ( nbequ : out integer32;
                sols : out OctoDobl_Complex_Solutions.Solution_List;
               -- tpow : in natural32;
                gamma : in OctoDobl_Complex_Numbers.Complex_Number;
                homcrd,rabin : in boolean := false ) is

    target,start : OctoDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    Artificial_Parameter_Homotopy_io.get(target,start,sols);
    if not rabin then
      if homcrd
       then OctoDobl_Projective_Transformation(target,start,sols);
      end if;
      nbequ := target'last;
      if not homcrd then
        OctoDobl_Homotopy.Create(target.all,start.all,1,gamma); -- tpow,gamma);
      else -- set tpow to one for homogeneous coordinates
        OctoDobl_Homotopy.Create(target.all,start.all,1,gamma);
        OctoDobl_Coefficient_Homotopy.Create(start.all,target.all,1,gamma);
      end if;
    else
      nbequ := target'last;
      OctoDobl_Homotopy.Create(target.all,start.all,1,gamma); -- tpow,gamma);
    end if;
  end OctoDobl_Reader;

  procedure DecaDobl_Reader
              ( nbequ : out integer32;
                sols : out DecaDobl_Complex_Solutions.Solution_List;
               -- tpow : in natural32;
                gamma : in DecaDobl_Complex_Numbers.Complex_Number;
                homcrd,rabin : in boolean := false ) is

    target,start : DecaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    Artificial_Parameter_Homotopy_io.get(target,start,sols);
    if not rabin then
      if homcrd
       then DecaDobl_Projective_Transformation(target,start,sols);
      end if;
      nbequ := target'last;
      if not homcrd then
        DecaDobl_Homotopy.Create(target.all,start.all,1,gamma); -- tpow,gamma);
      else -- set tpow to one for homogeneous coordinates
        DecaDobl_Homotopy.Create(target.all,start.all,1,gamma);
        DecaDobl_Coefficient_Homotopy.Create(start.all,target.all,1,gamma);
      end if;
    else
      nbequ := target'last;
      DecaDobl_Homotopy.Create(target.all,start.all,1,gamma); -- tpow,gamma);
    end if;
  end DecaDobl_Reader;

  procedure HexaDobl_Reader
              ( nbequ : out integer32;
                sols : out HexaDobl_Complex_Solutions.Solution_List;
               -- tpow : in natural32;
                gamma : in HexaDobl_Complex_Numbers.Complex_Number;
                homcrd,rabin : in boolean := false ) is

    target,start : HexaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    Artificial_Parameter_Homotopy_io.get(target,start,sols);
    if not rabin then
      if homcrd
       then HexaDobl_Projective_Transformation(target,start,sols);
      end if;
      nbequ := target'last;
      if not homcrd then
        HexaDobl_Homotopy.Create(target.all,start.all,1,gamma); -- tpow,gamma);
      else -- set tpow to one for homogeneous coordinates
        HexaDobl_Homotopy.Create(target.all,start.all,1,gamma);
        HexaDobl_Coefficient_Homotopy.Create(start.all,target.all,1,gamma);
      end if;
    else
      nbequ := target'last;
      HexaDobl_Homotopy.Create(target.all,start.all,1,gamma); -- tpow,gamma);
    end if;
  end HexaDobl_Reader;

  procedure Standard_Reader
              ( nbequ : out integer32;
                sols : out Standard_Complex_Solutions.Solution_List;
               -- tpow : in natural32 := 2;
                homcrd,rabin : in boolean := false ) is

    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Random_Numbers.Random1;

  begin
   -- Standard_Reader(nbequ,sols,tpow,gamma,homcrd,rabin);
    Standard_Reader(nbequ,sols,gamma,homcrd,rabin);
  end Standard_Reader;

  procedure DoblDobl_Reader
              ( nbequ : out integer32;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
               -- tpow : in natural32 := 2;
                homcrd,rabin : in boolean := false ) is

    gamma : constant DoblDobl_Complex_Numbers.Complex_Number
          := DoblDobl_Random_Numbers.Random1;

  begin
   -- DoblDobl_Reader(nbequ,sols,tpow,gamma,homcrd,rabin);
    DoblDobl_Reader(nbequ,sols,gamma,homcrd,rabin);
  end DoblDobl_Reader;

  procedure TripDobl_Reader
              ( nbequ : out integer32;
                sols : out TripDobl_Complex_Solutions.Solution_List;
               -- tpow : in natural32 := 2;
                homcrd,rabin : in boolean := false ) is

    gamma : constant TripDobl_Complex_Numbers.Complex_Number
          := TripDobl_Random_Numbers.Random1;

  begin
   -- TripDobl_Reader(nbequ,sols,tpow,gamma,homcrd,rabin);
    TripDobl_Reader(nbequ,sols,gamma,homcrd,rabin);
  end TripDobl_Reader;

  procedure QuadDobl_Reader
              ( nbequ : out integer32;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
               -- tpow : in natural32 := 2;
                homcrd,rabin : in boolean := false ) is

    gamma : constant QuadDobl_Complex_Numbers.Complex_Number
          := QuadDobl_Random_Numbers.Random1;

  begin
   -- QuadDobl_Reader(nbequ,sols,tpow,gamma,homcrd,rabin);
    QuadDobl_Reader(nbequ,sols,gamma,homcrd,rabin);
  end QuadDobl_Reader;

  procedure PentDobl_Reader
              ( nbequ : out integer32;
                sols : out PentDobl_Complex_Solutions.Solution_List;
               -- tpow : in natural32 := 2;
                homcrd,rabin : in boolean := false ) is

    gamma : constant PentDobl_Complex_Numbers.Complex_Number
          := PentDobl_Random_Numbers.Random1;

  begin
   -- PentDobl_Reader(nbequ,sols,tpow,gamma,homcrd,rabin);
    PentDobl_Reader(nbequ,sols,gamma,homcrd,rabin);
  end PentDobl_Reader;

  procedure OctoDobl_Reader
              ( nbequ : out integer32;
                sols : out OctoDobl_Complex_Solutions.Solution_List;
               -- tpow : in natural32 := 2;
                homcrd,rabin : in boolean := false ) is

    gamma : constant OctoDobl_Complex_Numbers.Complex_Number
          := OctoDobl_Random_Numbers.Random1;

  begin
   -- OctoDobl_Reader(nbequ,sols,tpow,gamma,homcrd,rabin);
    OctoDobl_Reader(nbequ,sols,gamma,homcrd,rabin);
  end OctoDobl_Reader;

  procedure DecaDobl_Reader
              ( nbequ : out integer32;
                sols : out DecaDobl_Complex_Solutions.Solution_List;
               -- tpow : in natural32 := 2;
                homcrd,rabin : in boolean := false ) is

    gamma : constant DecaDobl_Complex_Numbers.Complex_Number
          := DecaDobl_Random_Numbers.Random1;

  begin
   -- DecaDobl_Reader(nbequ,sols,tpow,gamma,homcrd,rabin);
    DecaDobl_Reader(nbequ,sols,gamma,homcrd,rabin);
  end DecaDobl_Reader;

  procedure HexaDobl_Reader
              ( nbequ : out integer32;
                sols : out HexaDobl_Complex_Solutions.Solution_List;
               -- tpow : in natural32 := 2;
                homcrd,rabin : in boolean := false ) is

    gamma : constant HexaDobl_Complex_Numbers.Complex_Number
          := HexaDobl_Random_Numbers.Random1;

  begin
   -- HexaDobl_Reader(nbequ,sols,tpow,gamma,homcrd,rabin);
    HexaDobl_Reader(nbequ,sols,gamma,homcrd,rabin);
  end HexaDobl_Reader;

  procedure Standard_Parameter_Reader
              ( nbequ,nbvar,idxpar : out integer32;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Parameter_Systems;

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    nbpar : integer32;

  begin
    Read_Parameter_Homotopy(lp,sols,nbequ,nbvar,nbpar);
    declare
      par : Standard_Integer_Vectors.Vector(1..nbpar);
    begin
      par := Define_Parameters(nbequ,nbvar,nbpar);
      idxpar := par(1);
    end;
    Standard_Homotopy.Create(lp.all,idxpar);
  end Standard_Parameter_Reader;

  procedure DoblDobl_Parameter_Reader
              ( nbequ,nbvar,idxpar : out integer32;
                sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Parameter_Systems;

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    nbpar : integer32;

  begin
    Read_Parameter_Homotopy(lp,sols,nbequ,nbvar,nbpar);
    declare
      par : Standard_Integer_Vectors.Vector(1..nbpar);
    begin
      par := Define_Parameters(nbequ,nbvar,nbpar);
      idxpar := par(1);
    end;
    DoblDobl_Homotopy.Create(lp.all,idxpar);
  end DoblDobl_Parameter_Reader;

  procedure TripDobl_Parameter_Reader
              ( nbequ,nbvar,idxpar : out integer32;
                sols : out TripDobl_Complex_Solutions.Solution_List ) is

    use TripDobl_Parameter_Systems;

    lp : TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    nbpar : integer32;

  begin
    Read_Parameter_Homotopy(lp,sols,nbequ,nbvar,nbpar);
    declare
      par : Standard_Integer_Vectors.Vector(1..nbpar);
    begin
      par := Define_Parameters(nbequ,nbvar,nbpar);
      idxpar := par(1);
    end;
    TripDobl_Homotopy.Create(lp.all,idxpar);
  end TripDobl_Parameter_Reader;

  procedure QuadDobl_Parameter_Reader
              ( nbequ,nbvar,idxpar : out integer32;
                sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Parameter_Systems;

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    nbpar : integer32;

  begin
    Read_Parameter_Homotopy(lp,sols,nbequ,nbvar,nbpar);
    declare
      par : Standard_Integer_Vectors.Vector(1..nbpar);
    begin
      par := Define_Parameters(nbequ,nbvar,nbpar);
      idxpar := par(1);
    end;
    QuadDobl_Homotopy.Create(lp.all,idxpar);
  end QuadDobl_Parameter_Reader;

  procedure PentDobl_Parameter_Reader
              ( nbequ,nbvar,idxpar : out integer32;
                sols : out PentDobl_Complex_Solutions.Solution_List ) is

    use PentDobl_Parameter_Systems;

    lp : PentDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    nbpar : integer32;

  begin
    Read_Parameter_Homotopy(lp,sols,nbequ,nbvar,nbpar);
    declare
      par : Standard_Integer_Vectors.Vector(1..nbpar);
    begin
      par := Define_Parameters(nbequ,nbvar,nbpar);
      idxpar := par(1);
    end;
    PentDobl_Homotopy.Create(lp.all,idxpar);
  end PentDobl_Parameter_Reader;

  procedure OctoDobl_Parameter_Reader
              ( nbequ,nbvar,idxpar : out integer32;
                sols : out OctoDobl_Complex_Solutions.Solution_List ) is

    use OctoDobl_Parameter_Systems;

    lp : OctoDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    nbpar : integer32;

  begin
    Read_Parameter_Homotopy(lp,sols,nbequ,nbvar,nbpar);
    declare
      par : Standard_Integer_Vectors.Vector(1..nbpar);
    begin
      par := Define_Parameters(nbequ,nbvar,nbpar);
      idxpar := par(1);
    end;
    OctoDobl_Homotopy.Create(lp.all,idxpar);
  end OctoDobl_Parameter_Reader;

  procedure DecaDobl_Parameter_Reader
              ( nbequ,nbvar,idxpar : out integer32;
                sols : out DecaDobl_Complex_Solutions.Solution_List ) is

    use DecaDobl_Parameter_Systems;

    lp : DecaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    nbpar : integer32;

  begin
    Read_Parameter_Homotopy(lp,sols,nbequ,nbvar,nbpar);
    declare
      par : Standard_Integer_Vectors.Vector(1..nbpar);
    begin
      par := Define_Parameters(nbequ,nbvar,nbpar);
      idxpar := par(1);
    end;
    DecaDobl_Homotopy.Create(lp.all,idxpar);
  end DecaDobl_Parameter_Reader;

  procedure HexaDobl_Parameter_Reader
              ( nbequ,nbvar,idxpar : out integer32;
                sols : out HexaDobl_Complex_Solutions.Solution_List ) is

    use HexaDobl_Parameter_Systems;

    lp : HexaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    nbpar : integer32;

  begin
    Read_Parameter_Homotopy(lp,sols,nbequ,nbvar,nbpar);
    declare
      par : Standard_Integer_Vectors.Vector(1..nbpar);
    begin
      par := Define_Parameters(nbequ,nbvar,nbpar);
      idxpar := par(1);
    end;
    HexaDobl_Homotopy.Create(lp.all,idxpar);
  end HexaDobl_Parameter_Reader;

  procedure Standard_Series_Newton
              ( sol : in Standard_Complex_Vectors.Vector;
                idx,nbequ : in integer32; nbterms,nbiters : in natural32;
                srv,eva : out Standard_Complex_Series_Vectors.Vector ) is

    hom : constant Standard_Complex_Poly_Systems.Poly_Sys(1..nbequ)
        := Standard_Homotopy.Homotopy_System;
    sys : Standard_CSeries_Poly_Systems.Poly_Sys(1..nbequ)
        := Series_and_Homotopies.Create(hom,idx);
    maxdeg : constant integer32 := integer32(nbterms);
    nit : constant integer32 := integer32(nbiters);

  begin
    Series_and_Predictors.Newton_Prediction(maxdeg,nit,sys,sol,srv,eva);
   -- Standard_Complex_Poly_Systems.Clear(hom); -- sharing, do not clear!
    Standard_CSeries_Poly_Systems.Clear(sys);
  end Standard_Series_Newton;

  procedure Standard_Series_Newton
              ( sol : in Standard_Complex_Solutions.Solution;
                idx,nbequ : in integer32; nbterms,nbiters : in natural32;
                srv,eva : out Standard_Complex_Series_Vectors.Vector ) is
  begin
    Standard_Series_Newton(sol.v,idx,nbequ,nbterms,nbiters,srv,eva);
  end Standard_Series_Newton;

  procedure DoblDobl_Series_Newton
              ( sol : in DoblDobl_Complex_Vectors.Vector;
                idx,nbequ : in integer32; nbterms,nbiters : in natural32;
                srv,eva : out DoblDobl_Complex_Series_Vectors.Vector ) is

    hom : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nbequ)
        := DoblDobl_Homotopy.Homotopy_System;
    sys : DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..nbequ)
        := Series_and_Homotopies.Create(hom,idx);
    maxdeg : constant integer32 := integer32(nbterms);
    nit : constant integer32 := integer32(nbiters);

  begin
    Series_and_Predictors.Newton_Prediction(maxdeg,nit,sys,sol,srv,eva);
   -- DoblDobl_Complex_Poly_Systems.Clear(hom); -- sharing, do not clear!
    DoblDobl_CSeries_Poly_Systems.Clear(sys);
  end DoblDobl_Series_Newton;

  procedure DoblDobl_Series_Newton
              ( sol : in DoblDobl_Complex_Solutions.Solution;
                idx,nbequ : in integer32; nbterms,nbiters : in natural32;
                srv,eva : out DoblDobl_Complex_Series_Vectors.Vector ) is
  begin
    DoblDobl_Series_Newton(sol.v,idx,nbequ,nbterms,nbiters,srv,eva);
  end DoblDobl_Series_Newton;

  procedure QuadDobl_Series_Newton
              ( sol : in QuadDobl_Complex_Vectors.Vector;
                idx,nbequ : in integer32; nbterms,nbiters : in natural32;
                srv,eva : out QuadDobl_Complex_Series_Vectors.Vector ) is

    hom : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nbequ)
        := QuadDobl_Homotopy.Homotopy_System;
    sys : QuadDobl_CSeries_Poly_Systems.Poly_Sys(1..nbequ)
        := Series_and_Homotopies.Create(hom,idx);
    maxdeg : constant integer32 := integer32(nbterms);
    nit : constant integer32 := integer32(nbiters);

  begin
    Series_and_Predictors.Newton_Prediction(maxdeg,nit,sys,sol,srv,eva);
   -- QuadDobl_Complex_Poly_Systems.Clear(hom); -- sharing, do not clear!
    QuadDobl_CSeries_Poly_Systems.Clear(sys);
  end QuadDobl_Series_Newton;

  procedure QuadDobl_Series_Newton
              ( sol : in QuadDobl_Complex_Solutions.Solution;
                idx,nbequ : in integer32; nbterms,nbiters : in natural32;
                srv,eva : out QuadDobl_Complex_Series_Vectors.Vector ) is
  begin
    QuadDobl_Series_Newton(sol.v,idx,nbequ,nbterms,nbiters,srv,eva);
  end QuadDobl_Series_Newton;

end Homotopy_Series_Readers;
