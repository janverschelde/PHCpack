with Homotopy_Series_Readers;           use Homotopy_Series_Readers;

package body Homotopy_Pade_Approximants is

  procedure Standard_Pade_Approximant
              ( sol : in Standard_Complex_Vectors.Vector;
                idx,nbequ,numdeg,dendeg : in integer32;
                nbiters : in natural32;
                srv,eva : out Standard_Dense_Series_Vectors.Vector;
                pv : out Standard_Pade_Approximants.Pade_Vector ) is

    nbterms : constant natural32 := natural32(numdeg+dendeg+1);

  begin
    Standard_Series_Newton(sol,idx,nbequ,nbterms,nbiters,srv,eva);
    pv := Standard_Pade_Approximants.Create(numdeg,dendeg,srv);
  end Standard_Pade_Approximant;

  procedure DoblDobl_Pade_Approximant
              ( sol : in DoblDobl_Complex_Vectors.Vector;
                idx,nbequ,numdeg,dendeg : in integer32;
                nbiters : in natural32;
                srv,eva : out DoblDobl_Dense_Series_Vectors.Vector;
                pv : out DoblDobl_Pade_Approximants.Pade_Vector ) is

    nbterms : constant natural32 := natural32(numdeg+dendeg+1);

  begin
    DoblDobl_Series_Newton(sol,idx,nbequ,nbterms,nbiters,srv,eva);
    pv := DoblDobl_Pade_Approximants.Create(numdeg,dendeg,srv);
  end DoblDobl_Pade_Approximant;

  procedure QuadDobl_Pade_Approximant
              ( sol : in QuadDobl_Complex_Vectors.Vector;
                idx,nbequ,numdeg,dendeg : in integer32;
                nbiters : in natural32;
                srv,eva : out QuadDobl_Dense_Series_Vectors.Vector;
                pv : out QuadDobl_Pade_Approximants.Pade_Vector ) is

    nbterms : constant natural32 := natural32(numdeg+dendeg+1);

  begin
    QuadDobl_Series_Newton(sol,idx,nbequ,nbterms,nbiters,srv,eva);
    pv := QuadDobl_Pade_Approximants.Create(numdeg,dendeg,srv);
  end QuadDobl_Pade_Approximant;

end Homotopy_Pade_Approximants;
