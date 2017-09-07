with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Multprec_Natural_Numbers;           use Multprec_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Arrays_of_Floating_Vector_Lists;    use Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Black_Box_Root_Counters;            use Black_Box_Root_Counters;
with Root_Counters_Output;               use Root_Counters_Output;

procedure ts_blkroco is

-- DESCRIPTION :
--   Development of writing the root counter output to string.

  procedure Count ( p : in out Poly_Sys ) is

  -- DESCRIPTION :
  --   Applies the black box root counters to the system p.

    neq : constant natural32 := natural32(p'length);
    tode,mhbz,setb : natural64;
    mivo,stmv : natural32;
    mptode : Natural_Number;
    deg : constant boolean := false;
    zz : Partition(1..neq);
    nz : natural32;
    stlb : double_float;
    lifsup : Link_to_Array_of_Lists; 
    mix,perm,iprm : Link_to_Vector;
    orgmcc,stbmcc : Mixed_Subdivision;
    rocotime : duration;

  begin
    Count_Roots(p,deg,tode,mptode,mhbz,setb,mivo,stmv,
                zz,nz,stlb,lifsup,mix,perm,iprm,orgmcc,stbmcc,rocotime);
    Write_Root_Counts
      (standard_output,deg,tode,mptode,nz,mhbz,setb,mivo,stmv,zz);
  end Count;

  procedure Main is

    lp : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(lp);
    new_line;
    put_line("Your polynomial system :"); put(lp.all);
    Count(lp.all);
  end Main;

begin
  Main;
end ts_blkroco;
