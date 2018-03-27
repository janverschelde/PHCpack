with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with Witness_Sets;
with Sampling_Machine;
with Sampling_Laurent_Machine;
with DoblDobl_Sampling_Machine;
with DoblDobl_Sampling_Laurent_Machine;
with QuadDobl_Sampling_Machine;
with QuadDobl_Sampling_Laurent_Machine;

package body Random_Test_Points is

  function Standard_Random_Point
             ( file : file_type;
               p : Standard_Complex_Poly_Systems.Poly_Sys;
               s : Standard_Complex_Solutions.Solution_List;
               d : natural32 )
             return Standard_Complex_Solutions.Solution_List is

    res : Standard_Complex_Solutions.Solution_List;
    n : constant natural32 := natural32(p'last);
    hyp : Standard_Complex_VecVecs.VecVec(1..integer32(d))
        := Witness_Sets.Random_Hyperplanes(d,n);
    startsol : Standard_Complex_Solutions.Solution(integer32(n));
    newsol : Standard_Complex_Solutions.Solution(integer32(n));

  begin
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Default_Tune_Sampler(0);
    Sampling_Machine.Default_Tune_Refiner;
    startsol := Standard_Complex_Solutions.Head_Of(s).all;
    startsol.t := Standard_Complex_Numbers.Create(0.0);
    Sampling_Machine.Sample(file,true,startsol,hyp,newsol);
    Standard_Complex_Solutions.Add(res,newsol);
    Sampling_Machine.Clear;
    return res;
  end Standard_Random_Point;

  function Standard_Random_Point
             ( file : file_type;
               p : Standard_Complex_Laur_Systems.Laur_Sys;
               s : Standard_Complex_Solutions.Solution_List;
               d : natural32 )
             return Standard_Complex_Solutions.Solution_List is

    res : Standard_Complex_Solutions.Solution_List;
    n : constant natural32 := natural32(p'last);
    hyp : Standard_Complex_VecVecs.VecVec(1..integer32(d))
        := Witness_Sets.Random_Hyperplanes(d,n);
    startsol : Standard_Complex_Solutions.Solution(integer32(n));
    newsol : Standard_Complex_Solutions.Solution(integer32(n));

  begin
    Sampling_Laurent_Machine.Initialize(p);
    Sampling_Laurent_Machine.Default_Tune_Sampler(0);
    Sampling_Laurent_Machine.Default_Tune_Refiner;
    startsol := Standard_Complex_Solutions.Head_Of(s).all;
    startsol.t := Standard_Complex_Numbers.Create(0.0);
    Sampling_Laurent_Machine.Sample(file,true,startsol,hyp,newsol);
    Standard_Complex_Solutions.Add(res,newsol);
    Sampling_Laurent_Machine.Clear;
    return res;
  end Standard_Random_Point;

  function DoblDobl_Random_Point
             ( file : file_type;
               p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               s : DoblDobl_Complex_Solutions.Solution_List;
               d : natural32 )
             return DoblDobl_Complex_Solutions.Solution_List is

    res : DoblDobl_Complex_Solutions.Solution_List;
    n : constant natural32 := natural32(p'last);
    hyp : DoblDobl_Complex_VecVecs.VecVec(1..integer32(d))
        := Witness_Sets.Random_Hyperplanes(d,n);
    startsol : DoblDobl_Complex_Solutions.Solution(integer32(n));
    newsol : DoblDobl_Complex_Solutions.Solution(integer32(n));
    ddzero : double_double := Double_Double_Numbers.Create(0.0);

  begin
    DoblDobl_Sampling_Machine.Initialize(p);
    DoblDobl_Sampling_Machine.Default_Tune_Sampler(0);
    DoblDobl_Sampling_Machine.Default_Tune_Refiner;
    startsol := DoblDobl_Complex_Solutions.Head_Of(s).all;
    startsol.t := DoblDobl_Complex_Numbers.Create(ddzero);
    DoblDobl_Sampling_Machine.Sample(file,true,startsol,hyp,newsol);
    DoblDobl_Complex_Solutions.Add(res,newsol);
    DoblDobl_Sampling_Machine.Clear;
    return res;
  end DoblDobl_Random_Point;

  function DoblDobl_Random_Point
             ( file : file_type;
               p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
               s : DoblDobl_Complex_Solutions.Solution_List;
               d : natural32 )
             return DoblDobl_Complex_Solutions.Solution_List is

    res : DoblDobl_Complex_Solutions.Solution_List;
    n : constant natural32 := natural32(p'last);
    hyp : DoblDobl_Complex_VecVecs.VecVec(1..integer32(d))
        := Witness_Sets.Random_Hyperplanes(d,n);
    startsol : DoblDobl_Complex_Solutions.Solution(integer32(n));
    newsol : DoblDobl_Complex_Solutions.Solution(integer32(n));
    ddzero : double_double := Double_Double_Numbers.Create(0.0);

  begin
    DoblDobl_Sampling_Laurent_Machine.Initialize(p);
    DoblDobl_Sampling_Laurent_Machine.Default_Tune_Sampler(0);
    DoblDobl_Sampling_Laurent_Machine.Default_Tune_Refiner;
    startsol := DoblDobl_Complex_Solutions.Head_Of(s).all;
    startsol.t := DoblDobl_Complex_Numbers.Create(ddzero);
    DoblDobl_Sampling_Laurent_Machine.Sample(file,true,startsol,hyp,newsol);
    DoblDobl_Complex_Solutions.Add(res,newsol);
    DoblDobl_Sampling_Laurent_Machine.Clear;
    return res;
  end DoblDobl_Random_Point;

  function QuadDobl_Random_Point
             ( file : file_type;
               p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               s : QuadDobl_Complex_Solutions.Solution_List;
               d : natural32 )
             return QuadDobl_Complex_Solutions.Solution_List is

    res : QuadDobl_Complex_Solutions.Solution_List;
    n : constant natural32 := natural32(p'last);
    hyp : QuadDobl_Complex_VecVecs.VecVec(1..integer32(d))
        := Witness_Sets.Random_Hyperplanes(d,n);
    startsol : QuadDobl_Complex_Solutions.Solution(integer32(n));
    newsol : QuadDobl_Complex_Solutions.Solution(integer32(n));
    qdzero : quad_double := Quad_Double_Numbers.Create(0.0);

  begin
    QuadDobl_Sampling_Machine.Initialize(p);
    QuadDobl_Sampling_Machine.Default_Tune_Sampler(0);
    QuadDobl_Sampling_Machine.Default_Tune_Refiner;
    startsol := QuadDobl_Complex_Solutions.Head_Of(s).all;
    startsol.t := QuadDobl_Complex_Numbers.Create(qdzero);
    QuadDobl_Sampling_Machine.Sample(file,true,startsol,hyp,newsol);
    QuadDobl_Complex_Solutions.Add(res,newsol);
    QuadDobl_Sampling_Machine.Clear;
    return res;
  end QuadDobl_Random_Point;

  function QuadDobl_Random_Point
             ( file : file_type;
               p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
               s : QuadDobl_Complex_Solutions.Solution_List;
               d : natural32 )
             return QuadDobl_Complex_Solutions.Solution_List is

    res : QuadDobl_Complex_Solutions.Solution_List;
    n : constant natural32 := natural32(p'last);
    hyp : QuadDobl_Complex_VecVecs.VecVec(1..integer32(d))
        := Witness_Sets.Random_Hyperplanes(d,n);
    startsol : QuadDobl_Complex_Solutions.Solution(integer32(n));
    newsol : QuadDobl_Complex_Solutions.Solution(integer32(n));
    qdzero : quad_double := Quad_Double_Numbers.Create(0.0);

  begin
    QuadDobl_Sampling_Laurent_Machine.Initialize(p);
    QuadDobl_Sampling_Laurent_Machine.Default_Tune_Sampler(0);
    QuadDobl_Sampling_Laurent_Machine.Default_Tune_Refiner;
    startsol := QuadDobl_Complex_Solutions.Head_Of(s).all;
    startsol.t := QuadDobl_Complex_Numbers.Create(qdzero);
    QuadDobl_Sampling_Laurent_Machine.Sample(file,true,startsol,hyp,newsol);
    QuadDobl_Complex_Solutions.Add(res,newsol);
    QuadDobl_Sampling_Laurent_Machine.Clear;
    return res;
  end QuadDobl_Random_Point;

end Random_Test_Points;
