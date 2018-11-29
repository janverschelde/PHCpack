with Standard_Random_Numbers;

package body Homotopy_Continuation_Parameters is

  function Default_Values return Parameters is

    res : Parameters;

  begin
    res.alpha := 1.0E-4;
    res.sbeta := 0.5;
    res.pbeta := 0.5;
    res.gamma := Standard_Random_Numbers.Random1;
    res.tolres := 1.0E-12;
    res.epsilon := 1.0E-12;
    res.numdeg := 4;
    res.dendeg := 4;
    res.maxsize := 0.1;
    res.minsize := 1.0E-6;
    res.corsteps := 4;
    res.maxsteps := 1000;
    return res;
  end Default_Values;

end Homotopy_Continuation_Parameters;
