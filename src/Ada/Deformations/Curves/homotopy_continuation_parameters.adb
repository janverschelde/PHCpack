with unchecked_deallocation;
with Standard_Random_Numbers;

package body Homotopy_Continuation_Parameters is

  function Default_Values return Parameters is

    res : Parameters;

  begin
    res.alpha := 1.0E-3;
    res.pbeta := 0.5;
    res.cbeta := 0.005;
    res.gamma := Standard_Random_Numbers.Random1;
    res.tolres := 1.0E-8;
    res.epsilon := 1.0E-12;
    res.numdeg := 5;
    res.dendeg := 1;
    res.maxsize := 0.1;
    res.minsize := 1.0E-6;
    res.corsteps := 4;
    res.maxsteps := 1000;
    return res;
  end Default_Values;

  procedure Clear ( pars : in out Link_to_Parameters ) is

    procedure free is
      new unchecked_deallocation(Parameters,Link_to_Parameters);

  begin
    if pars /= null
     then free(pars);
    end if;
  end Clear;

-- STORAGE of an instance of the parameters :

  parameters_instance : Link_to_Parameters;

  procedure Construct ( pars : in Parameters ) is
  begin
    parameters_instance := new Parameters'(pars);
  end Construct;

  function Retrieve return Link_to_Parameters is
  begin
    return parameters_instance;
  end Retrieve;

  procedure Destruct is
  begin
    if parameters_instance /= null
     then clear(parameters_instance);
    end if;
  end Destruct;

end Homotopy_Continuation_Parameters;
