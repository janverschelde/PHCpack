package body Semaphore is

  protected body Binary_Semaphore is

    function Is_Available return boolean is
    begin
      return available;
    end Is_Available;

    entry Seize when available is
    begin
      available := false;
    end Seize;

    procedure Release is
    begin
      available := true;
    end Release;

  end Binary_Semaphore;

  function Available ( s : Lock ) return boolean is
  begin
    return s.b.Is_Available;
  end Available;

  procedure Request ( s : in out Lock ) is
  begin
    s.b.Seize;
  end Request;

  procedure Release ( s : in out Lock ) is
  begin
    s.b.Release;
  end Release;

end Semaphore;
