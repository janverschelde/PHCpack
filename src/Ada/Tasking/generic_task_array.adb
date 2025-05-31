with Identification;

procedure Generic_Task_Array ( p : in integer ) is

  task type worker;

  task body worker is

    idnbr : constant integer := Identification.Number;

  begin
    do_job(idnbr);
  end worker;

  workers : array(1..p) of worker;

begin
  null;
end Generic_Task_Array;
