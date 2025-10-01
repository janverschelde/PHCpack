with System.Multiprocessors;

function Number_of_Cores return integer32 is

  res : constant integer32
      := integer32(System.Multiprocessors.Number_of_CPUs);

begin
  return res;
end Number_of_Cores;
