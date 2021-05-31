with System.Multiprocessors;

function Number_of_Cores return integer32 is

 -- function core_count return integer32;
 -- pragma import(C, core_count, "corecount");
 -- res : constant integer32 := core_count;

  res : constant integer32
      := integer32(System.Multiprocessors.Number_of_CPUs);

begin
  return res;
end Number_of_Cores;
