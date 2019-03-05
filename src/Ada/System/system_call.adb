with System;                             use System;

package body System_Call is

  procedure Call ( Command: in string ) is

    function system(command: address) return natural;
   -- pragma Interface(C, system); -- obsolescent feature
    pragma Import(C, system);
    cmd: constant string := command & ASCII.NUL;
    ret: Natural;

  begin
    ret := system(cmd'address);
    if ret /= 0 then
      raise System_Error;
    end if;
  end Call;

end System_Call;
