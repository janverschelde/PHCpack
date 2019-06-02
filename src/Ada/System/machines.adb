with System_Call;
with text_io;                       use text_io;

package body Machines is

  function getpid return integer;
 -- pragma interface(C, getpid); -- obsolescent featur
  pragma import(C, getpid);

  function Process_ID return integer is
  begin
    return getpid;
  end Process_ID;

  function Integer_to_String ( i : integer ) return string is

    answer : string(1..20);
    nb : natural := 0;
    tmp : integer := i;

    function Number_to_Character ( n : integer ) return character is
    begin
      case n is
        when 0 => return '0';
        when 1 => return '1';
        when 2 => return '2';
        when 3 => return '3';
        when 4 => return '4';
        when 5 => return '5';
        when 6 => return '6';
        when 7 => return '7';
        when 8 => return '8';
        when 9 => return '9';
        when others => return ' ';
      end case;
    end Number_to_Character;

  begin
    if tmp = 0
     then nb := nb + 1;
	  answer(nb) := Number_to_Character(0);
    end if;
    while tmp /= 0 loop
      nb := nb + 1;
      answer(nb) := Number_to_Character(tmp mod 10);
      tmp := tmp / 10;
    end loop;
    declare
      res : string(1..nb);
    begin
      for j in res'range loop
	res(j) := answer(nb-j+1);
      end loop;
      return res;
    end;
  end Integer_to_String;

  function Process_ID return string is
  begin
    return Integer_to_String(getpid);
  end Process_ID;

  function User_Name ( pid : string ) return string is

    temp : file_type;
    name : string(1..80);
    last : natural;

  begin
    System_Call.Call("whoami > /tmp/user_name" & pid);
    Open(temp,in_file,"/tmp/user_name" & pid);
    get_line(temp,name,last);
    Close(temp);
    System_Call.Call("rm /tmp/user_name" & pid);
    return name(1..last);
  exception
    when others => return "???";
  end User_Name;

  function Architecture ( pid : string ) return string is

    temp : file_type;
    answer : string(1..80);
    last : natural;

  begin
    System_Call.Call("uname -a > /tmp/arch_type" & pid);
    Open(temp,in_file,"/tmp/arch_type" & pid);
    get_line(temp,answer,last);
    Close(temp);
    System_Call.Call("rm /tmp/arch_type" & pid);
    return answer(1..last);
  exception
    when others => return "???";
  end Architecture;

  function Architecture ( pid : string; machine : string ) return string is

    temp : file_type;
    answer : string(1..80);
    last : natural;

  begin
    System_Call.Call("rsh " &  machine & " uname -a > /tmp/arch_type" & pid);
    Open(temp,in_file,"/tmp/arch_type" & pid);
    get_line(temp,answer,last);
    Close(temp);
    System_Call.Call("rm /tmp/arch_type" & pid);
    return answer(1..last);
  exception
    when others => return "???";
  end Architecture;

  function Host_Name ( pid : string ) return string is

    temp : file_type;
    answer : string(1..80);
    last : natural;

  begin
    System_Call.Call("hostname > /tmp/host_name" & pid);
    Open(temp,in_file,"/tmp/host_name" & pid);
    get_line(temp,answer,last);
    Close(temp);
    System_Call.Call("rm /tmp/host_name" & pid);
    return answer(1..last);
  exception
    when others => return "???";
  end Host_Name;

  function Date ( pid : string ) return string is
   
    temp : file_type;
    answer : string(1..80);
    last : natural;

  begin
    System_Call.Call("date > /tmp/date" & pid);
    Open(temp,in_file,"/tmp/date" & pid);
    get_line(temp,answer,last);
    Close(temp);
    System_Call.Call("rm /tmp/date" & pid);
    return answer(1..last);
  exception
    when others => return "???";
  end Date;

end Machines;
