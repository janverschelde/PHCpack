with text_io;                         use text_io;
with Interfaces.C.Strings;            use Interfaces.C.Strings;

function hello ( id : integer; action,message : chars_ptr ) return integer is

  package integer_io is new text_io.integer_io(integer);
  use integer_io;

  v_action : constant string := Value(action);
  v_message : constant string := Value(message);

begin
 -- put("The id : "); put(id,1); new_line;
 -- put("The action : "); put(v_action); new_line;
 -- put("The message : "); put(v_message); new_line;
  put(id,1); put(" ");
  put(v_action); put(" ");
  put(v_message); new_line;
  return 0;
end hello;
