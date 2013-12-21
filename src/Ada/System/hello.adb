with text_io; use text_io;

procedure hello is

  type natural32 is range 0..2**31-1;
  type integer64 is range -2**63..2**63-1;
  n : natural32;
  a,b : integer64;
  i : integer;
  package integer64_io is new text_io.integer_io(integer64);
  use integer64_io;
  package natural32_io is new text_io.integer_io(natural32);
  use natural32_io;
  package integer_io is new text_io.integer_io(integer);
  use integer_io;

begin
  put_line("Hello world!");
  put("Give ordinary natural : "); get(n);
  put("your integer : "); put(n); new_line;
  put("Give ordinary integer : "); get(i);
  put("your integer : "); put(i); new_line;
  put("Give integer number a : "); get(a);
  put("your integer a : "); put(a); new_line;
  put("Give integer number b : "); get(b);
  put("your integer b : "); put(b); new_line;
  put(" a + b : "); put(a+b); new_line;
end hello;
