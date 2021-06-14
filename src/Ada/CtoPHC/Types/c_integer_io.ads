with text_io;                      use text_io;
with Interfaces.C;

package C_Integer_io is new text_io.integer_io(Interfaces.C.int);
