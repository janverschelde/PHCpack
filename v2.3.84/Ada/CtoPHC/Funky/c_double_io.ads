with text_io;                      use text_io;
with Interfaces.C;

package C_Double_io is new text_io.float_io(Interfaces.C.double);
