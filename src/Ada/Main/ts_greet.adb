with text_io;                            use text_io;
with Greeting_Banners;

procedure ts_greet is

-- DESCRIPTION :
--   Prints the version string of PHCpack.

begin
  put_line(Greeting_Banners.Version);
  Greeting_Banners.How_to_Cite;
end ts_greet;
