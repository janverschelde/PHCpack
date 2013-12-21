with text_io;                            use text_io;
with Standard_Puiseux_Certificates;      use Standard_Puiseux_Certificates;

package Standard_Puiseux_Certificates_io is

-- DESCRIPTION :
--   Provides output of the first two terms of a Puiseux series expansion.

  procedure put ( g : in Germ );
  procedure put ( file : in file_type; g : in Germ );

  -- DESCRIPTION :
  --   Writes the information in the germ g as the first two 
  --   leading terms of a Puiseux series expansion.

  procedure put ( g : in List_of_Germs );
  procedure put ( file : in file_type; g : in List_of_Germs );

  -- DESCRIPTION :
  --   Writes all germs in the list to standard output or to file.

end Standard_Puiseux_Certificates_io;
