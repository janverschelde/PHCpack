with text_io;                            use text_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Standard_Dictionary_Solutions_io is

-- DESCRIPTION :
--   This package converts a list of solutions into a dictionary format
--   accessible to Python.  The keys in the dictionary are the names of
--   the variables and attributes, of type string.  The values in the
--   dictionary are one of the following types:  natural (multiplicity),
--   float (for err, rco, and res), or otherwise complex.

  -- FORMAT of one solution x = 2.3 - 1.7*I, I = sqrt(-1):
  -- >>> s = {'x':complex(2.3,-1.7),'err':1.0e-10}
  --  >>> s['x']
  -- (2.3-1.7j)
  -- >>> s['err']
  -- 1e-10
  -- FORMAT of a list of solutions is formatted as a list of dictionaries:
  -- >>> L = [{'x':complex(1.0,-2),'err':1.0e-10},
  --          {'x':complex(2,0),'err':1.0e-15}]
  -- >>> L[0]
  -- {'x': (1-2j), 'err': 1e-10}
  -- >>> L[1]['x']
  -- (2+0j)

  procedure put ( s : in Solution );
  procedure put ( file : in file_type; s : in Solution );
  procedure put ( sols : in Solution_List );
  procedure put ( file : in file_type; sols : in Solution_List );

  -- DESCRIPTION :
  --   Write the dictionary format (see above) of the solution (list)
  --   to standard output or to file.

--  function Write ( s : Solution ) return string;
--  function Write ( sols : Solution_List ) return string;

  -- DESCRIPTION :
  --   Returns a string writing the solution(s) in dictionary format.

  procedure get ( sols : out Solution_List );
  procedure get ( file : in file_type; sols : out Solution_List );

  -- DESCRIPTION :
  --   Reads solutions in dictionary format from standard input
  --   or from file, in the format produced by put.

end Standard_Dictionary_Solutions_io;
