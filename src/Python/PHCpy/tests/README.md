This folder contains test scripts to help the phcpy development.

| file name      |              description                     |
|----------------|----------------------------------------------|
| testimport.py  | does the import phcpy work in a script?      |
| haas.py        | illustrate solver on many cores              |
| apollonius.py  | computes all circles tangent to 3 circles    |
| fourbar.py     | design of a 4-bar mechanism                  |
| fourlines.py   | two lines meeting four lines in 3-space      |
| fourspheres.py | tangent lines to 4 mutually touching spheres |
| showpaths.py   | plots four paths in a simple homotopy        |
| showpoles.py   | shows the poles closest to four paths        |
| minors2x2.py   | decomposition of adjacent 2-by-2 minors      |
| sevenbar.py    | irreducible decomposition of 7-bar design    |
| touchcircle.py | diagonal homotopy to touch a circle          |

The apollonius.py, fourbar.py, and fourlines.py correspond to
the use cases in the tutorial of the phcpy documentation.
The plots by showpaths.py are made with aposteriori step size control,
and plots by showpoles.py are done with apriori step size control.
The adjacent 2-by-2 minors problem runs the monodromy algorithm
to decompose a positive dimensional set into irreducible factors.
The design of a 7-bar mechanism is formulated as a Laurent system.
