# Links

http://www.dei.unipd.it/~fisch/ricop/OR2/

# TODO
## Basic
 - basic plotter with GNUplot
 - utilities
   - make cmd parser safe if not enough commands are submitted
## Important
 - install CPLEX
 - read on TSP 
 - read on simplex
## Extra
 - read on C profilers
 - investigate make vs CMake and other methods to compile more complex software
 - 

# Misc
odd behavior when printing floats: if called with "%d" the result is wrong, but no error is raised. Maybe a compiler flag?

# CLI
either get_opt from the C standard or argp from GNU
## ArgP
You use its structure to handle inputs and the use a switch system 
It handles printing, parsing most stuff and things like help and version
https://girishjoshi.io/post/glibc-argument-parsing-argp/
https://nongnu.askapache.com/argpbook/step-by-step-into-argp.pdf

# TSPLib
http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/
http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/TSPFAQ.html
http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp95.pdf


# Utility

## Log
[Best pracitces](https://dev.to/raysaltrelli/logging-best-practices-obo)
[C code snippet](https://tuttlem.github.io/2012/12/08/simple-logging-in-c.html)