# Links

http://www.dei.unipd.it/~fisch/ricop/OR2/

# TODO
 - CLI parser for args
   - setup argp
   - options
     - filename
     - verbosity
     - seed
     - n_threads
 - setup git and share with colleagues 
 - install TSPlib
 - reader for TSPlib
   - DIMENSIONS: #nodes
   - then the nodes themselves
 - basic plotter with GNUplot
 - ensure Debugging works correctly
 
 - install CPLEX
 - read on TSP 
 - read on simplex
 - read on C profilers

 - investigare make VS Cmake and other methods to compile more complex software
 - 

# CLI
either get_opt from the C standard or argp from GNU
i chose argp
## ArgP
you use it's structure to handle inputs and the use a switch system to make it
it handles printing, parsing most stuff and things like help and version
https://girishjoshi.io/post/glibc-argument-parsing-argp/
https://nongnu.askapache.com/argpbook/step-by-step-into-argp.pdf
