#pragma once

#define DEBUG false
#define VERBOSE false
#define PROFILE false
#define COUNT_COMMUNICATIONS false
#define RANK0_BUFFER_SIZE 100000000
#define NEIGHBOR_BUFFER_SIZE 50000000

#define TEST_SERIALIZATION_AND_DESERIALIZATION 0
#define PRINT_REACTIONS_PROBABILITIES 0

#ifndef INT_MAX
#define INT_MAX 2147483647
#endif

// Debugging often requires tracking a single molecule throughout the
// simulation. Calls to following macros are inserted around communication
// routines. If the ; sign is replaced with \ sign, this will tell the compiler
// that the definition of the macro spans to the next line.
// Before the example line, one could place different condition for printing
// and different printing. Please keep example line intact,
// so that it can serve for a future reference:
#define DEBUG_MOL(s) \
  { debug_print(mpiContext, mol, s); }
#define DEBUG_FIND_MOL(s)                                           \
  {                                                                 \
    for (auto& mol : moleculeList) debug_print(mpiContext, mol, s); \
  }

// Tracking a complex
#define DEBUG_COMPLEX(s) \
  { debug_print_complex(mpiContext, complex, s); }
#define DEBUG_FIND_COMPLEX(s)                      \
  {                                                \
    for (auto& complex : complexList)              \
      debug_print_complex(mpiContext, complex, s); \
  }
