\page development Development Guidelines
\brief Guidelines for development of the FPR code

# Development Guidelines

Updated 02/2019

---

- Write using C++11, conforming to completeness of GCC 4.9.0

---

## Code Style

### General Style

  - Follow the [Webkit style guide](https://webkit.org/code-style-guidelines/) for actual code style
  - Code column limit of 120.
  - If using `clang-format` (STRONGLY recommended), use this command line flag: 

        -style={BasedOnStyle : Webkit, IndentWidth : 4, ColumnLimit : 120, CommentPragmas: '/\*(.+\n.+)+\*/'}

### Brackets

  - Bracket usage follows Webkit style:
    - New line before brackets in functions only:

          // correct
          void some_func(Args...)
          {
              ...
          }

          // incorrect
          void some_func(Args...) {
              ...
          }

    - All other brackets start at the end of the line (e.g. if-else, classes, structs, enums)
    - If an if-else block is only one line, do not use brackets:

          // correct
          if (true)
              some_func();

          // correct
          if (true)
              some_func();
          else {
              some_func();
              some_other_func();
          }

          // incorrect
          if (true) {
              some_func();
          }

    - An exception is if the one line of code is broken into two due to the column limit:

          if (true) {
              some_function_which_spans_more_than_120_columns(arg1, arg2, arg3,
                  arg4);
          }

## Naming Convention

  - Variables follow Camel case (compound phrases wherein each word or abbreviation is capitalized, except the first word):
    - example: `requiredBonds`
  - Classes and Structs follow Pascal Case (subset of Camel case, but the first word is capitalized):
    - example: `MolTemplate`
  - Functions (other than those which return bools) should be named with snake case
    - example: `calculate_association_angles()`
  - Bools of all types should be name so that they are immediately distinguishable in meaning
    - example: `isBound` clearly indicates that the molecule is bound, where as `!isBound`  indicates the opposite

## Uniform initialization 

  - Uniform initialization, i.e. bracket initialization, should be used at all times, with some exceptions. Why?
    - Distinguishes between initialization, assignment, and copy constructor calls

          int x; // uninitialized variable
          int x = 2; // initialized variable with value 2
          int y = x; // calls a copy constructor, not an assignment
          x = y; // an assignment, calls the copy operator=

          int x{ 2 }; // initialized variable with value 2
          int y{ x }; // initialized variable with value 2

      Note that it is very easy to distinguish initialization when using uniform initialization
    - Some exceptions (`std::initializer_list`).
      - Since GCC 4.9.0 support is required, object references and arrays/vectors cannot be initialized with brace initializaton. For example:
          
            // GOOD
            Molecule& mol = moleculeList[0];
            std::array<double, 2> arr = some_function_returning_array(params);

            // BAD
            Molecule& mol{ moleculeList[0] };
            std::array<double, 2> arr{ some_function_returning_array(params) };
