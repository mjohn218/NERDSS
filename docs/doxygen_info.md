\page doxygen Using Doxygen
\brief Just a quick tutorial on Doxygen

# Using Doxygen

## Intro
  - For a file to be parsed for documentation, it must first be declared as a file with a header at the top, as such:
    \code
    /*! \file header.hpp
     * \brief A description of the header
     */
    \endcode
    - Without this, nothing from the file will appear in the documentation
  - Has good support for markdown, including code blocks, with three versions
    \code
    /*!
     * # Header
     *   - Spaces
     *      
     *      This will be handled as a code block
     *
     *   - Tildes
     *     ~~~
     *     So will this
     *     ~~~
     *   - Ticks
     *     ```
     *     And this
     *     ```
     */
    \endcode
      - In line code with ticks, will work
    
## Comment Style
  - Two main types of comments: object and after variable.

### Object Comments
  - Comments an object, such as a function or a class.
  - Style:
    - Traditional:
      \code
      /*!
       * \brief Here is a brief description.
       *
       * A more extensive description is here
       */
      \endcode
    - JavaDoc (must be declared in the doxygen file with `JAVADOC_AUTOBRIEF = YES`
      \code
      /**
       * The first complete sentence, ending in a period, will automatically
       * be used as the brief description.
       * With anything after the period being used for the extensive description.
       */
      \endcode
    - C++ Style Blocks:
      \code
      /// Follows the same idea as JavaDoc, with the first full sentence being used
      /// as the brief description.
      /// The extensive description will then follow.
      \endcode
      - Also,
        \code
        ///! Can also use exclamation marks, if you want to avoid using the required two lines.
        \endcode
        
  - `brief` starts a brief description of the object, which will appear at the top of the doxygen page.
  - A more extensive description, if needed, starts two lines below the brief description and will appear lower on the page.
  - If commenting a function, you can use the `param` keys to comment the purpose of parameters, using `[in]` to denote parameters and `[out]` to denote returned variables:
    \code
    int sum(int x, int y)
    {
        /*!
         * \brief Adds two integers
         *
         * \param x[in] first integer to add
         * \param y[in] second integer to add
         * \param [out] sum of x and y
         */
         ...
    }
    \endcode
  - Function documentation may appear in one of two places:
    1. If a class member function, it will appear under the class of which it is a member.
    2. If a global function, it will appear under Files->Globals->Functions
  - If commenting an object, the `struct` or `class`, whichever appropriate, should follow the exclamation mark, as such:
    \code
    class A {
        /*! \class A
         * \brief A random class
         */
         ...
    };
    \endcode

### After Variable Comments
  - After variable comments should really only be used in classes/structs/enums (basically something with a namespace)
  - Several different styles:
    \code
    int x; //!< Brief description

    int y; ///< Brief description

    int z; /*! Detailed description */

    int a; /**< Detailed descripton */

    int b; //!< Detailed description
           //!< when spanning two lines
    in c; /// Another detailed description
          /// spanning two lines
    \endcode
  - After variable comments can also be used within a function signature (though it looks really ugly and obfuscates the code, in my opinion):
    \code
    int sum(int x /**< [in] first integer */, int y /**< [in] second integer */)
    {
        ...
    }
    \endcode

### Full example

  \code
  /*! \class Coordinate
   * \brief Holds xyz coordinates.
   */
  class Coordinate {
    private:
      double x; //!< x-axis coordinate
      double y; //!< y-axis coordinate
      double z; //!< z-axis coordinate

    public:
      int get_z() {return z;} /*! Returns the z coordinate */
  };
  \endcode

### Some available keywords:
  - Note that all the following commands should have a leading backslash `\`
  - `class` - document a class
  - `struct` - document a struct
  - `enum` - document a enumeration
  - `fn` - document a function (doesn't need to be used often, if the comment is within the function itself)
  - `file` - document a file (see above)
  - `namespace` - document a namespace
