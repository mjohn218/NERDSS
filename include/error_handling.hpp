/*! @file error_handling.hpp
 * \brief Custom exception classes
 *
 * ### Created on 4/1/18 by Matthew Varga
 * ### Purpose
 * ***
 * Contains custom error classes to, for example, pass local arguments to catch blocks.
 *
 * ### Notes
 * ***
 */
#pragma once
#include <iomanip>
#include <iostream>
#include <stdexcept>

/*! \defgroup Text
 * \brief Functions to alter text or output.
 */

class Exception : public std::exception {
private:
    const char* passedMsg;

public:
    const char* what() const noexcept { return this->passedMsg; };
    Exception() = default;
    Exception(const char* _passedMsg)
        : passedMsg(_passedMsg)
    {
    }
};

template <class ArgType> class ArgException : public std::invalid_argument {
    /* This is a simple class to throw an argument along with an exception message.
     */
private:
    //    const char* passedMsg;
    ArgType _arg1;
    ArgType _arg2;

public:
    //    const char* what() const noexcept { return this->passedMsg; }
    ArgType get_arg1() const { return _arg1; };
    ArgType get_arg2() const { return _arg2; };

    ArgException() = default;
    ArgException(const std::string& __s, const ArgType _passedArg)
        : invalid_argument(__s)
        , _arg1(_passedArg)
    {
    }
    ArgException(const std::string& __s, const ArgType __arg1, const ArgType __arg2)
        : invalid_argument(__s)
        , _arg1(__arg1)
        , _arg2(__arg2)
    {
    }
};

template <class ArgType> class ArgOOFException : public std::out_of_range {
    /* This is a simple class to throw an argument along with an exception message.
     */
private:
    //    const char* passedMsg;
    ArgType _arg1;
    ArgType _arg2;

public:
    //    const char* what() const noexcept { return this->passedMsg; }
    ArgType get_arg1() const { return _arg1; };
    ArgType get_arg2() const { return _arg2; };

    ArgOOFException() = default;
    ArgOOFException(const std::string& __s, const ArgType _passedArgOOF)
        : out_of_range(__s)
        , _arg1(_passedArgOOF)
    {
    }
    ArgOOFException(const std::string& __s, const ArgType __arg1, const ArgType __arg2)
        : out_of_range(__s)
        , _arg1(__arg1)
        , _arg2(__arg2)
    {
    }
};

template <class X, class Y> void involvement_err(X _func, Y _line)
{
    // error of unrecognized involvement of iface in a reaction
    std::cerr << "INVOLVE_ERR: Unrecognized involvement in function (" << _func << ", line " << _line
              << "). Exiting with code 12\n";
    exit(12);
}

template <class X, class Y> void gen_read_err(X _func, Y _line)
{
    // generic read error
    std::cerr << "READ_ERR: Cannot read from file in function (" << _func << ", line " << _line
              << "). Exiting with code 11.\n";
    exit(2);
}

template <class X, class Y> void crd_read_err(X _func, Y _line)
{
    // coordinate specific read error
    std::cerr << "CRD_READ_ERR: cannot read coordinates in function (" << _func << ", line " << _line
              << "). Exiting with code 11.\n";
    exit(3);
}

template <class X, class Y> void not_in_cont_err(X _func, Y _line)
{
    // error specific to not finding a value in a container within which it absolutely should be
    std::cerr << "INVALID_CONT_VALUE_ERR: cannot find value in container, within which it "
                 "absolutely should be, in function ("
              << _func << ", line " << _line << "). Exiting with code 11.\n";
    exit(4);
}

template <class X, class Y> void unbalanced_rxn_err(X _func, Y _line)
{
    std::cerr << "UNBALANCED_RXN_ERR: The reaction is unbalanced. Error in function (" << _func << ", line" << _line
              << ". Exiting with code 5.\n";
    exit(5);
}

template <class X, class Y> void invalid_rxn(std::string&& message, X _func, Y _line)
{
    std::cerr << "INVALID_RXN_ERR (" << _func << ", " << _line << "): " << message << '\n';
    exit(6);
}
