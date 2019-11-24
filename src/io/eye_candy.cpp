#include "io/io.hpp"

#include <chrono>
#include <ctime>
#include <iomanip>

std::ostream& bon(std::ostream& os)
{
    /*! \ingroup Text
     * \brief Turns bolding of iostream on (only works for macOS/Unix
     */

    return os << "\e[1m";
}

std::ostream& boff(std::ostream& os)
{
    /*! \ingroup Text
     * \brief Turns bolding of iostream off (only works for macOS/Unix
     */
    return os << "\e[0m";
}

std::ostream& llinebreak(std::ostream& os)
{
    /*! \ingroup Text
     * \brief Just a 20 character long linebreak of dashes
     */

    return os << std::setw(50) << std::setfill('-') << ' ' << std::setfill(' ') << '\n';
};

std::ostream& linebreak(std::ostream& os)
{
    /*! \ingroup Text
     * \brief Just a 10 character long linebreak of dashes
     */

    return os << std::setw(20) << std::setfill('-') << ' ' << std::setfill(' ') << '\n';
}
