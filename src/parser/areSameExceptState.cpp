#include "parser/parser_functions.hpp"

bool areSameExceptState(const ParsedMol::IfaceInfo& iface1, const ParsedMol::IfaceInfo& iface2)
{
    return (iface1.ifaceName == iface2.ifaceName) && (iface1.isBound == iface2.isBound)
        && (iface1.state != iface2.state) && (iface1.speciesIndex == iface2.speciesIndex);
}
