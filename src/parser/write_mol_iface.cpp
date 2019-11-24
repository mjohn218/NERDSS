#include "parser/parser_functions.hpp"

std::string write_mol_iface(std::string mol, std::string iface)
{
    std::string out = mol + "(" + iface + ")";
    return out;
}
