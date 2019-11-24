#include "parser/parser_functions.hpp"

void create_conjugate_reaction_itrs(std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns)
{
    forwardRxns.back().conjBackRxnIndex = backRxns.size() - 1;
    backRxns.back().conjForwardRxnIndex = forwardRxns.size() - 1;
}
