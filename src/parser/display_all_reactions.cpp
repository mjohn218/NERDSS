#include "io/io.hpp"
#include "parser/parser_functions.hpp"

void display_all_reactions(const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns, const std::vector<CreateDestructRxn>& createDestructRxns)
{
    for (auto& forwardRxn : forwardRxns) {
        std::cout << "Forward Rxn: " << &forwardRxn - &forwardRxns[0] << std::endl;
        forwardRxn.display();
        // std::cout << linebreak;
        if (forwardRxn.isReversible) {
            std::cout << "\nBack Rxn: " << &forwardRxn - &forwardRxns[0] << std::endl;
            //            forwardRxn.conjBackRxn->display();
            backRxns[forwardRxn.conjBackRxnIndex].display();
            // std::cout << linebreak;
        } //else
        // std::cout << linebreak;
    }

    if (!createDestructRxns.empty()) {
        std::cout << "Creation and Destruction reactions" << std::endl;
        for (auto& oneRxn : createDestructRxns) {
            std::cout << "Create/Destruct Rxn: " << &oneRxn - &createDestructRxns[0] << std::endl;
            oneRxn.display();
            // if ((&oneRxn - &createDestructRxns[0]) + 1 != createDestructRxns.size())
            //     std::cout << linebreak;
        }
    }
}
