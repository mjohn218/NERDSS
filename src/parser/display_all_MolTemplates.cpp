#include "io/io.hpp"
#include "parser/parser_functions.hpp"

void display_all_MolTemplates(const std::vector<MolTemplate>& molTemplates)
{
    for (auto& oneTemp : molTemplates) {
        oneTemp.display();
        // if (&oneTemp - &molTemplates[0] + 1 != molTemplates.size())
        //     std::cout << linebreak;
    }
}
