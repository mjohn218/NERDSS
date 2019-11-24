#include "parser/parser_functions.hpp"

void remove_comment(std::string& line)
{
    // getline(file, line, '#') doesn't work, because it just keeps reading until it finds one
    auto commentPos = line.find('#');
    if (commentPos != std::string::npos)
        line.erase(commentPos, std::string::npos);
}
