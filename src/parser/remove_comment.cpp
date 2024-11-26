#include "parser/parser_functions.hpp"

void remove_comment(std::string& line)
{
    // getline(file, line, '#') doesn't work, because it just keeps reading until it finds one
    auto commentPos = line.find('#');
    if (commentPos != std::string::npos)
        line.erase(commentPos, std::string::npos);
}

// std::string create_tmp_line(const std::string &line) {
//   std::string tmpLine{line};
//   std::transform(tmpLine.begin(), tmpLine.end(), tmpLine.begin(), ::tolower);
//   tmpLine.erase(std::remove_if(tmpLine.begin(), tmpLine.end(),
//                                [](unsigned char x) { return std::isspace(x); }),
//                 tmpLine.end());
//   remove_comment(tmpLine);
//   return tmpLine;
// }