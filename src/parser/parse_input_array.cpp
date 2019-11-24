#include "parser/parser_functions.hpp"

std::vector<double> parse_input_array(std::string& line)
{
    // remove the brackets
    std::vector<char> brackets { '[', ']', '(', ')' };
    for (auto& elem : brackets) {
        auto bracketPos = line.find(elem);
        if (bracketPos != std::string::npos)
            line.erase(bracketPos, 1);
    }

    // parse the comma delimited values
    size_t dItr = 0;
    std::vector<std::string> strVals;
    std::vector<double> tmpVals;
    while ((dItr = line.find(',')) != std::string::npos) {
        strVals.emplace_back(line.substr(0, dItr));
        line.erase(0, dItr + 1);
    }
    strVals.emplace_back(line.substr(0, std::string::npos));

    for (auto& value : strVals) {
        try {
            tmpVals.emplace_back(stod(value));
        } catch (const std::exception& e) {
            std::string tmpValue { value };
            std::transform(tmpValue.begin(), tmpValue.end(), tmpValue.begin(), ::tolower);
            bool isPi { tmpValue == "m_pi" || tmpValue == "pi" };
            bool isNull { tmpValue == "nan" };
            if (isPi)
                tmpVals.emplace_back(M_PI);
            else if (isNull)
                tmpVals.emplace_back(std::numeric_limits<double>::quiet_NaN());
            else
                throw "Error, cannot read angles value " + value;
        } catch (const std::string& msg) {
            std::cout << msg << '\n';
            exit(1);
        }
    }

    return tmpVals;
}
