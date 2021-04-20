#include "parser/parser_functions.hpp"

void read_internal_coordinates(std::ifstream& molFile, MolTemplate& molTemplate)
{
    std::string iface;
    double x { 0 };
    double y { 0 };
    double z { 0 };

    auto initialPos = molFile.tellg();
    while (molFile >> iface >> x >> y >> z) {
        std::string tmpIface { iface };
        std::transform(tmpIface.begin(), tmpIface.end(), tmpIface.begin(), ::tolower);
        if (tmpIface == "com")
            molTemplate.comCoord = Coord { x, y, z };
        else {
            molTemplate.interfaceList.emplace_back(iface, Coord { x, y, z });
        }
        std::cout << iface << ": [" << x << "nm, " << y << "nm, " << z << "nm]" << std::endl;
        molFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // ignore the comment

        initialPos = molFile.tellg();
    }
    molFile.clear();
    molFile.seekg(initialPos);
}
