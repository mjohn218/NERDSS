/*! \file class_parameter.cpp
 * \ingroup Parser
 * ### Created on 5/25/18 by Matthew Varga
 * ### Purpose
 * ***
 * Class member functions for class Parameter
 *
 * ### Notes
 * ***
 *
 * ### TODO List
 * ***
 */
#include "classes/class_Membrane.hpp"
#include "classes/class_Parameters.hpp"
#include "io/io.hpp"
#include "parser/parser_functions.hpp"

double Parameters::dt = 0.0;
std::vector<long long int> Parameters::lastUpdateTransition {};

// this is so we can compare parsed keywords with the enumerations (we need them as strings)
std::map<const std::string, ParamKeyword> parmKeywords = {
    { "nummoltypes", ParamKeyword::numMolTypes }, { "numtotalspecies", ParamKeyword::numTotalSpecies },
    { "nitr", ParamKeyword::nItr },
    { "fromrestart", ParamKeyword::fromRestart }, { "timewrite", ParamKeyword::timeWrite },
    { "trajwrite", ParamKeyword::trajWrite }, { "timestep", ParamKeyword::timeStep },
    { "numtotalcomplex", ParamKeyword::numTotalComplex },
    { "mass", ParamKeyword::mass }, { "restartwrite", ParamKeyword::restartWrite },
    { "pdbwrite", ParamKeyword::pdbWrite },
    { "overlapseplimit", ParamKeyword::overlapSepLimit }, { "name", ParamKeyword::name },
    { "checkpoint", ParamKeyword::checkPoint }, { "scalemaxdisplace", ParamKeyword::scaleMaxDisplace },
    { "transitionwrite", ParamKeyword::transitionWrite }, { "clusteroverlapcheck", ParamKeyword::clusterOverlapCheck }
};

void Parameters::set_value(std::string value, ParamKeyword keywords)
{
    /*! \ingroup Parser
     * \brief Sets the parameters based on the enumeration key.
     * @param value value of the parameter as a string
     * @param keywords the keyword parsed from the input file to match to the enumeration Keywords
     */
    double nit, checkit, transitionit;
    try {
        auto key = static_cast<std::underlying_type<ParamKeyword>::type>(keywords);
        switch (key) {
        case 0:
            this->numMolTypes = std::stoi(value);
            std::cout << "Read in numMolTypes: " << this->numMolTypes << std::endl;
            break;
        case 1:
            this->numTotalSpecies = std::stoi(value);
            std::cout << "Read in numTotalSpecies: " << this->numTotalSpecies << std::endl;
            break;
        case 2:
            nit = std::stod(value);
            this->nItr = (long long)(nit); //std::stoi(value);
            std::cout << "Read in nItr: " << this->nItr << " timeSteps" << std::endl;
            break;
        case 3:
            this->fromRestart = read_boolean(value);
            std::cout << "Read in fromRestart: " << std::boolalpha << this->fromRestart << std::endl;
            break;
        case 4:
            this->timeWrite = std::stoi(value);
            std::cout << "Read in timeWrite: " << this->timeWrite << " timeSteps" << std::endl;
            break;
        case 5:
            this->trajWrite = std::stoi(value);
            std::cout << "Read in trajWrite: " << this->trajWrite << " timeSteps" << std::endl;
            break;
        case 6:
            this->timeStep = std::stod(value);
            std::cout << "Read in timeStep: " << this->timeStep << " us" << std::endl;
            break;
        case 7:
            this->numTotalComplex = std::stoi(value);
            std::cout << "Read in numTotalComplex: " << this->numTotalComplex << std::endl;
            break;
        case 8:
            this->mass = std::stod(value);
            std::cout << "Read in mass: " << this->mass << std::endl;
            break;
        case 10:
            this->restartWrite = std::stoi(value);
            std::cout << "Read in restartWrite: " << this->restartWrite << " timeSteps" << std::endl;
            break;
        case 11:
            this->pdbWrite = std::stoi(value);
            std::cout << "Read in pdbWrite: " << this->pdbWrite << " timeSteps" << std::endl;
            break;
        case 12:
            this->overlapSepLimit = std::stod(value);
            std::cout << "Read in overlapSepLimit: " << this->overlapSepLimit << " nm" << std::endl;
            break;
        case 13:
            this->name = value;
            std::cout << "Read in name: " << value << std::endl;
            break;
        case 14:
            checkit = std::stod(value);
            this->checkPoint = (long long)(checkit);
            std::cout << "Read in checkPoint: " << this->checkPoint << " timeSteps" << std::endl;
            break;
        case 15:
	        this->scaleMaxDisplace = std::stod(value);
	        std::cout << "Read in scaleMaxDisplace: " << this->scaleMaxDisplace << std::endl;
            break;
        case 16:
            transitionit = std::stod(value);
            this->transitionWrite = (long long)(transitionit);
            std::cout << "Read in transitionWrite: " << this->transitionWrite << " timeSteps" << std::endl;
            break;
        case 17:
            this->clusterOverlapCheck = read_boolean(value);
            std::cout << "Read in clusterOverlapCheck: " << std::boolalpha << this->clusterOverlapCheck << std::endl;
            break;
        default:
            throw std::invalid_argument("Not a valid keyword.");
        }
    } catch (std::invalid_argument& e) {
        std::cout << e.what() << '\n';
        exit(1);
    }
}

void Parameters::parse_paramFile(std::ifstream& paramFile)
{
    /*! \ingroup Parser
     * \brief Main function to parse the parameters block of an input file
     */

    while (!paramFile.eof()) {
        std::string line;
        // see if we're at the end of the parameter block
        auto startPos = paramFile.tellg();
        getline(paramFile, line);

        // remove spaces
        line.erase(
            std::remove_if(line.begin(), line.end(), [](unsigned char x) { return std::isspace(x); }), line.end());
        std::string tmpLine { line };
        std::transform(tmpLine.begin(), tmpLine.end(), tmpLine.begin(), ::tolower);

        if (tmpLine == "endparameters") {
            paramFile.seekg(startPos);
            return;
        }
        // if the line starts with a comment, skip it
        if (tmpLine[0] == '#')
            continue;
        else
            remove_comment(line);

        bool gotValue { false };
        std::string buffer;
        for (auto lineItr = line.begin(); lineItr != line.end(); ++lineItr) {
            if (std::isalnum(*lineItr))
                buffer += std::tolower(static_cast<char>(*lineItr));
            else if (*lineItr == '=') {
                auto keyFind = parmKeywords.find(buffer);
                line.erase(line.begin(), lineItr + 1); // + 1 removes the '=' sign. could make this erase(remove_if)

                // find the value type from the keyword and then set that parameter
                if (keyFind != parmKeywords.end()) {
                    // std::cout << "Keyword found: " << keyFind->first << '\n';
                    this->set_value(line, keyFind->second);
                    gotValue = true;
                    break;
                } else {
                    std::cout << "Warning, ignoring unknown keyword " << buffer << '\n';
                    break;
                }
            }
        }
    }
}

void Parameters::display()
{
    std::cout << "Number of iterations: " << nItr << " timesteps\n";
    std::cout << "Timestep: " << timeStep << " us\n";
    std::cout << "Timestep log interval (timeWrite): " << timeWrite << " timesteps\n";
    std::cout << "Restart file write interval: " << restartWrite << " timesteps\n";
    std::cout << "Coordinate write interval (trajWrite): " << trajWrite << " timesteps\n";
    std::cout << "PDB Coordinate write interval: " << pdbWrite << " timesteps\n";
    std::cout << "Checkpoint write interval: " << checkPoint << " timesteps\n";
    std::cout << "overlapSepLimit: " << overlapSepLimit << " nm\n";
    std::cout << "Transition matrix write interval: " << transitionWrite << " timesteps\n";
    std::cout << "ClusterOverlapCheck: " << clusterOverlapCheck << "\n";

    std::cout << "Molecule specific parameters:\n";
    std::cout << "Number of unique molecule types: " << numMolTypes << '\n';
    std::cout << "Total number of unique interfaces and states, including product states: " << numTotalSpecies << '\n';
    std::cout << "Total number of complexes in system at start: " << numTotalComplex << '\n';
    std::cout << "Total number of units (molecules + interfaces) in system at start: " << numTotalUnits << '\n';
    std::cout << "Maximum allowed number of unique 2D reactions: " << max2DRxns << '\n';
}
