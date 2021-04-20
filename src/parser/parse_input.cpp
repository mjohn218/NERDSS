#include "classes/class_Membrane.hpp"
#include "io/io.hpp"
#include "parser/parser_functions.hpp"
#include "system_setup/system_setup.hpp"

// this is so we can compare parsed keywords with the enumerations (we need them as strings)
std::map<const std::string, BoundaryKeyword> bcKeywords = {
    { "implicitlipid", BoundaryKeyword::implicitLipid }, { "waterbox", BoundaryKeyword::waterBox },
    { "xbctype", BoundaryKeyword::xBCtype }, { "ybctype", BoundaryKeyword::yBCtype },
    { "zbctype", BoundaryKeyword::zBCtype }, { "issphere", BoundaryKeyword::isSphere }, { "spherer", BoundaryKeyword::sphereR }
};

void Membrane::set_value_BC(std::string value, BoundaryKeyword keywords)
{
    try {
        auto key = static_cast<std::underlying_type<BoundaryKeyword>::type>(keywords);
        switch (key) {
        case 0:
            this->implicitLipid = read_boolean(value);
            break;
        case 1:
            this->waterBox = WaterBox(parse_input_array(value));
            this->isBox = true;
            std::cout << "Read in waterBox: "
                      << "[" << waterBox.x << " nm, " << waterBox.y << " nm, " << waterBox.z << " nm]" << std::endl;
            break;
        case 2:
            this->xBCtype = value;
            std::cout << "Read in xBCtype: "
                      << value << std::endl;
            break;
        case 3:
            this->yBCtype = value;
            std::cout << "Read in yBCtype: "
                      << value << std::endl;
            break;
        case 4:
            this->zBCtype = value;
            std::cout << "Read in zBCtype: "
                      << value << std::endl;
            break;
        case 5:
            this->isSphere = read_boolean(value);
            std::cout << "Read in isSphere: " << std::boolalpha << this->isSphere << std::endl;
            break;
        case 6:
            this->sphereR = std::stod(value);
            std::cout << "Read in sphereR: " << this->sphereR << " nm" << std::endl;
            this->isSphere = true;
            break;
        default:
            throw std::invalid_argument("Not a valid keyword.");
        }
    } catch (std::invalid_argument& e) {
        std::cout << e.what() << '\n';
        exit(1);
    }
}

void Membrane::display()
{
    std::cout << " isSphere? " << std::boolalpha << isSphere << std::endl;
    std::cout << " sphere Radius " << sphereR << std::endl;
    if (isBox == true) {
        std::cout << " BOX geometry, dimensions: " << std::endl;
        std::cout << waterBox.x << ' ' << waterBox.y << ' ' << waterBox.z << std::endl;
    }
    std::cout << " hasImplicitLipid? " << std::boolalpha << implicitLipid << std::endl;
}

void Membrane::create_water_box()
{
    waterBox.x = 2 * sphereR;
    waterBox.y = 2 * sphereR;
    waterBox.z = 2 * sphereR;
}

std::string create_tmp_line(const std::string& line)
{
    std::string tmpLine { line };
    std::transform(tmpLine.begin(), tmpLine.end(), tmpLine.begin(), ::tolower);
    tmpLine.erase(
        std::remove_if(tmpLine.begin(), tmpLine.end(), [](unsigned char x) { return std::isspace(x); }), tmpLine.end());
    remove_comment(tmpLine);
    return tmpLine;
}

void parse_input(std::string& fileName, Parameters& params, std::map<std::string, int>& observableList,
    std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns,
    std::vector<CreateDestructRxn>& createDestructRxns, std::vector<MolTemplate>& molTemplateList, Membrane& membraneObject)
{
    bool hasParsedMol = false;
    std::ifstream inputFile { fileName };
    if (!inputFile) {
        std::cerr << "Reaction file cannot be opened. Exiting..." << std::endl;
        exit(1);
    }

    std::vector<std::string> providedObs {};
    // the start and end blocks are hardcoded since regex has proved itself woefully non-portable,
    // especially with the Intel compiler (icpc)
    while (!inputFile.eof()) {
        std::string line;
        getline(inputFile, line);
        std::string tmpLine { create_tmp_line(line) };

        if (skipLine(tmpLine)) {
            continue;
        } else if (tmpLine == "startparameters") {
            // read in parameters
            std::cout << "Parsing simulation parameters..." << '\n';
            params.parse_paramFile(inputFile);
        } else if (tmpLine == "startboundaries") {
            // read in boundaries
            std::cout << "Parsing simulation boundary conditions..." << '\n';
            while (getline(inputFile, line)) {
                tmpLine = create_tmp_line(line);
                if (skipLine(tmpLine))
                    continue;
                else if (tmpLine == "endboundaries")
                    break;
                else {
                    line.erase(std::remove_if(line.begin(), line.end(), [](unsigned char x) { return std::isspace(x); }), line.end());
                    std::transform(line.begin(), line.end(), line.begin(), ::tolower); //make all lower case

                    // if the line starts with a comment, skip it
                    if (line[0] == '#')
                        continue;
                    else
                        remove_comment(line);

                    bool gotValue { false };
                    std::string buffer;
                    for (auto lineItr = line.begin(); lineItr != line.end(); ++lineItr) {
                        if (std::isalnum(*lineItr))
                            buffer += std::tolower(static_cast<char>(*lineItr));
                        else if (*lineItr == '=') {
                            auto keyFind = bcKeywords.find(buffer);
                            line.erase(line.begin(), lineItr + 1); // + 1 removes the '=' sign. could make this erase(remove_if)
                            // find the value type from the keyword and then set that parameter
                            if (keyFind != bcKeywords.end()) {
                                // std::cout << "Keyword found: " << keyFind->first << '\n';
                                membraneObject.set_value_BC(line, keyFind->second);
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
        } else if (tmpLine == "startmolecules") {
            hasParsedMol = true;
            // get the molecule names and copy numbers from the parameter input file
            std::cout << "Parsing mol file information..." << '\n';
            std::vector<std::string> providedMols {};
            std::vector<std::string> providedMols_temp {};
            std::vector<MolTemplate> molTemplateList_temp {};
            std::vector<ParsedMolNumState> providedNumState {};
            std::vector<int> providedNums {};
            while (getline(inputFile, line)) {
                tmpLine = create_tmp_line(line);
                if (skipLine(tmpLine))
                    continue;
                else if (tmpLine == "endmolecules")
                    break;
                else {
                    line.erase(std::remove_if(line.begin(), line.end(), [](unsigned char x) { return std::isspace(x); }), line.end());
                    remove_comment(line);
                    bool gotValue { false };
                    std::string buffer;
                    for (auto lineItr = line.begin(); lineItr != line.end(); ++lineItr) {
                        if (std::isalnum(*lineItr))
                            buffer += static_cast<char>(*lineItr);
                        else if (*lineItr == ':') {
                            providedMols.emplace_back(buffer);
                            line.erase(line.begin(), lineItr + 1); // + 1 removes the ':' sign. could make this erase(remove_if)
                            // find the value type from the keyword and then set that parameter
                            ParsedMolNumState tmpMolNumState {};
                            tmpMolNumState = parse_number_bngl(line);
                            providedNums.emplace_back(tmpMolNumState.totalCopyNumbers);
                            providedNumState.emplace_back(tmpMolNumState);
                            gotValue = true;
                            break;
                        }
                    }
                    if (gotValue == false) {
                        std::cout << "Please provide the copy number of each molecule in INP file in this format-- molName:100" << std::endl;
                        exit(0);
                    }
                }
            }

            // parse the .mol files for each molecule
            for (int providedMolIndex = 0; providedMolIndex < providedMols.size(); providedMolIndex++) {
                molTemplateList.emplace_back(parse_molFile(providedMols[providedMolIndex])); // do the actual parsing
                molTemplateList.back().molTypeIndex = &providedMols[providedMolIndex] - &providedMols[0];
                molTemplateList.back().copies = providedNums[providedMolIndex];
                molTemplateList.back().startingNumState = providedNumState[providedMolIndex];
                ++MolTemplate::numMolTypes;
                if (molTemplateList.back().isImplicitLipid == false) {
                    params.numTotalComplex += molTemplateList.back().copies;
                    params.numTotalUnits
                        += (molTemplateList.back().copies * (molTemplateList.back().interfaceList.size() + 1));
                    if (molTemplateList.back().isLipid)
                        params.numLipids += molTemplateList.back().copies;
                } else {
                }
            }

            MolTemplate::numEachMolType = std::vector<int>(MolTemplate::numMolTypes);
            params.numMolTypes = MolTemplate::numMolTypes;

            // Give a buffer for the trajectory file (since VMD can't understand trajectories with frames that have
            // different numbers of atoms)
            // TODO: Think about the size of buffer
            //            params.numTotalUnits += params.numTotalUnits / 100;

            // Make sure the writing of restart file and timestep information is in unison
            if (params.restartWrite % params.timeWrite != 0) {
                int tmp = int(params.restartWrite / params.timeWrite);
                if (tmp == 0)
                    tmp = 1;
                params.restartWrite = tmp * params.timeWrite;
            }

            // set the isPoint and isRod for molecule
            determine_shape_molecule(molTemplateList);
        } else if (tmpLine == "startreactions") {
            if (hasParsedMol == false) {
                std::cerr << "Error, molecule section must be before reaction section in input file, exiting...\n";
                exit(1);
            }
            // index of the last template. should. starts out as number of interface states
            int totSpecies { (molTemplateList.back().interfaceList.back().stateList.back().index) }; // last iface index
            int numProvidedRxns { 0 };

            auto linePos = inputFile.tellg();
            std::cout << "Parsing reactions...\n";
            while (getline(inputFile, line)) {
                tmpLine = create_tmp_line(line);
                // if the line is an ignored line, ignore it. if it's an end of block line, break
                if (skipLine(tmpLine)) {
                    linePos = inputFile.tellg();
                    continue;
                } else if (tmpLine == "endreactions")
                    break;
                else {
                    // if it's neither, parse the line
                    inputFile.seekg(linePos);
                    parse_reaction(inputFile, totSpecies, numProvidedRxns, molTemplateList, forwardRxns, backRxns,
                        createDestructRxns, observableList, membraneObject);
                }
                linePos = inputFile.tellg();
            }

            RxnBase::totRxnSpecies = totSpecies + 1; // had to increment this by 1, it was not correct!

        } else if (tmpLine == "startobservables") {
            // TODO: make this use parse_molecule_bngl() and create a vector of ParsedMol
            // Need to change the function to allow for molecules with no explicit interfaces
            std::cout << "Gathering observables...\n";
            while (getline(inputFile, line)) {
                tmpLine = create_tmp_line(line);
                tmpLine.erase(
                    std::remove_if(tmpLine.begin(), tmpLine.end(), [](unsigned char x) { return std::isspace(x); }),
                    tmpLine.end());
                if (skipLine(tmpLine)) {
                    // if the line is a comment line, ignore it.
                    continue;
                } else if (tmpLine == "endobservables") {
                    // if end of parameter block, break out from loop
                    break;
                } else {
                    remove_comment(line);
                    providedObs.emplace_back(line); // save these to parse AFTER reactions
                }
            }
        }
    }

    // TODO: TEMPORARY COUPLED RXN
    for (unsigned rxnItr { 0 }; rxnItr < forwardRxns.size(); ++rxnItr) {
        if (forwardRxns[rxnItr].isCoupled) {
            if (forwardRxns[rxnItr].coupledRxn.absRxnIndex != -1) {
                // for the case destruction coupled to disscociation
                for (const auto& createDestructRxn : createDestructRxns) {
                    if (createDestructRxn.absRxnIndex == forwardRxns[rxnItr].coupledRxn.absRxnIndex) {
                        forwardRxns[rxnItr].coupledRxn.rxnType = ReactionType::destruction;
                        forwardRxns[rxnItr].coupledRxn.relRxnIndex = createDestructRxn.relRxnIndex;
                        forwardRxns[rxnItr].coupledRxn.label = createDestructRxn.rxnLabel;
                        backRxns[forwardRxns[rxnItr].conjBackRxnIndex].coupledRxn = forwardRxns[rxnItr].coupledRxn;
                        break;
                    }
                }

                // for the case uniStateChange coupled to disscociation
                for (const auto& uniStateChangeRxn : forwardRxns) {
                    if (uniStateChangeRxn.rxnType == ReactionType::uniMolStateChange) {
                        if (uniStateChangeRxn.absRxnIndex == forwardRxns[rxnItr].coupledRxn.absRxnIndex) {
                            forwardRxns[rxnItr].coupledRxn.rxnType = ReactionType::uniMolStateChange;
                            forwardRxns[rxnItr].coupledRxn.relRxnIndex = uniStateChangeRxn.relRxnIndex;
                            forwardRxns[rxnItr].coupledRxn.label = uniStateChangeRxn.rxnLabel;
                            backRxns[forwardRxns[rxnItr].conjBackRxnIndex].coupledRxn = forwardRxns[rxnItr].coupledRxn;
                            break;
                        }
                    }
                }
            } else {
                // if the input is a label
                // for the case destruction coupled to disscociation
                for (const auto& createDestructRxn : createDestructRxns) {
                    if (createDestructRxn.rxnLabel == forwardRxns[rxnItr].coupledRxn.label) {
                        forwardRxns[rxnItr].coupledRxn.rxnType = ReactionType::destruction;
                        forwardRxns[rxnItr].coupledRxn.relRxnIndex = createDestructRxn.relRxnIndex;
                        forwardRxns[rxnItr].coupledRxn.absRxnIndex = createDestructRxn.absRxnIndex;
                        backRxns[forwardRxns[rxnItr].conjBackRxnIndex].coupledRxn = forwardRxns[rxnItr].coupledRxn;
                        break;
                    }
                }

                // for the case uniStateChange coupled to disscociation
                for (const auto& uniStateChangeRxn : forwardRxns) {
                    if (uniStateChangeRxn.rxnType == ReactionType::uniMolStateChange) {
                        if (uniStateChangeRxn.rxnLabel == forwardRxns[rxnItr].coupledRxn.label) {
                            forwardRxns[rxnItr].coupledRxn.rxnType = ReactionType::uniMolStateChange;
                            forwardRxns[rxnItr].coupledRxn.relRxnIndex = uniStateChangeRxn.relRxnIndex;
                            forwardRxns[rxnItr].coupledRxn.absRxnIndex = uniStateChangeRxn.absRxnIndex;
                            backRxns[forwardRxns[rxnItr].conjBackRxnIndex].coupledRxn = forwardRxns[rxnItr].coupledRxn;
                            break;
                        }
                    }
                }
            }
        }
    }

    // populate the list of reactions on each state
    populate_reaction_lists(forwardRxns, backRxns, createDestructRxns, molTemplateList);
    for (auto& oneTemp : molTemplateList) {
        // here's where we remove duplicate reaction partners
        std::sort(oneTemp.rxnPartners.begin(), oneTemp.rxnPartners.end());
        oneTemp.rxnPartners.erase(
            std::unique(oneTemp.rxnPartners.begin(), oneTemp.rxnPartners.end()), oneTemp.rxnPartners.end());

        for (auto& oneIface : oneTemp.interfaceList) {
            for (auto& oneState : oneIface.stateList) {
                std::sort(oneState.rxnPartners.begin(), oneState.rxnPartners.end());
                oneState.rxnPartners.erase(
                    std::unique(oneState.rxnPartners.begin(), oneState.rxnPartners.end()), oneState.rxnPartners.end());
            }
        }
    }

    params.numTotalSpecies = RxnBase::totRxnSpecies; // total number of all species possible in the system (includes
        // reactant and product states).

    // TODO TEMPORARY
    params.isNonEQ = createDestructRxns.size() > 0;

    std::cout << '\n'
              << "Input file parsing complete\n";
    // std::cout << "Simulation Parameters\n";
    // params.display();
    // std::cout << "Molecule Information\n";
    // display_all_MolTemplates(molTemplateList);
    // std::cout << "Reactions\n";
    // display_all_reactions(forwardRxns, backRxns, createDestructRxns);
}

void parse_input_for_add(std::string& fileName, Parameters& params, std::map<std::string, int>& observableList,
    std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns,
    std::vector<CreateDestructRxn>& createDestructRxns, std::vector<MolTemplate>& molTemplateList, Membrane& membraneObject, int numDoubleBeforeAdd)
{
    bool hasParsedMol = false;
    int addForwardRxnNum { static_cast<int>(forwardRxns.size()) }; //this is the origin number of forwardRxns for add case
    int addBackRxnNum { static_cast<int>(backRxns.size()) }; //this is the origin number of backRxns for add case
    int addCreateDestructRxnNum { static_cast<int>(createDestructRxns.size()) }; //this is the origin number of createDestructRxns for add case
    int addMolTemplateListNum { static_cast<int>(molTemplateList.size()) };

    std::ifstream inputFile { fileName };
    if (!inputFile) {
        std::cerr << "Add reaction file cannot be opened. Exiting..." << std::endl;
        exit(1);
    }

    std::vector<std::string> providedObs {};
    // the start and end blocks are hardcoded since regex has proved itself woefully non-portable,
    // especially with the Intel compiler (icpc)
    while (!inputFile.eof()) {
        std::string line;
        getline(inputFile, line);
        std::string tmpLine { create_tmp_line(line) };

        if (skipLine(tmpLine)) {
            continue;
        } else if (tmpLine == "startparameters") {
            // read in parameters
            std::cout << "Parsing simulation parameters..." << '\n';
            params.parse_paramFile(inputFile);
        } else if (tmpLine == "startboundaries") {
            // read in boundaries
            std::cout << "Parsing simulation boundary conditions..." << '\n';
            //params.parse_paramFile(inputFile);
            while (getline(inputFile, line)) {
                tmpLine = create_tmp_line(line);
                if (skipLine(tmpLine))
                    continue;
                else if (tmpLine == "endboundaries")
                    break;
                else {
                    line.erase(std::remove_if(line.begin(), line.end(), [](unsigned char x) { return std::isspace(x); }), line.end());
                    std::transform(line.begin(), line.end(), line.begin(), ::tolower); //make all lower case

                    // if the line starts with a comment, skip it
                    if (line[0] == '#')
                        continue;
                    else
                        remove_comment(line);

                    bool gotValue { false };
                    std::string buffer;
                    for (auto lineItr = line.begin(); lineItr != line.end(); ++lineItr) {
                        if (std::isalnum(*lineItr))
                            buffer += std::tolower(static_cast<char>(*lineItr));
                        else if (*lineItr == '=') {
                            auto keyFind = bcKeywords.find(buffer);
                            line.erase(line.begin(), lineItr + 1); // + 1 removes the '=' sign. could make this erase(remove_if)
                            // find the value type from the keyword and then set that parameter
                            if (keyFind != bcKeywords.end()) {
                                // std::cout << "Keyword found: " << keyFind->first << '\n';
                                //this->set_value(line, keyFind->second);
                                membraneObject.set_value_BC(line, keyFind->second);
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
        } else if (tmpLine == "startmolecules") {
            // get the molecule names from the parameter input file
            std::cout << "Parsing mol file information..." << '\n';
            std::vector<std::string> providedMols {};
            std::vector<std::string> providedMols_temp {};
            std::vector<MolTemplate> molTemplateList_temp {};
            std::vector<ParsedMolNumState> providedNumState {};
            std::vector<int> providedNums {};
            while (getline(inputFile, line)) {
                tmpLine = create_tmp_line(line);
                if (skipLine(tmpLine))
                    continue;
                else if (tmpLine == "endmolecules")
                    break;
                else {
                    line.erase(std::remove_if(line.begin(), line.end(), [](unsigned char x) { return std::isspace(x); }), line.end());
                    remove_comment(line);
                    bool gotValue { false };
                    std::string buffer;
                    for (auto lineItr = line.begin(); lineItr != line.end(); ++lineItr) {
                        if (std::isalnum(*lineItr))
                            buffer += static_cast<char>(*lineItr);
                        else if (*lineItr == ':') {
                            providedMols.emplace_back(buffer);
                            line.erase(line.begin(), lineItr + 1); // + 1 removes the ':' sign. could make this erase(remove_if)
                            // find the value type from the keyword and then set that parameter
                            ParsedMolNumState tmpMolNumState {};
                            tmpMolNumState = parse_number_bngl(line);
                            providedNums.emplace_back(tmpMolNumState.totalCopyNumbers);
                            providedNumState.emplace_back(tmpMolNumState);
                            gotValue = true;
                            break;
                        }
                    }
                    if (gotValue == false) {
                        std::cout << "Please provide the copy number of each molecule in INP file in this format-- molName:100" << std::endl;
                        exit(0);
                    }
                }
            }

            // parse the .mol files for each molecule
            for (int providedMolIndex = 0; providedMolIndex < providedMols.size(); providedMolIndex++) {
                molTemplateList.emplace_back(parse_molFile(providedMols[providedMolIndex])); // do the actual parsing
                molTemplateList.back().molTypeIndex = &providedMols[providedMolIndex] - &providedMols[0] + addMolTemplateListNum;
                molTemplateList.back().copies = providedNums[providedMolIndex];
                molTemplateList.back().startingNumState = providedNumState[providedMolIndex];
                ++MolTemplate::numMolTypes;
                if (molTemplateList.back().isImplicitLipid == false) {
                    params.numTotalComplex += molTemplateList.back().copies;
                    params.numTotalUnits
                        += (molTemplateList.back().copies * (molTemplateList.back().interfaceList.size() + 1));
                    if (molTemplateList.back().isLipid)
                        params.numLipids += molTemplateList.back().copies;
                } else {
                }
                MolTemplate::numEachMolType.emplace_back(0);
            }

            params.numMolTypes = MolTemplate::numMolTypes;

            // Give a buffer for the trajectory file (since VMD can't understand trajectories with frames that have
            // different numbers of atoms)
            // TODO: Think about the size of buffer
            //            params.numTotalUnits += params.numTotalUnits / 100;

            // Make sure the writing of restart file and timestep information is in unison
            if (params.restartWrite % params.timeWrite != 0) {
                int tmp = int(params.restartWrite / params.timeWrite);
                if (tmp == 0)
                    tmp = 1;
                params.restartWrite = tmp * params.timeWrite;
            }
            determine_shape_molecule(molTemplateList);
        } else if (tmpLine == "startreactions") {
            // TODO: make sure this can only be read in after molecule templates
            // index of the last template. should. starts out as number of interface states
            int totSpecies { (molTemplateList.back().interfaceList.back().stateList.back().index) }; // last iface index
            int numProvidedRxns { 0 };

            auto linePos = inputFile.tellg();
            std::cout << "Parsing reactions...\n";
            while (getline(inputFile, line)) {
                tmpLine = create_tmp_line(line);
                // if the line is an ignored line, ignore it. if it's an end of block line, break
                if (skipLine(tmpLine)) {
                    linePos = inputFile.tellg();
                    continue;
                } else if (tmpLine == "endreactions")
                    break;
                else {
                    // if it's neither, parse the line
                    inputFile.seekg(linePos);

                    parse_reaction(inputFile, totSpecies, numProvidedRxns, molTemplateList, forwardRxns, backRxns,
                        createDestructRxns, observableList, membraneObject);
                }
                linePos = inputFile.tellg();
            }

            RxnBase::totRxnSpecies = totSpecies + 1; // had to increment this by 1, it was not correct!

        } else if (tmpLine == "startobservables") {
            // TODO: make this use parse_molecule_bngl() and create a vector of ParsedMol
            // Need to change the function to allow for molecules with no explicit interfaces
            std::cout << "Gathering observables...\n";
            while (getline(inputFile, line)) {
                tmpLine = create_tmp_line(line);
                tmpLine.erase(
                    std::remove_if(tmpLine.begin(), tmpLine.end(), [](unsigned char x) { return std::isspace(x); }),
                    tmpLine.end());
                if (skipLine(tmpLine)) {
                    // if the line is a comment line, ignore it.
                    continue;
                } else if (tmpLine == "endobservables") {
                    // if end of parameter block, break out from loop
                    break;
                } else {
                    remove_comment(line);
                    providedObs.emplace_back(line); // save these to parse AFTER reactions
                }
            }
        }
    }

    // TODO: TEMPORARY COUPLED RXN
    for (unsigned rxnItr { 0 }; rxnItr < forwardRxns.size(); ++rxnItr) {
        if (forwardRxns[rxnItr].isCoupled) {
            // if the input is a asRxnIndex
            if (forwardRxns[rxnItr].coupledRxn.absRxnIndex != -1) {
                // for the case destruction coupled to disscociation
                for (const auto& createDestructRxn : createDestructRxns) {
                    if (createDestructRxn.absRxnIndex == forwardRxns[rxnItr].coupledRxn.absRxnIndex) {
                        forwardRxns[rxnItr].coupledRxn.rxnType = ReactionType::destruction;
                        forwardRxns[rxnItr].coupledRxn.relRxnIndex = createDestructRxn.relRxnIndex;
                        forwardRxns[rxnItr].coupledRxn.label = createDestructRxn.rxnLabel;
                        backRxns[forwardRxns[rxnItr].conjBackRxnIndex].coupledRxn = forwardRxns[rxnItr].coupledRxn;
                        break;
                    }
                }

                // for the case uniStateChange coupled to disscociation
                for (const auto& uniStateChangeRxn : forwardRxns) {
                    if (uniStateChangeRxn.rxnType == ReactionType::uniMolStateChange) {
                        if (uniStateChangeRxn.absRxnIndex == forwardRxns[rxnItr].coupledRxn.absRxnIndex) {
                            forwardRxns[rxnItr].coupledRxn.rxnType = ReactionType::uniMolStateChange;
                            forwardRxns[rxnItr].coupledRxn.relRxnIndex = uniStateChangeRxn.relRxnIndex;
                            forwardRxns[rxnItr].coupledRxn.label = uniStateChangeRxn.rxnLabel;
                            backRxns[forwardRxns[rxnItr].conjBackRxnIndex].coupledRxn = forwardRxns[rxnItr].coupledRxn;
                            break;
                        }
                    }
                }
            } else {
                // if the input is a label
                // for the case destruction coupled to disscociation
                for (const auto& createDestructRxn : createDestructRxns) {
                    if (createDestructRxn.rxnLabel == forwardRxns[rxnItr].coupledRxn.label) {
                        forwardRxns[rxnItr].coupledRxn.rxnType = ReactionType::destruction;
                        forwardRxns[rxnItr].coupledRxn.relRxnIndex = createDestructRxn.relRxnIndex;
                        forwardRxns[rxnItr].coupledRxn.absRxnIndex = createDestructRxn.absRxnIndex;
                        backRxns[forwardRxns[rxnItr].conjBackRxnIndex].coupledRxn = forwardRxns[rxnItr].coupledRxn;
                        break;
                    }
                }

                // for the case uniStateChange coupled to disscociation
                for (const auto& uniStateChangeRxn : forwardRxns) {
                    if (uniStateChangeRxn.rxnType == ReactionType::uniMolStateChange) {
                        if (uniStateChangeRxn.rxnLabel == forwardRxns[rxnItr].coupledRxn.label) {
                            forwardRxns[rxnItr].coupledRxn.rxnType = ReactionType::uniMolStateChange;
                            forwardRxns[rxnItr].coupledRxn.relRxnIndex = uniStateChangeRxn.relRxnIndex;
                            forwardRxns[rxnItr].coupledRxn.absRxnIndex = uniStateChangeRxn.absRxnIndex;
                            backRxns[forwardRxns[rxnItr].conjBackRxnIndex].coupledRxn = forwardRxns[rxnItr].coupledRxn;
                            break;
                        }
                    }
                }
            }
        }
    }

    //here wew need to add three new inputs for pop_react_lists
    populate_reaction_lists_for_add(forwardRxns, backRxns, createDestructRxns, molTemplateList, addForwardRxnNum, addBackRxnNum, addCreateDestructRxnNum);

    for (auto& oneTemp : molTemplateList) {
        // here's where we remove duplicate reaction partners
        std::sort(oneTemp.rxnPartners.begin(), oneTemp.rxnPartners.end());
        oneTemp.rxnPartners.erase(
            std::unique(oneTemp.rxnPartners.begin(), oneTemp.rxnPartners.end()), oneTemp.rxnPartners.end());

        for (auto& oneIface : oneTemp.interfaceList) {
            for (auto& oneState : oneIface.stateList) {
                std::sort(oneState.rxnPartners.begin(), oneState.rxnPartners.end());
                oneState.rxnPartners.erase(
                    std::unique(oneState.rxnPartners.begin(), oneState.rxnPartners.end()), oneState.rxnPartners.end());
            }
        }
    }

    params.numTotalSpecies = RxnBase::totRxnSpecies; // total number of all species possible in the system (includes
        // reactant and product states).

    // TODO TEMPORARY
    params.isNonEQ = createDestructRxns.size() > 0;

    std::cout << '\n'
              << "Add input file parsing complete\n";
    // std::cout << "Simulation Parameters\n";
    // params.display();
    // std::cout << "Molecule Information\n";
    // display_all_MolTemplates(molTemplateList);
    // std::cout << "Reactions\n";
    // display_all_reactions(forwardRxns, backRxns, createDestructRxns);
}
