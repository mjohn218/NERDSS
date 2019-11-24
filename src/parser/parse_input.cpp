#include "io/io.hpp"
#include "parser/parser_functions.hpp"
#include "classes/class_Membrane.hpp"

// this is so we can compare parsed keywords with the enumerations (we need them as strings)
std::map<const std::string, BoundaryKeyword> bcKeywords = {
    { "implicitlipid", BoundaryKeyword::implicitLipid }, { "waterbox", BoundaryKeyword::waterBox},
    { "xbctype", BoundaryKeyword::xBCtype },{ "ybctype", BoundaryKeyword::yBCtype },
{ "zbctype", BoundaryKeyword::zBCtype }, { "sphere", BoundaryKeyword::sphere }
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
            break;
        case 2:
            this->xBCtype = value;
            break;
        case 3:
            this->yBCtype = value;
            break;
        case 4:
            this->zBCtype = value;
            break;
        case 5:
            this->sphereRadius = stod(value);
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
		 std::vector<CreateDestructRxn>& createDestructRxns, std::vector<MolTemplate>& molTemplateList, Membrane &membraneObject)
{
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
            std::cout << bon << "Parsing simulation parameters..." << boff << '\n';
            params.parse_paramFile(inputFile);
        } else if (tmpLine == "startboundaries") {
            // read in boundaries
            std::cout << bon << "Parsing simulation boundary conditions..." << boff << '\n';
            //params.parse_paramFile(inputFile);
	         while (getline(inputFile, line)) {
                tmpLine = create_tmp_line(line);
                if (skipLine(tmpLine))
                    continue;
                else if (tmpLine == "endboundaries")
                    break;
                else {
		              line.erase(std::remove_if(line.begin(), line.end(), [](unsigned char x) { return std::isspace(x); }),line.end());
		              std::transform(line.begin(), line.end(), line.begin(), ::tolower);//make all lower case
		    
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
				                    std::cout << "Keyword found: " << keyFind->first << '\n';
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
            std::cout << linebreak << bon << "Parsing mol file information..." << boff << '\n';
            std::vector<std::string> providedMols {};
            std::vector<std::string> providedMols_temp {};
            std::vector<MolTemplate> molTemplateList_temp {};
            while (getline(inputFile, line)) {
                tmpLine = create_tmp_line(line);
                if (skipLine(tmpLine))
                    continue;
                else if (tmpLine == "endmolecules")
                    break;
                else {
                    line.erase(std::remove_if(line.begin(), line.end(), [](unsigned char x) { return std::isspace(x); }),line.end());
                    remove_comment(line);
                    providedMols.emplace_back(line); 
                }
            }
            
            // parse the .mol files for each molecule
            for (auto& mol : providedMols) {
                molTemplateList.emplace_back(parse_molFile(mol)); // do the actual parsing
                molTemplateList.back().molTypeIndex = &mol - &providedMols[0];
                ++MolTemplate::numMolTypes;
                if (molTemplateList.back().isImplicitLipid == false) {
                    params.numTotalComplex += molTemplateList.back().copies;
                    params.numTotalUnits
                        += (molTemplateList.back().copies * (molTemplateList.back().interfaceList.size() + 1));
                    if (molTemplateList.back().isLipid)
                        params.numLipids += molTemplateList.back().copies;
                }
                else {
                    
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
                        createDestructRxns, observableList);
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
            for (const auto& createDestructRxn : createDestructRxns) {
                if (createDestructRxn.absRxnIndex == forwardRxns[rxnItr].coupledRxn.absRxnIndex) {
                    forwardRxns[rxnItr].coupledRxn.rxnType = ReactionType::destruction;
                    forwardRxns[rxnItr].coupledRxn.relRxnIndex = createDestructRxn.relRxnIndex;
                    backRxns[forwardRxns[rxnItr].conjBackRxnIndex].coupledRxn = forwardRxns[rxnItr].coupledRxn;
                    rxnItr = forwardRxns.size();
                    break;
                }
            }
        }
    }

    // TODO: PARSE OBSERVABLES HERE
    //    for (const auto& obs : providedObs)
    //        observableList.emplace_back(parse_observable(obs, molTemplateList, forwardRxns, createDestructRxns));

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

    std::cout << '\n' << llinebreak << bon << "Input file parsing complete\n" << boff << llinebreak;
    std::cout << bon << "Simulation Parameters\n" << boff << linebreak;
    params.display();
    std::cout << llinebreak << bon << "Molecule Information\n" << boff << linebreak;
    display_all_MolTemplates(molTemplateList);
    std::cout << llinebreak << bon << "Reactions\n" << boff << linebreak;
    display_all_reactions(forwardRxns, backRxns, createDestructRxns);
    std::cout << llinebreak;
}
