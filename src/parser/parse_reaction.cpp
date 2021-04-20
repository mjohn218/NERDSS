#include "error_handling.hpp"
#include "io/io.hpp"
#include "parser/parser_functions.hpp"
#include <cmath>

void parse_reaction(std::ifstream& reactionFile, int& totSpecies, int& numProvidedRxns,
    std::vector<MolTemplate>& molTemplateList, std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns,
    std::vector<CreateDestructRxn>& createDestructRxns, std::map<std::string, int>& observablesList, Membrane& membraneObject)
{
    /* NOTE: need to edit both this and enum class RxnKeyword if you want to add keywords  */
    std::map<const std::string, RxnKeyword> rxnKeywords = { { "onrate3dka", RxnKeyword::onRate3Dka }, { "onrate3dmacro", RxnKeyword::onRate3DMacro },
        { "offratekb", RxnKeyword::offRatekb }, { "offratemacro", RxnKeyword::offRateMacro }, { "norm1", RxnKeyword::norm1 }, { "norm2", RxnKeyword::norm2 },
        { "sigma", RxnKeyword::sigma }, { "assocangles", RxnKeyword::assocAngles }, { "onmem", RxnKeyword::onMem },
        { "rate", RxnKeyword::rate }, { "iscoupled", RxnKeyword::isCoupled }, { "isobserved", RxnKeyword::isObserved },
        { "observelabel", RxnKeyword::observeLabel }, { "bindradsamecom", RxnKeyword::bindRadSameCom },
        { "irrevringclosure", RxnKeyword::irrevRingClosure },
        { "creationradius", RxnKeyword::creationRadius }, { "loopcoopfactor", RxnKeyword::loopCoopFactor },
        { "length3dto2d", RxnKeyword::length3Dto2D }, { "rxnlabel", RxnKeyword::rxnLabel }, { "coupledrxnlabel", RxnKeyword::coupledRxnLabel }, { "kcat", RxnKeyword::kcat }, { "excludevolumebound", RxnKeyword::excludeVolumeBound } };

    // Parse reaction
    ParsedRxn parsedRxn;

    // read the reaction itself
    std::string reaction;
    getline(reactionFile, reaction);
    remove_comment(reaction);
    std::cout << "Parsing new reaction with form:\n"
              << reaction << '\n';

    // temporary storage containers
    std::string productSide;
    std::string reactantSide;
    std::vector<std::pair<std::string, int>> reactantSpecies;
    std::vector<std::pair<std::string, int>> productSpecies;

    { // break into reactant and product side using <-> or ->
        // remove spaces
        reaction.erase(
            std::remove_if(reaction.begin(), reaction.end(), [](unsigned char x) { return std::isspace(x); }),
            reaction.end());
        if (reaction.find("<->") != std::string::npos) {
            std::cout << "Reversible reaction." << std::endl;
            size_t position = reaction.find("<->");
            reactantSide = reaction.substr(0, position);
            productSide = reaction.substr(position + 3, std::string::npos); // +1 is to remove the delimiter
            parsedRxn.isReversible = true;
        } else if (reaction.find("->") != std::string::npos) {
            std::cout << "Irreversible reaction." << std::endl;
            size_t position = reaction.find("->");
            reactantSide = reaction.substr(0, position);
            productSide = reaction.substr(position + 2, std::string::npos); // +2 is to remove delimiter
        } else
            invalid_rxn(std::string("Missing reaction arrow."), __func__, __LINE__);
    }

    { // break into species based on '+'
        size_t position { 0 };
        int speciesIndex { 0 };
        while ((position = reactantSide.find('+')) != std::string::npos) {
            reactantSpecies.emplace_back(reactantSide.substr(0, position), speciesIndex);
            reactantSide.erase(0, position + 1);
            ++speciesIndex;
        }
        // std::string::substr(position, position + length of substr). I think this might accidentally work
        reactantSpecies.emplace_back(reactantSide.substr(0, std::string::npos), speciesIndex);

        speciesIndex = 0;
        while ((position = productSide.find('+')) != std::string::npos) {
            productSpecies.emplace_back(productSide.substr(0, position), speciesIndex);
            productSide.erase(0, position + 1);
            ++speciesIndex;
        }
        productSpecies.emplace_back(productSide.substr(0, std::string::npos), speciesIndex);
    }

    // check to make sure the reactants side only has two molecules
    if (reactantSpecies.size() > 2) {
        std::cerr << "Error, invalid reaction.\n";
        std::cerr << "Reaction " << reaction << " has more than two reacting molecules.\n";
        exit(1);
    }

    // Check for the reaction type
    // TODO: make this account for more than one species annihilating?
    if ((reactantSpecies.size() == 1) && (productSpecies.size() == 1)) {
        /* searches for which side is the null side
         */
        std::string tmpLSide { reactantSpecies[0].first };
        std::string tmpRSide { productSpecies[0].first };
        std::transform(tmpLSide.begin(), tmpLSide.end(), tmpLSide.begin(), ::tolower);
        std::transform(tmpRSide.begin(), tmpRSide.end(), tmpRSide.begin(), ::tolower);
        bool lside { tmpLSide == "0" || tmpLSide == "null" };
        bool rside { tmpRSide == "0" || tmpRSide == "null" };
        if (!lside && !rside) {
            parsedRxn.rxnType = ReactionType::uniMolStateChange;
            std::cout << "UniMolStateChange reaction detected\n";
        } else if (lside && !rside) {
            // Creation reaction
            parsedRxn.rxnType = ReactionType::zerothOrderCreation;
            std::cout << "ZerothOrderCreation reaction detected\n";
        } else if (rside && !lside) {
            // Destruction reaction
            parsedRxn.rxnType = ReactionType::destruction;
            std::cout << "Destruction reaction detected\n";
        } else {
            std::cerr << "FATAL ERROR: Ccannot determine reaction type. Please check before moving on.\n";
            exit(1);
        }
    } else if ((reactantSpecies.size()) == 1 && (productSpecies.size() == 2)) {
        parsedRxn.rxnType = ReactionType::uniMolCreation;
        std::cout << "UniMolCreation reaction detected.\n";
    } else {
        parsedRxn.rxnType = ReactionType::bimolecular;
        std::cout << "Bimolecular reaction detected.\n";
    }

    // the following two large blocks of text parsing can't be split off into one subroutine, because of the use of
    // the local structs ParsedMol and ParsedRxn
    if (parsedRxn.rxnType != ReactionType::zerothOrderCreation) {
        // TODO: need to pass some integer to tell ParsedMol which species the molecule is in (via speciesIndex)
        std::cout << "Parsing reactants...\n";
        std::vector<std::pair<std::string, int>> reactantMols;
        for (auto& specie : reactantSpecies) {
            size_t position { 0 };
            std::string tmpSpecie { specie.first }; // so i don't erase the real thing
            while ((position = tmpSpecie.find('.')) != std::string::npos) {
                reactantMols.emplace_back(tmpSpecie.substr(0, position), specie.second);
                tmpSpecie.erase(0, position + 1);
            }
            reactantMols.emplace_back(
                tmpSpecie.substr(0, std::string::npos), specie.second); // so it actually includes the last molecule
            // if there was only one molecule, make sure it makes its way into the vector
            if (reactantMols.empty() && !specie.first.empty())
                reactantMols.push_back(specie);
        }

        for (auto& oneMol : reactantMols) {
            parsedRxn.reactantList.push_back(parse_molecule_bngl(totSpecies, false, oneMol));
            parsedRxn.reactantList.back().set_molTypeIndex(molTemplateList);
        }
    }

    if (parsedRxn.rxnType != ReactionType::destruction) {
        std::cout << "Parsing products...\n";
        std::vector<std::pair<std::string, int>> productMols;
        int speciesIndex { 0 };
        for (auto& specie : productSpecies) {
            size_t position { 0 };
            std::string tmpSpecie { specie.first }; // so i don't erase the real thing
            if (tmpSpecie.find('.') != std::string::npos) {
                while ((position = tmpSpecie.find('.')) != std::string::npos) {
                    productMols.emplace_back(tmpSpecie.substr(0, position), speciesIndex);
                    tmpSpecie.erase(0, position + 1);
                    ++speciesIndex;
                }
                productMols.emplace_back(
                    tmpSpecie.substr(0, std::string::npos), speciesIndex); // so it actually includes the last molecule
            } else {
                // this is to account for bimolecular state changes that have no bound species
                productMols.emplace_back(
                    tmpSpecie.substr(0, std::string::npos), speciesIndex); // so it actually includes the last molecule
                ++speciesIndex;
            }
            // if there was only one molecule, make sure it makes its way into the vector
            if (productMols.empty() && !specie.first.empty())
                productMols.push_back(specie);
        }

        for (auto& oneMol : productMols) {
            parsedRxn.productList.push_back(parse_molecule_bngl(totSpecies, true, oneMol));
            parsedRxn.productList.back().set_molTypeIndex(molTemplateList);
        }
    }

    // Get the reaction parameters
    std::string line;
    std::string tmpLine { line }; // to avoid altering the original line, just in case
    while (line != "endreactions") {
        auto initialPos = reactionFile.tellg();
        getline(reactionFile, line);
        line.erase(
            std::remove_if(line.begin(), line.end(), [](unsigned char x) { return std::isspace(x); }), line.end());

        if ((line.find('#') != std::string::npos) && line[0] != '#')
            remove_comment(line);

        if (line[0] == '#') {
            initialPos = reactionFile.tellg();
            continue;
        } else if ((line.find("->") != std::string::npos) || (line.find("<->") != std::string::npos)) {
            // if the next line contains a reaction, seek back to where you started in the file and stop looking for
            // parameter keywords
            reactionFile.seekg(initialPos);
            std::cout << "Finished reading reaction parameters, constructing reaction.\n";
            break;
        }

        std::string buffer;
        for (auto lineItr = line.begin(); lineItr != line.end(); ++lineItr) {
            if (std::isalnum(*lineItr))
                buffer += std::tolower(static_cast<char>(*lineItr));
            else if (*lineItr == '=') {
                auto keyFind = rxnKeywords.find(buffer);
                if (keyFind == rxnKeywords.end()) {
                    std::cerr << buffer + " is an invalid argument for the reactions block. Exiting...";
                    exit(1);
                }

                line.erase(line.begin(), lineItr + 1); // hard coded for character length '='
                parsedRxn.set_value(line, keyFind->second);
                break;
            }
        }
    }

    parsedRxn.norm1.calc_magnitude();
    parsedRxn.norm2.calc_magnitude();

    // get the indices of each iface in the reactants and products
    for (auto& reactant : parsedRxn.reactantList) {
        determine_iface_indices(
            &reactant - &parsedRxn.reactantList[0], totSpecies, reactant, parsedRxn, forwardRxns, molTemplateList);
        check_for_valid_states(&reactant - &parsedRxn.reactantList[0], reactant, parsedRxn, molTemplateList);
    }

    for (auto& product : parsedRxn.productList)
        determine_iface_indices(
            &product - &parsedRxn.productList[0], totSpecies, product, parsedRxn, forwardRxns, molTemplateList);

    // determine the indices for the reactants and products of the real reaction here
    if (parsedRxn.rxnType == ReactionType::bimolecular) {
        parsedRxn.determine_reactants();

        if (parsedRxn.rxnProducts[0].speciesIndex != 0 && parsedRxn.rxnProducts[1].speciesIndex != 1) {
            std::swap(parsedRxn.rxnProducts[0], parsedRxn.rxnProducts[1]);
        }
    } else if (parsedRxn.rxnType == ReactionType::uniMolStateChange) {
        parsedRxn.determine_reactants();
    } else if (parsedRxn.rxnType == ReactionType::zerothOrderCreation) {
        parsedRxn.determine_creation_products(molTemplateList);
    } else if (parsedRxn.rxnType == ReactionType::uniMolCreation) {
        parsedRxn.determine_creation_products(molTemplateList);
    } else if (parsedRxn.rxnType == ReactionType::destruction) {
        for (auto& reactant : parsedRxn.reactantList) {
            for (auto& iface : reactant.interfaceList) {
                iface.molTypeIndex = reactant.molTypeIndex;
                parsedRxn.intReactantList.emplace_back(iface.absIndex);
                parsedRxn.rxnReactants.emplace_back(iface);
            }
        }
    }

    // TODO: temporary. add to observables list
    if (parsedRxn.isObserved) {
        if (parsedRxn.rxnType == ReactionType::zerothOrderCreation) {
            int startNum { molTemplateList[parsedRxn.rxnProducts[0].molTypeIndex].copies };
            observablesList.emplace(std::make_pair(parsedRxn.observeLabel, startNum));
        } else if (parsedRxn.rxnType == ReactionType::uniMolCreation) {
            int startNum { molTemplateList[parsedRxn.rxnProducts[1].molTypeIndex].copies };
            observablesList.emplace(std::make_pair(parsedRxn.observeLabel, startNum));
        } else if (parsedRxn.rxnType == ReactionType::destruction) {
            int startNum { molTemplateList[parsedRxn.rxnReactants[0].molTypeIndex].copies };
            observablesList.emplace(std::make_pair(parsedRxn.observeLabel, startNum));
        } else if (parsedRxn.rxnType == ReactionType::uniMolStateChange) {
            int molTypeIndex { parsedRxn.rxnReactants[0].molTypeIndex };
            int relIfaceIndex { parsedRxn.rxnReactants[0].relIndex };
            if (parsedRxn.rxnProducts[0].state
                == molTemplateList[molTypeIndex].interfaceList[relIfaceIndex].stateList[0].iden) {
                int startNum { molTemplateList[molTypeIndex].copies };
                observablesList.emplace(std::make_pair(parsedRxn.observeLabel, startNum));
            } else
                observablesList.emplace(std::make_pair(parsedRxn.observeLabel, 0));
        } else {
            observablesList.emplace(std::make_pair(parsedRxn.observeLabel, 0));
        }
    }

    if (parsedRxn.rxnType == ReactionType::bimolecular) {
        int react1ProType { parsedRxn.rxnReactants[0].molTypeIndex };
        std::string react1 { molTemplateList[react1ProType].molName + '(' + parsedRxn.rxnReactants[0].ifaceAndStateName
            + "!1)" };
        int react2ProType { parsedRxn.rxnReactants[1].molTypeIndex };
        std::string react2 { molTemplateList[react2ProType].molName + '(' + parsedRxn.rxnReactants[1].ifaceAndStateName
            + "!1)" };
        parsedRxn.productName = react1 + '.' + react2;
    }
    if (parsedRxn.rxnType == ReactionType::biMolStateChange) {
        parsedRxn.productName = "twoStates";
    }
    if (parsedRxn.rxnType == ReactionType::uniMolStateChange) {
        parsedRxn.productName = "oneStates";
    }
    if (parsedRxn.rxnType == ReactionType::zerothOrderCreation) {
        parsedRxn.productName = "zerothOrderCreation";
    }
    if (parsedRxn.rxnType == ReactionType::destruction) {
        parsedRxn.productName = "destruction";
    }
    if (parsedRxn.rxnType == ReactionType::uniMolCreation) {
        parsedRxn.productName = "uniMolCreation";
    }

    if (parsedRxn.rxnType == ReactionType::bimolecular || parsedRxn.rxnType == ReactionType::biMolStateChange) {
        /*for all bimolecular reactions, they might take place in 2D, need to assign 3Dto2D length. By default,
	  if it is not read in, it will be set to 2*sigma*/

        if (parsedRxn.length3Dto2D == -1)
            parsedRxn.length3Dto2D = 2.0 * parsedRxn.bindRadius;
    }

    if (parsedRxn.rxnType == ReactionType::bimolecular || parsedRxn.rxnType == ReactionType::biMolStateChange) {
        // set the micro rate according to the macro rate
        // for 3D reaction, Dz != 0 for each reactant; for 3D->2D, Dz == 0 for one reactant. determine this first
        bool to2D = false;
        bool is2D = true;
        bool isSelf = false; // is a A(a)+A(a) -> A(a!).A(a!) ?
        for (auto& reactant : parsedRxn.reactantList) {
            if (std::abs(molTemplateList[reactant.molTypeIndex].D.z - 0.0) < 1E-10) {
                to2D = true;
            } else {
                is2D = false;
            }
        }
        if (parsedRxn.rxnType == ReactionType::bimolecular) {
            if (parsedRxn.isSymmetric == true) {
                isSelf = true; // is a A(a)+A(a) -> A(a!).A(a!)
            }
        }

        if (isSelf == true) { //Self reaction
            if (std::isnan(parsedRxn.onRate3Dka) == true) {
                if (std::isnan(parsedRxn.onRate3DMacro) == false) {
                    //convert macro rate to micro rate, uM-1s-1 = 1/0.602214076 nm3/us, get Dtot first
                    double Dtot { 0.0 };
                    for (auto& reactant : parsedRxn.reactantList) {
                        Dtot += (1 / 3.0 * molTemplateList[reactant.molTypeIndex].D.x + 1 / 3.0 * molTemplateList[reactant.molTypeIndex].D.y + 1 / 3.0 * molTemplateList[reactant.molTypeIndex].D.z);
                    }
                    parsedRxn.onRate3Dka = 0.5 / (1.0 / (2.0 * parsedRxn.onRate3DMacro / 0.602214076) - 1.0 / (4.0 * M_PI * parsedRxn.bindRadius * Dtot));
                }
                if (std::isnan(parsedRxn.offRatekb) == true && parsedRxn.isReversible == true) {
                    if (std::isnan(parsedRxn.offRateMacro) == false) {
                        parsedRxn.offRatekb = parsedRxn.offRateMacro * parsedRxn.onRate3Dka * 0.602214706 / parsedRxn.onRate3DMacro;
                    }
                }
            }
        } else {
            if (is2D == true) { //2D reaction
                if (std::isnan(parsedRxn.onRate3Dka) == true) {
                    if (std::isnan(parsedRxn.onRate3DMacro) == false) {
                        //convert macro rate to micro rate, uM-1s-1 = 1/0.602214076 nm3/us, get Dtot first
                        double Dtot { 0.0 };
                        for (auto& reactant : parsedRxn.reactantList) {
                            Dtot += (1 / 2.0 * molTemplateList[reactant.molTypeIndex].D.x + 1 / 2.0 * molTemplateList[reactant.molTypeIndex].D.y);
                        }
                        double tempVariable { 0.0 };
                        double maxN { 0.0 };
                        for (auto& reactant : parsedRxn.reactantList) {
                            if (molTemplateList[reactant.molTypeIndex].copies > maxN) {
                                maxN = molTemplateList[reactant.molTypeIndex].copies;
                            }
                        }
                        double area { 0.0 };
                        area = 1.0 * membraneObject.waterBox.x * membraneObject.waterBox.y;
                        double b { 0.0 };
                        double sigma { parsedRxn.bindRadius };
                        b = 2.0 * pow(area / (M_PI * maxN) + parsedRxn.bindRadius * parsedRxn.bindRadius, 0.5);
                        tempVariable = 4.0 * log(b / sigma) / pow(1.0 - pow(sigma / b, 2.0), 2.0) - 2.0 / (1.0 - pow(sigma / b, 2.0)) - 1.0;
                        parsedRxn.onRate3Dka = 2.0 * parsedRxn.bindRadius / (1.0 / (1.0 * parsedRxn.onRate3DMacro / 0.602214076) - 1.0 / (8.0 * M_PI * Dtot) * tempVariable);
                    }
                    if (std::isnan(parsedRxn.offRatekb) == true && parsedRxn.isReversible == true) {
                        if (std::isnan(parsedRxn.offRateMacro) == false) {
                            parsedRxn.offRatekb = parsedRxn.offRateMacro * parsedRxn.onRate3Dka * 0.602214706 / (2.0 * parsedRxn.bindRadius) / parsedRxn.onRate3DMacro;
                        }
                    }
                }
            } else {
                if (to2D == false) { // 3D reaction
                    if (std::isnan(parsedRxn.onRate3Dka) == true) {
                        if (std::isnan(parsedRxn.onRate3DMacro) == false) {
                            //convert macro rate to micro rate, uM-1s-1 = 1/0.602214076 nm3/us, get Dtot first
                            double Dtot { 0.0 };
                            for (auto& reactant : parsedRxn.reactantList) {
                                Dtot += (1 / 3.0 * molTemplateList[reactant.molTypeIndex].D.x + 1 / 3.0 * molTemplateList[reactant.molTypeIndex].D.y + 1 / 3.0 * molTemplateList[reactant.molTypeIndex].D.z);
                            }
                            parsedRxn.onRate3Dka = 1.0 / (1.0 / (parsedRxn.onRate3DMacro / 0.602214076) - 1.0 / (4.0 * M_PI * parsedRxn.bindRadius * Dtot));
                        }
                        if (std::isnan(parsedRxn.offRatekb) == true && parsedRxn.isReversible == true) {
                            if (std::isnan(parsedRxn.offRateMacro) == false) {
                                parsedRxn.offRatekb = parsedRxn.offRateMacro * parsedRxn.onRate3Dka * 0.602214706 / parsedRxn.onRate3DMacro;
                            }
                        }
                    }
                } else { // 3D->2D
                    if (std::isnan(parsedRxn.onRate3Dka) == true) {
                        if (std::isnan(parsedRxn.onRate3DMacro) == false) {
                            //convert macro rate to micro rate, uM-1s-1 = 1/0.602214076 nm3/us, get Dtot first
                            double Dtot { 0.0 };
                            for (auto& reactant : parsedRxn.reactantList) {
                                Dtot += (1 / 3.0 * molTemplateList[reactant.molTypeIndex].D.x + 1 / 3.0 * molTemplateList[reactant.molTypeIndex].D.y + 1 / 3.0 * molTemplateList[reactant.molTypeIndex].D.z);
                            }
                            parsedRxn.onRate3Dka = 0.5 / (1.0 / (2.0 * parsedRxn.onRate3DMacro / 0.602214076) - 1.0 / (4.0 * M_PI * parsedRxn.bindRadius * Dtot));
                        }
                        if (std::isnan(parsedRxn.offRatekb) == true && parsedRxn.isReversible == true) {
                            if (std::isnan(parsedRxn.offRateMacro) == false) {
                                parsedRxn.offRatekb = parsedRxn.offRateMacro * parsedRxn.onRate3Dka * 0.602214706 / parsedRxn.onRate3DMacro;
                            }
                        }
                    }
                }
            }
        }
    }
    if (parsedRxn.isCoupled == true) {
        //create a probability for this reaction to be Michaelis Menten, if kcat exists and kb exists. Otherwise prob=1.
        if (std::isnan(parsedRxn.kcat) == false) {

            double kcat = parsedRxn.kcat; // kcat is not used after this function in the code, only probCoupled is.
            if (std::isnan(parsedRxn.offRatekb) == true)
                parsedRxn.offRatekb = 0; //set this to zero, so dissociation occurs with rate kcat, and coupled Rxn always occurs.
            parsedRxn.coupledRxn.probCoupled = kcat / (kcat + parsedRxn.offRatekb);
            parsedRxn.offRatekb += kcat; // the rate that the bound complex undergoes reaction is due to both catalysis and unbinding. Will then choose which one to perform based on probCoupled.
            std::cout << "Allow Michaelis-Menten for this reaction. Probability of catalysis: kcat/(kcat+offRatekb): " << parsedRxn.coupledRxn.probCoupled << " total rate: " << parsedRxn.offRatekb << std::endl;
        } else {
            std::cout << " No kcat rate specified, will NOT perform coupled reaction! Must specify rate using kcat " << std::endl;
        }
    }
    // Done parsing, set up the reaction(s) and place into their correct vectors
    if (parsedRxn.willBeMultipleRxns) {
        // if there is an interface that is missing a state (on purpose), we need to create a different forward
        // reaction for each of those missing states
        auto newRxns = parsedRxn.split_into_multiple_reactions(totSpecies, forwardRxns, molTemplateList);
        --numProvidedRxns; // first remove the provded reaction that was split
        for (auto& oneRxn : newRxns) {
            // Make sure everything that we need to be defined was provided
            auto status = oneRxn.isComplete(molTemplateList);
            if (!status.first) {
                std::cerr << status.second << " [" << reaction << "].\n";
                exit(1);
            } else
                std::cout << status.second << " [" << reaction << "].\n";

            oneRxn.create_other_iface_lists(molTemplateList);
            // check if this is a repeat reaction and create new RateState if so
            if (!oneRxn.check_for_conditional_rates(totSpecies, forwardRxns, backRxns))
                oneRxn.assemble_reactions(forwardRxns, backRxns, createDestructRxns, molTemplateList);
        }
    } else {
        // if all states are explicit, create the lists of ancillary ifaces
        parsedRxn.create_other_iface_lists(molTemplateList);
        // check if this is a repeat reaction and create new RateState if so
        if (!parsedRxn.check_for_conditional_rates(totSpecies, forwardRxns, backRxns)) {
            // otherwise, just make the one reaction (or two, if its reversible)
            // Make sure everything that we need to be defined was provided
            auto status = parsedRxn.isComplete(molTemplateList);
            if (!status.first) {
                std::cerr << "Error, " << status.second << " [" << reaction << "].\n";
                exit(1);
            } else
                std::cout << status.second << " [" << reaction << "].\n";
            parsedRxn.assemble_reactions(forwardRxns, backRxns, createDestructRxns, molTemplateList);
        }
    }

    // add the MolTemplate to the other MolTemplate's list of reaction partners
    // sort and remove duplicates later
    if (parsedRxn.rxnType == ReactionType::bimolecular || parsedRxn.rxnType == ReactionType::biMolStateChange) {
        MolTemplate& reactTemp1 = molTemplateList[forwardRxns.back().reactantListNew.at(0).molTypeIndex];
        MolTemplate& reactTemp2 = molTemplateList[forwardRxns.back().reactantListNew.at(1).molTypeIndex];
        reactTemp1.rxnPartners.push_back(reactTemp2.molTypeIndex);
        reactTemp2.rxnPartners.push_back(reactTemp1.molTypeIndex);
    }
    std::cout << std::endl;
}
