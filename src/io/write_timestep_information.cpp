#include "io/io.hpp"
#include "tracing.hpp"
#include <chrono>
#include <ctime>

void write_timestep_information(long long int simItr, std::ofstream& outFile, std::ofstream& molecTypesFile,
    std::ofstream& textTimeStatFile, const Parameters& params, std::vector<std::vector<int>>& molecTypesList,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, std::vector<MolTemplate>& molTemplateList)
{
    // TRACE();
    Complex::currNumberComTypes = 0; // reset and then recount complex types
    int HEIGHT { Complex::numberOfComplexes };
    int WIDTH = 2 + 2 * molTemplateList.size(); // max number of different molecules in a complex is Nprotypes

    int** adjmat;
    // Allocate memory for complex info
    adjmat = new int*[HEIGHT];
    for (int i = 0; i < HEIGHT; i++)
        adjmat[i] = new int[WIDTH];

    for (int i = 0; i < HEIGHT; i++) {
        for (int j = 0; j < WIDTH; j++) {
            adjmat[i][j] = 0;
        }
    }

    int** complexstat; // define and allocate unique complex ids/positions
    complexstat = new int*[Complex::numberOfComplexes];
    for (int i = 0; i < Complex::numberOfComplexes; i++)
        complexstat[i] = new int[2];

    for (int i = 0; i < Complex::numberOfComplexes; i++) {
        for (int j = 0; j < 2; j++) {
            complexstat[i][j] = 0;
        }
    }

    for (int i = 0; i < HEIGHT; i++) {
        bool isDupeComplex { false };
        if (Complex::currNumberComTypes == 0) {
            complexstat[Complex::currNumberComTypes][0] = i;
            Complex::currNumberComTypes += 1;
        }

        int numUniqueMols { 0 };
        int uniqueComID { 1 };
        for (int j { 0 }; j < complexList[i].memberList.size(); j++) {
            int existingMolIndex { 0 };
            isDupeComplex = false;
            existingMolIndex = 0;
            for (int k = 0; k < numUniqueMols; k++) {
                if (adjmat[i][2 * (k + 1)] == moleculeList[i].molTypeIndex) {
                    isDupeComplex = true;
                    existingMolIndex = k;
                }
            }

            if (isDupeComplex) {
                adjmat[i][2 * (existingMolIndex + 1) + 1] += 1;
            } else {
                numUniqueMols += 1;
                adjmat[i][2 * numUniqueMols] = moleculeList[i].molTypeIndex;
                adjmat[i][2 * numUniqueMols + 1] += 1;
            }
        }

        adjmat[i][1] = numUniqueMols;
        for (int k = 1; k < 2 * (numUniqueMols) + 2; k++) { // calculate unique compid
            uniqueComID = uniqueComID * (adjmat[i][k] + 1);
        }
        adjmat[i][0] = uniqueComID;

        isDupeComplex = false;
        for (int j = 0; j < Complex::currNumberComTypes; j++) {
            if (adjmat[i][0] == adjmat[complexstat[j][0]][0]) {
                isDupeComplex = true;
                break;
            }
        }
        if (!isDupeComplex) {
            complexstat[Complex::currNumberComTypes][0] = i;
            Complex::currNumberComTypes += 1;
        }
    }

    outFile << simItr << "\t" << simItr * params.timeStep * 1e-6 << "\t";

    for (int i = 0; i < Complex::currNumberComTypes; i++) {

        bool uniqmolexists { false };
        int compindex = complexstat[i][0];
        int numuniqmol = adjmat[compindex][1];

        for (int j = 0; j < Complex::currNumberMolTypes; j++) {
            if (adjmat[compindex][0] == molecTypesList[j][0]) {
                uniqmolexists = 1;
                break;
            }
        }

        if (uniqmolexists == 0) { // this is for keeping track of unique types of molecules produced

            for (int k = 0; k < 2 * adjmat[compindex][1] + 2; k++) {
                molecTypesList[Complex::currNumberMolTypes][k] = adjmat[compindex][k];
            }
            for (int k = 1; k < numuniqmol + 1; k++) {
                molecTypesFile << "P";
                molecTypesFile << adjmat[compindex][2 * k];
                molecTypesFile << ":";
                molecTypesFile << adjmat[compindex][2 * k + 1];
            }
            molecTypesFile << '\n';
            Complex::currNumberMolTypes += 1;
        }

        for (int j = 0; j < HEIGHT; j++) {
            if (adjmat[j][0] == adjmat[compindex][0]) {
                complexstat[i][1] += 1;
            }
        }

        outFile << complexstat[i][1];

        for (int j = 0; j < molTemplateList.size(); j++) {
            for (int k = 1; k < numuniqmol + 1; k++) {
                if (adjmat[compindex][2 * k] == j) {
                    textTimeStatFile << "P";
                    textTimeStatFile << adjmat[compindex][2 * k];
                    textTimeStatFile << ":";
                    textTimeStatFile << adjmat[compindex][2 * k + 1];
                }
            }
        }

        outFile << "\t";
        textTimeStatFile << "\t";
    }
    outFile << '\n';
    textTimeStatFile << std::endl;

    // De-Allocate memory to prevent memory leak
    for (int i = 0; i < HEIGHT; i++)
        delete[] adjmat[i];
    delete[] adjmat;

    for (int i = 0; i < Complex::numberOfComplexes; i++)
        delete[] complexstat[i];
    delete[] complexstat;
}
