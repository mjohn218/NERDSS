#include "classes/class_Molecule_Complex.hpp"
#include "io/io.hpp"
#include "tracing.hpp"
#include <iostream>

using namespace std;

void print_dimers(std::vector<Complex>& complexList, std::ofstream& outfile, int it, Parameters params,
    std::vector<MolTemplate>& molTemplateList)
{
    // TRACE();
    int i, j;
    int size;
    int p1;
    double index; // might exceed max integer value.

    /*Rename for simpler typing*/
    int nTypes = params.numMolTypes;

    double mult[nTypes];
    int loc, flag;

    std::vector<double> assemblylist;
    std::vector<int> histogram;
    std::vector<int> complexrep;
    for (j = 0; j < nTypes; j++) {
        mult[j] = 1;
        for (i = 0; i < j; i++)
            mult[j] = mult[j] * (MolTemplate::numEachMolType[i] + 1);
        // cout <<"mult factor for type: "<<j<<" is: "<<mult[j]<<endl;
    }
    assemblylist.reserve(nTypes);
    complexrep.reserve(nTypes);
    histogram.reserve(nTypes);

    if (complexList.size() == 0) {
        for (unsigned i { 0 }; i < molTemplateList.size(); ++i)
            outfile << "0\t0" << '\t';
        return;
    }

    // cout <<"Ncomplexes: "<<Nc<<endl;
    //Create the first complex type.
    i = 0;
    while (complexList[i].isEmpty == true)
        i++;
    int i1 = i; // might not be zero if this element is empty

    index = 0;
    for (j = 0; j < nTypes; j++) {
        index += complexList[i].numEachMol[j] * mult[j];
    }

    // new assembly type.

    assemblylist.push_back(index);
    complexrep.push_back(i); // example complex with this composition
    histogram.push_back(1); // created a new assembly type with one so far
    // cout <<"New assembly type: "<<index<<" ffrom complex: "<<i<<" size of assemblylist: "<<assemblylist.size()<<endl;

    for (i = i1 + 1; i < complexList.size(); i++) {
        if (complexList[i].isEmpty == false) {

            index = 0;
            for (j = 0; j < nTypes; j++) {
                index += complexList[i].numEachMol[j] * mult[j];
            }
            flag = 0;
            for (int a = 0; a < assemblylist.size(); a++) {
                if (index == assemblylist[a]) {
                    loc = a;
                    flag = 1;
                    a = assemblylist.size(); // break from loop, found your index already.
                }
            }
            if (flag == 0) {
                // new assembly type.

                assemblylist.push_back(index);
                complexrep.push_back(i); // example complex with this composition
                histogram.push_back(1); // created a new assembly type with one so far
                //   cout <<"New assembly type: "<<index<<" ffrom complex: "<<i<<" size of assemblylist:
                //   "<<assemblylist.size()<<endl;
            } else {
                histogram[loc] += 1;
            }

        } // some elements are empty, skip them
    }
    // cout <<" Assembly types: "<<assemblylist.size()<<endl;
    // cout <<" histogram of last one: "<<histogram[assemblylist.size()-1]<<endl;
    /*Calculated histograms of the assemblies, now write them out.*/

    int monomers[nTypes];
    int dimers[nTypes];
    for (int a = 0; a < nTypes; a++) {
        monomers[a] = 0;
        dimers[a] = 0;
    }
    for (int a = 0; a < assemblylist.size(); a++) {
        int c1 = complexrep[a];
        if (complexList[c1].memberList.size() == 1) {
            /*This is a monomer*/
            // outfile<<histogram[a]<<'\t';
            for (j = 0; j < nTypes; j++) {
                if (complexList[c1].numEachMol[j] == 1) {
                    monomers[j] = histogram[a];
                }
            }
        } else if (complexList[c1].memberList.size() == 2) {
            /*This is a dimer, either homo (2 copies) or hetero (1 copy)*/
            for (j = 0; j < nTypes; j++) {
                if (complexList[c1].numEachMol[j] == 1 || complexList[c1].numEachMol[j] == 2) {
                    dimers[j] += histogram[a]; // SUM OVER ALL POSSIBLE DIMERS FOR PROTEIN J
                }
            }
        }
        // outfile <<infofilenames[j].substr(0,3)<<": "<<complexList[c1].numEachMol[j]<<". ";

        //	  outfile<<endl;

    } // loop over all assemblies in the system
    outfile << (it - params.itrRestartFrom) * params.timeStep * 1E-6 + params.timeRestartFrom << '\t';
    for (j = 0; j < nTypes; j++) {
        outfile << monomers[j] << '\t' << dimers[j] << '\t'; // endl;
    }
    outfile << endl;
}
