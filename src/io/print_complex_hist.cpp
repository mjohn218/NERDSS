#include "classes/class_Molecule_Complex.hpp"
#include "io/io.hpp"
#include "tracing.hpp"
#include <iostream>

double print_complex_hist(std::vector<Complex>& complexList, std::ofstream& outfile, int it, Parameters params,
    std::vector<MolTemplate>& molTemplateList, int nImplicitLipids)
{
    // TRACE();
    int i { 0 };
    int j { 0 };
    int size { 0 };
    int p1 { 0 };
    double index { 0 }; // might exceed max integer value.
    /*Rename for simpler typing*/
    int nTypes = params.numMolTypes;

    double mult[nTypes + 1]; // add plus one in case of surface interactions
    int loc { 0 };
    int flag { 0 };

    std::vector<double> assemblylist;
    std::vector<int> histogram;
    std::vector<int> complexrep;
    // cout <<" in calc_complex_hist: "<<nTypes <<std::endl;

    for (j = 0; j < nTypes; j++) {
        mult[j] = 1;
        for (i = 0; i < j; i++)
            mult[j] = mult[j] * (MolTemplate::numEachMolType[i] + 1);
        // cout <<"mult factor for type: "<<j<<" is: "<<mult[j]<<std::endl;
    }
    int jend = nTypes;
    mult[jend] = 1;
    for (i = 0; i < jend; i++)
        mult[jend] = mult[jend]
            * (MolTemplate::numEachMolType[i] + 1); // add this in for possible surface interactions, linksToSurface

    assemblylist.reserve(nTypes);
    complexrep.reserve(nTypes);
    histogram.reserve(nTypes);

    if (complexList.size() == 0) {
        outfile << "iter: " << it << '\n';
        outfile << "NA\n";
        return 0.0;
    }

    // cout <<"Ncomplexes: "<<Nc<<std::endl;
    /*Create the first complex type.*/
    i = 0;
    while (complexList[i].isEmpty)
        i++;
    int i1 = i; // might not be zero if this element is empty

    index = 0;
    for (j = 0; j < nTypes; j++) {
        index += complexList[i].numEachMol[j] * mult[j];
    }
    //    std::cout <<" Complex, LINKS TO SURFACE: "<<complexList[i].linksToSurface<<std::endl;
    index += complexList[i].linksToSurface * mult[jend];
    // new assembly type.

    assemblylist.push_back(index);
    complexrep.push_back(i); // example complex with this composition
    histogram.push_back(1); // created a new assembly type with one so far
    // cout <<"New assembly type: "<<index<<" ffrom complex: "<<i<<" size of assemblylist:
    // "<<assemblylist.size()<<std::endl;

    for (i = i1 + 1; i < complexList.size(); i++) {
        // smoe elements of this array are empty, do not count them
        if (!complexList[i].isEmpty) {

            index = 0;
            for (j = 0; j < nTypes; j++) {
                index += complexList[i].numEachMol[j] * mult[j];
            }
            index += complexList[i].linksToSurface * mult[jend];
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
                //   "<<assemblylist.size()<<std::endl;
            } else {
                histogram[loc] += 1;
            }

        } // done if empty
    } // done over all complexes
    // cout <<" Assembly types: "<<assemblylist.size()<<std::endl;
    // cout <<" histogram of last one: "<<histogram[assemblylist.size()-1]<<std::endl;
    /*Calculated histograms of the assemblies, now write them out.*/
    double meanComplexSize = 0.0; // This will count mean complex size over all complexes >1 protein
    int numComplexTypes = 0;
    int totProteins = 0;
    outfile << "iter: " << it << std::endl;
    for (int a = 0; a < assemblylist.size(); a++) {
        int c1 = complexrep[a];
        if (complexList[c1].memberList.size() == 1) {
            //int proIndex=complexList[c1].memberList[0];
            for (j = 0; j < nTypes; j++) {
                if (complexList[c1].numEachMol[j] != 0 && molTemplateList[j].isImplicitLipid == true) {
                    //instead of saying 1 IL, print out N copies of the IL.
                    histogram[a] = nImplicitLipids;
                }
            }
        }
        // cout <<"Representative complex: "<<c1<<std::endl;
        outfile << histogram[a] << '\t';
        totProteins = 0;
        for (j = 0; j < nTypes; j++) {
            if (complexList[c1].numEachMol[j] != 0) {
                outfile << molTemplateList[j].molName << ": " << complexList[c1].numEachMol[j] << ". ";
                totProteins += complexList[c1].numEachMol[j];
            }
        }
        if (complexList[c1].linksToSurface > 0) {
            outfile << molTemplateList[0].molName
                    << ": " << complexList[c1].linksToSurface << ". ";
        }

        // here we want to record the complex's coord
        //if (histogram[a] == 1) {
        //    outfile << "?" << complexList[c1].comCoord.x << '\t' << complexList[c1].comCoord.y << '\t' << complexList[c1].comCoord.z << "?";
        //}

        if (assemblylist[a] == 0)
            outfile << "PI1: 1. ";
        outfile << std::endl;
        if (totProteins > 1) {
            numComplexTypes += histogram[a];
            meanComplexSize += histogram[a] * totProteins;
        }
    }

    if (meanComplexSize != 0) {
        // this is also = NtotPro_inAssemblies/NAssemblies, so the numerator is all proteins that are not monomers
        meanComplexSize = meanComplexSize / (1.0 * numComplexTypes);
    }

    return meanComplexSize;
}
