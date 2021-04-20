#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "tracing.hpp"

bool determine_parent_complex(int pro1Index, int pro2Index, int newComIndex, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList)
{
    // TRACE();
    // TODO: clean this up
    int c1 { moleculeList[pro1Index].myComIndex };
    std::vector<int> origlist = complexList[c1].memberList;
    moleculeList[pro2Index].myComIndex = newComIndex;

    //    bool keepSameComplex{ false };
    // determine if any two molecules are doubly bound
    bool boundTwice { false };
    for (auto& mp : origlist) {
        std::vector<int> tmppartvec = moleculeList[mp].bndpartner;
        if (std::find(tmppartvec.begin(), tmppartvec.end(), pro1Index) != tmppartvec.end()
            && std::find(tmppartvec.begin(), tmppartvec.end(), pro2Index) != tmppartvec.end()) {
            // std::cout << "Molecule " << mp << " is bound to both dissociating molecules. Keeping as one complex.\n";
            boundTwice = true;
            break;
        }
        for (unsigned i { 0 }; i < tmppartvec.size(); i++) {
            for (unsigned j { 0 }; j < tmppartvec.size(); j++) {
                if (i != j && tmppartvec[i] != -1 && tmppartvec[j] != -1) {
                    if (tmppartvec[i] == tmppartvec[j]) {
                        boundTwice = true;
                    }
                }
            }
        }

        // if (boundTwice)
        //     std::cout << "Molecule " << mp << " is doubly bound!" << std::endl;
    }

    if (boundTwice) {
        moleculeList[pro1Index].myComIndex = c1;
        moleculeList[pro2Index].myComIndex = c1;
        return true;
    } else {
        std::vector<int> tmpc1 { pro1Index };
        std::vector<int> tmpc2 { pro2Index };
        // Now check each protein bound to a dissociating protein to determine which complex they'll belong to after
        // dissociation
        std::vector<int> recheck;
        for (auto memMol : complexList[moleculeList[pro1Index].myComIndex].memberList) {
            if (memMol != pro1Index && memMol != pro2Index) {
                int com_flag { 0 };
                for (auto& ppart : moleculeList[memMol].bndpartner) {
                    if (ppart == pro1Index) {
                        tmpc1.push_back(memMol);
                        com_flag++;
                    } else if (ppart == pro2Index) {
                        tmpc2.push_back(memMol);
                        com_flag++;
                    }
                }
                //                if (com_flag > 1 || com_flag == 0)
                if (com_flag != 1)
                    recheck.push_back(memMol);
            }
        }

        // recheck flagged proteins
        { // make sure to only recheck each protein once
            std::sort(recheck.begin(), recheck.end());
            auto last = std::unique(recheck.begin(), recheck.end());
            recheck.erase(last, recheck.end());
        }
        //        std::cout << "HERE DETERMINE PARENT COMPLEX\n";
        // what this SHOULD do is take each protein in the bndlist of the protein to recheck (a)
        // and check each of those proteins' (b) bndlists (c) for target proteins for the dissociating proteins
        // if they're not found, it goes in another layer and looks at the bndlists of (c),
        // and so on
        for (auto i = recheck.begin(); i != recheck.end();) {
            int it { 0 };
            int molIndex { *i };
            int maxit { 1000 };
            int com_flag { 0 };
            int no { 0 };
            bool doneChecking { false };
            // std::cout <<" NEED to check protein: "<<molIndex<<" boundpartners: ";

            std::vector<int> testList { moleculeList[molIndex].bndpartner };
            std::vector<int> checkedList; // these are the moleculeList that have been checked
            // for(int pp=0;pp<testList.size();pp++)
            // 	      std::cout <<testList[pp]<<'\t';
            // 	    std::cout<<std::endl;

            while (!doneChecking) {
                testList.erase(std::remove(testList.begin(), testList.end(), -1), testList.end());
                std::vector<int> testlist2;
                for (unsigned j { 0 }; j < testList.size(); j++) {
                    int ppart { testList[j] };
                    checkedList.push_back(ppart);
                    for (auto& ppart2 : moleculeList[ppart].bndpartner) {
                        testlist2.push_back(ppart2);
                        if (ppart != molIndex) {
                            if (ppart2 == pro1Index) {
                                tmpc1.push_back(molIndex);
                                //std::cout << "base " << molIndex << " belongs to complex c1." << std::endl;
                                com_flag++;
                            } else if (ppart2 == pro2Index) {
                                tmpc2.push_back(molIndex);
                                //std::cout << "base " << molIndex << " belongs to complex c2." << std::endl;
                                com_flag++;
                            }
                        }
                        if (com_flag > 1) {
                            // if both pro1Index and pro2Index are found, com_flag > 1 and its a closed loop
                            moleculeList[pro1Index].myComIndex = c1;
                            moleculeList[pro2Index].myComIndex = c1;
                            return true;
                        } else if (com_flag == 0) {
                            // neither pro1Index or pro2Index found, move to the next set of test subjects
                            //std::cout <<" did not find either for prote "<<molIndex<<std::endl;
                            no++;
                        } else {
                            no = 0;
                            break;
                        }
                    }
                    if (no == 0)
                        break;
                }
                if (no > 0) {
                    // the new test subject list is all the binding partners of the previous set of test subjects
                    testList.clear();
                    for (auto& mp : testlist2) {
                        if (mp != *i && find(checkedList.begin(), checkedList.end(), mp) == checkedList.end()) {
                            testList.push_back(mp);
                        }
                    }
                    std::sort(testList.begin(), testList.end());
                    testList.erase(std::unique(testList.begin(), testList.end()), testList.end());
                    //                                        testList.erase(testList.begin()); // I cannot remember why
                    //                                        this is here
                } else {
                    //std::cout <<" done checking protein: "<<molIndex<<std::endl;
                    doneChecking = true;
                    recheck.erase(i);
                }
                it++;
                if (it > maxit) {
                    std::cerr << "Stuck in infinite loop..." << std::endl;
                    exit(1);
                }
            }
        }

        //        std::cout << "HERE DETERMINE PARENT COMPLEX2\n";
        // Check for shared members between the old complex and new one
        for (auto& mp : tmpc1) {
            if (find(tmpc2.begin(), tmpc2.end(), mp) != tmpc2.end()) {
                moleculeList[pro1Index].myComIndex = c1;
                moleculeList[pro2Index].myComIndex = c1;
                return true;
            }

            for (auto& ppart : moleculeList[mp].bndpartner) {
                if (find(tmpc2.begin(), tmpc2.end(), ppart) != tmpc2.end()) {
                    // std::cout << "Complex is a closed loop.\n";
                    moleculeList[pro1Index].myComIndex = c1;
                    moleculeList[pro2Index].myComIndex = c1;
                    return true;
                }
            }
        }

        if (tmpc1.size() + tmpc2.size() != complexList[c1].memberList.size()) {
            std::cout << " complex sizes don't match the parent! " << tmpc1.size() << ' ' << tmpc2.size() << " original size: " << complexList[c1].memberList.size() << std::endl;
            std::cout << " pros in c1: " << std::endl;
            for (auto memMol : tmpc1)
                std::cout << memMol << '\t';
            std::cout << " pros in c2: " << std::endl;
            for (auto memMol : tmpc2)
                std::cout << memMol << '\t';
            std::cout << "display molecules " << std::endl;
            for (auto mp : tmpc1)
                moleculeList[mp].display_all();
            std::cout << " pros in c2: " << std::endl;
            for (auto mp : tmpc2)
                moleculeList[mp].display_all();

            std::cerr
                << "ERROR: Combined size of dissociated complexes does not match the parent complex. Exiting...\n";
            exit(1);
        }

        complexList[c1].memberList.swap(tmpc1);
        complexList[newComIndex].memberList.swap(tmpc2);
        complexList[newComIndex].index = newComIndex;
        std::sort(complexList[c1].memberList.begin(), complexList[c1].memberList.end());
        std::sort(complexList[newComIndex].memberList.begin(), complexList[newComIndex].memberList.end());

        for (auto& mp : complexList[c1].memberList)
            moleculeList[mp].myComIndex = c1;
        for (auto& mp : complexList[newComIndex].memberList)
            moleculeList[mp].myComIndex = newComIndex;
    }

    return false;
}
