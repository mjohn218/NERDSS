#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include <unordered_set>

bool determine_parent_complex_IL(int pro1Index, int pro2Index, int newComIndex, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, int ILindexMol)
{
    int c1 { moleculeList[pro1Index].myComIndex };
    std::vector<int> origlist = complexList[c1].memberList;
    moleculeList[pro2Index].myComIndex = newComIndex;
    std::unordered_set<int> toDoSet;
    std::unordered_set<int> doneSet;
    toDoSet.insert(pro1Index);
    while(toDoSet.empty() == false){
        auto itr = toDoSet.begin();
        int i = (*itr);
        doneSet.insert(i);
        toDoSet.erase(itr);
        for(auto & p : moleculeList[i].bndpartner){
            if (p == ILindexMol) continue;
            if(p == pro2Index){
                moleculeList[pro1Index].myComIndex = c1;
                moleculeList[pro2Index].myComIndex = c1;
                return true;
            }
            if(doneSet.find(p) == doneSet.end()) toDoSet.insert(p);
        }
    }
    complexList[c1].memberList.clear();
    complexList[c1].memberList.insert(complexList[c1].memberList.end(), doneSet.begin(), doneSet.end());
    complexList[newComIndex].memberList.clear();
    for(auto & m : origlist){
        if(doneSet.find(m) == doneSet.end()) complexList[newComIndex].memberList.push_back(m);
    }
    complexList[newComIndex].index = newComIndex;
    for(auto & m : complexList[c1].memberList){
        moleculeList[m].myComIndex = c1;
    }
    for(auto & m : complexList[newComIndex].memberList){
        moleculeList[m].myComIndex = newComIndex;
    }
    return false;
}
