#include "io/io.hpp"
#include "tracing.hpp"
#include <chrono>
#include <ctime>
#include <iomanip>

void write_transition(double time, std::ofstream& transitionFile, const std::vector<MolTemplate>& molTemplateList)
{
    transitionFile << "time: " << time << std::endl;
    transitionFile << "transion matrix for each mol type: " << std::endl;
    for (auto& molTemp : molTemplateList){
        if(molTemp.countTransition == true){
            transitionFile << molTemp.molName << std::endl;
            for (int indexOne = 0; indexOne < molTemp.transitionMatrixSize; ++indexOne) {
                for (int indexTwo = 0; indexTwo < molTemp.transitionMatrixSize; ++indexTwo) {
                    transitionFile <<' '<< molTemp.transitionMatrix[indexOne][indexTwo];
                }
                transitionFile << std::endl;
            }
        }
    }
    transitionFile << "lifetime for each mol type: " << std::endl;
    for (auto& molTemp : molTemplateList){
        if(molTemp.countTransition == true){
            transitionFile << molTemp.molName << std::endl;
            for (int indexOne = 0; indexOne < molTemp.transitionMatrixSize; ++indexOne) {
                transitionFile << "size of the cluster:" << indexOne+1 << std::endl;
                for (auto& oneTime : molTemp.lifeTime[indexOne]) {
                    transitionFile <<' '<< oneTime;
                }
                transitionFile << std::endl;
            }
        }
    }
}
