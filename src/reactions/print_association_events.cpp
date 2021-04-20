#include "boundary_conditions/reflect_functions.hpp"
#include "classes/class_Rxns.hpp"
#include "io/io.hpp"
#include "reactions/association/association.hpp"
#include "reactions/association/functions_for_spherical_system.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"
#include <cmath>
#include <iomanip>

void print_association_events(copyCounters& counterArrays, std::ofstream& outfile, int it, Parameters params)
{

    /*Bundle values for n>maxSingles*/
    int spacing = 10;
    int maxSingles = 10; //keep all n integers less than this value.
    int arraySize = counterArrays.eventArraySize; //final array element contains all additions above (arraySize-maxSingles)*spacing
    outfile << "time (s): " << (it - params.itrRestartFrom) * params.timeStep * 1E-6 + params.timeRestartFrom << std::endl;
    std::string text;
    bool start = true;
    for (int i = 0; i < arraySize; i++) {
        if (counterArrays.events3D[i] > 0) {
            if (i == 0)
                text = "3D dimers";
            else if (i < maxSingles)
                text = "3D n=" + std::to_string(i);
            else {
                int a = i - maxSingles + 1;
                int n = a * spacing;
                if (i == arraySize - 1)
                    text = "3D n>" + std::to_string(n); //this is the max index stored.
                else
                    text = "3D n=" + std::to_string(n) + " to " + std::to_string(n + spacing - 1);
            }

            outfile << text << ": " << counterArrays.events3D[i] << '\n';
        }
    }
    //now 3D to 2D
    for (int i = 0; i < arraySize; i++) {
        if (counterArrays.events3Dto2D[i] > 0) {
            if (i == 0)
                text = "3Dto2D dimers";
            else if (i < maxSingles)
                text = "3Dto2D n=" + std::to_string(i);
            else {
                int a = i - maxSingles + 1;
                int n = a * spacing;
                if (i == arraySize - 1)
                    text = "3Dto2D n>" + std::to_string(n); //this is the max index stored.
                else
                    text = "3Dto2D n=" + std::to_string(n) + " to " + std::to_string(n + spacing - 1);
            }

            outfile << text << ": " << counterArrays.events3Dto2D[i] << '\n';
        }
    }
    //now 2D
    for (int i = 0; i < arraySize; i++) {
        if (counterArrays.events2D[i] > 0) {
            if (i == 0)
                text = "2D dimers";
            else if (i < maxSingles)
                text = "2D n=" + std::to_string(i);
            else {
                int a = i - maxSingles + 1;
                int n = a * spacing;
                if (i == arraySize - 1)
                    text = "2D n>" + std::to_string(n); //this is the max index stored.
                else
                    text = "2D n=" + std::to_string(n) + " to " + std::to_string(n + spacing - 1);
            }

            outfile << text << ": " << counterArrays.events2D[i] << '\n';
        }
    }
}
