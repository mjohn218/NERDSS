#include "classes/class_Parameters.hpp"
#include <iostream>
#include <sstream>
#include <string>

void parse_command(int argc, char* argv[], Parameters& params, std::string& paramFileName, std::string& restartFileName, std::string& addFileName, unsigned int& seed)
{
    std::cout << "Command: " << std::string(argv[0]) << std::flush;
    for (int flagItr { 1 }; flagItr < argc; ++flagItr) {
        std::string flag { argv[flagItr] };
        std::cout << ' ' << flag << std::flush;
        if (flag == "-f") {
            paramFileName = argv[flagItr + 1];
            std::cout << ' ' << std::string(argv[flagItr + 1]) << std::flush;
            ++flagItr;
        } else if (flag == "-s" || flag == "--seed") {
            std::stringstream iss(argv[flagItr + 1]);
            unsigned tmpseed;
            if (iss >> tmpseed) {
                seed = tmpseed;
                std::cout << ' ' << seed << std::flush;
            } else {
                std::cerr << "Error reading seed, exiting.\n";
                exit(1);
            }
            ++flagItr;
            std::cout << '\n';
        } else if (flag == "--debug-force-dissoc") {
            params.debugParams.forceDissoc = true;
        } else if (flag == "--debug-force-assoc") {
            params.debugParams.forceAssoc = true;
        } else if (flag == "--print-system-info") {
            params.debugParams.printSystemInfo = true;
        } else if (flag == "-r" || flag == "--restart") {
	  if(params.rank<0){//for serial jobs
	    restartFileName = std::string(argv[flagItr + 1]);
            std::cout << ' ' << std::string(argv[flagItr + 1]) << std::flush;
	  }else{//for parallel jobs
	    restartFileName = "restart.dat";
	    char rankChar[10];
	    sprintf(rankChar, "%d",params.rank);
	    restartFileName+=rankChar;
	    std::cout << ' ' << restartFileName << std::flush;
	  }
	  params.fromRestart = true;
	  ++flagItr;
        } else if (flag == "-a" || flag == "--add") {
            addFileName = std::string(argv[flagItr + 1]);
            std::cout << ' ' << std::string(argv[flagItr + 1]) << std::flush;
            ++flagItr;
        } else if (flag == "-v") {
            params.debugParams.verbosity = 1;
        } else if (flag == "-vv") {
            params.debugParams.verbosity = 2;
        } else {
            std::cout << " ignored " << std::endl;
        }
    }
}
