/*

Terms:
  Shared complex - a complex that spreads over multiple ranks
  Complex portion - part of a complex belonging to a rank
  Holder - a rank that holds a portion of shared complex
  Complex owner rank, or complex owner, or owner - rank that is responsible for:
   - collecting requests from other ranks holding portions of the complex,
   - making decisions who associates
   - updating the complex
   - propagating decisions to holders of portions of the complex.

Communication among ranks for synchronization
  Each rank communicates with:
   - owners of its portion of complexes
   - holders of its complexes.

Communication parameters
  - complexList is a vector of Complex;
    Complex.ownerRank:
    Each complex (or portion of complex) keeps the owner rank number
  - myComplexes is an unordered map of MPIComplexInfo;
    MPIComplexInfo.minRank and MPIComplexInfo.maxRank:
    each complex owner keeps holder rank numbers

Example:
  Legend:
    - "x" - complex owner
    - { } - list or rank numbers
    - Owners - owners of portions of complexes from current rank
    - Holders - range of holders of portions of complexes that current rank owns
  Current rank: Rank 0    Rank 1      Rank 2      Rank 3 ...
  Complex 1:        ----------x-----------
  Complex 2:                                ---------x---
  Complex 3:    -----x---------
  Owners:      {1,0}      {1,0}       {1,3}       {3}
  Holders:     {0-1}      {0-2}                   {2-3}

All ranks do the same:

  (as a rank that tries to associate to the shared complex):
  Each time a molecule tries to associate with a complex,
  instead of associating (as in the case of the serial run),
  request reaction approval from owner,
  sending to the owner rotation possibilities
  and the number of molecules that tried to associate.

  (as owner):
  Collect reaction requests and possible rotations
  Overlap rotation possibilities and make random decision
  proportional to the number of molecules on each rank.
  Propagate decision to holders.

  (as a rank that tried to associate to the shared complex):
  Receive decision
  React if approved
  Rotate and translate all complex molecules

Usage:
  molecule of a complex enters ghosted zone:
  if last molecule of a complex exited  =>
  complex potentially not shared any more,
  else, update complex owner (or creates one); TODO: optimization:
  keep neighbor ranks sharing a complex informed,
  so that they don't need to communicate with the owner about this.

*/

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

// Action types:
#define COMPLEX_REQUEST_ASSOCIATE 1  // molecule trying to associate
#define COMPLEX_ASSOCIATE 2          // owner rank decided to allow associate
#define COMPLEX_MOVE 3               // owner rank decides and propagates
#define COMPLEX_ADD_RANK 4           // owner rank updates complex holders list
#define COMPLEX_REMOVE_RANK 5        // owner rank updates complex holders list
#define COMPLEX_UPDATE 6  // owner rank propagates updates to complex holders

// Complex message decoding:
// 1. complex ID
// 2. Complex action type:
//    COMPLEX_REQUEST_ASSOCIATE:
//      3. id complex
//      4. id molecule trying to associate
//      5. required translation (x, y, z)
//      6. required rotation angles
//    COMPLEX_REQUEST_NOTHING
//    COMPLEX_ASSOCIATE:
//      3. Allow portion-holder to associate
//      4. Translate coordinates (x, y, z)
//      5. Angles to rotate
//    COMPLEX_MOVE:
//      3. Translate coordinates (x, y, z)
//      4. Angles to rotate
//    COMPLEX_ADD_RANK
//    COMPLEX_REMOVE_RANK
//    COMPLEX_UPDATE:
//      3. fromRank
//      4. toRank
//      5. ownerRank
//      6. mass
//      7. center of mass

// TODO: Complex rankMin, rankMax?
void introduce_complex_portion() {}

// TODO:
void delete_complex_portion() {}

// TODO:
void share_complex() {}

// Initializes MPIComplexInfo struct for one owned complex:
void add_my_complex(MpiContext& mpiContext, int id, int minRank, int maxRank) {
  MPIComplexInfo* mTemp = new MPIComplexInfo();
  mTemp->minRank = minRank;
  mTemp->maxRank = maxRank;
  mpiContext.myComplexes[id] = mTemp;
}

// Deleting MPIComplexInfo struct for a given complex ID:
void delete_my_complex(MpiContext& mpiContext, int id) {
  // Deleting memory reserved for struct MPIComplexInfo:
  delete mpiContext.myComplexes[id];
  // Erasing entry:
  mpiContext.myComplexes.erase(id);
}

// Finding the MPIComplexInfo struct for a given complex ID:
MPIComplexInfo* find_my_complex(MpiContext& mpiContext, int id) {
  auto itr = mpiContext.myComplexes.find(44);  // e.g. 44
  if (itr == mpiContext.myComplexes.end())
    return nullptr;
  else
    return itr->second;
}

// Before each send and receive of complex events between ranks,
// make sure the counters are set to 0:
void reset_statistical_data(MpiContext& mpiContext) {
  for (auto& it : mpiContext.myComplexes) {
    it.second->iRankAssociationTries.clear();
    it.second->nRankAssociationTries.clear();
    it.second->nAssociationTries = 0;
  }
}

// Communication between ranks over complex related events:
void send_and_receive(MpiContext& mpiContext) {
  // Reserve memory for temporary arrays for receiving data from ranks:
  std::vector<unsigned char*> arrayRank;
  arrayRank.reserve(mpiContext.nprocs);
  for (int iRank = 0; iRank < mpiContext.nprocs; iRank++) {
    unsigned char* tempArrayRank;
    if (iRank != mpiContext.rank) {
      tempArrayRank = (unsigned char*)malloc(
          10000);  // TODO: determine size (receive size first)
    } else {
      tempArrayRank = nullptr;  // no need for a rank to send himself
    }
    arrayRank.push_back(tempArrayRank);
  }

  // Start sending and receiving complex requests in the background
  for (int iRank = 0; iRank < mpiContext.nprocs; iRank++) {
    if (iRank != mpiContext.rank) {  // no need for a rank to send himself
      // Sending data to complex owners
      // TODO: Optimize: Send only to complex owners and listen only to complex
      // owners+holders of portion of my complexes
      MPI_Isend(mpiContext.toRank[iRank].data(),
                mpiContext.toRank[iRank].size(), MPI_CHAR, iRank, 0,
                MPI_COMM_WORLD, &mpiContext.requestsSend[iRank]);

      // Receiving requests as complex owner
      // and receiving answers from complex owners
      MPI_Irecv(arrayRank[iRank], 10000, MPI_CHAR, iRank, 0, MPI_COMM_WORLD,
                &mpiContext.requestsRecv[iRank]);
    }
  }

  // Waiting for receiving to finish, updating fromRank,
  // and freeing up reserved memory:
  for (int iRank = 0; iRank < mpiContext.nprocs; iRank++)
    if (iRank != mpiContext.rank) {
      int length;  // Length of received message
      MPI_Status status;
      MPI_Wait(&mpiContext.requestsRecv[iRank],
               &status);  // waiting for receiving to finish
      MPI_Get_count(&status, MPI_CHAR, &length);
      // If anything received into arrayRank, put it into vector
      // fromRank[iRank]:
      if (length) {
        // if(length) std::cout << "rank " << mpiContext.rank << " received " <<
        // length << " bytes from " << iRank << std::endl;
        mpiContext.fromRank[iRank].clear();
        for (int nChr = 0; nChr < length; nChr++)
          mpiContext.fromRank[iRank].push_back(arrayRank[iRank][nChr]);
      }
      // Free reserved memory:
      free(arrayRank[iRank]);
    }

  // Waiting for sending to finish and clearing toRank:
  for (int iRank = 0; iRank < mpiContext.nprocs; iRank++) {
    if (iRank != mpiContext.rank) {  // no rank sending itself
      // Waiting for sending to finish:
      MPI_Wait(&mpiContext.requestsSend[iRank], MPI_STATUS_IGNORE);
      // Clearing the vector of bytes prepared from this rank for iRank,
      // as it is already sent:
      mpiContext.toRank[iRank].clear();
    }
  }
}

// Initialize sending and receiving vectors before first use,
// and clear after if needed:
void init_send_recv_vectors(MpiContext& mpiContext) {
  // Reset toRank and fromRank vectors:
  mpiContext.toRank.clear();
  mpiContext.fromRank.clear();
  for (int i = 0; i < mpiContext.nprocs; i++) {
    std::vector<unsigned char> temp;
    mpiContext.toRank.push_back(temp);
    mpiContext.fromRank.push_back(temp);
  }
}

// Send COMPLEX_ADD_RANK request
// to the owner of a complex with given id:
void add_rank_to_complex(MpiContext& mpiContext, int iRank, int iOwner,
                         int idComplex) {
  if (mpiContext.rank == iRank) {
    PUSH_TO(idComplex, iOwner);
    PUSH_TO(COMPLEX_ADD_RANK, iOwner);
  }
}

// Send COMPLEX_REMOVE_RANK request
// to the owner of a complex with given id:
void remove_rank_from_complex(MpiContext& mpiContext, int iRank, int iOwner,
                              int idComplex) {
  if (mpiContext.rank == iRank) {
    PUSH_TO(idComplex, iOwner);
    PUSH_TO(COMPLEX_REMOVE_RANK, iOwner);
  }
}

// Collect requests from other ranks for complexes this rank owns
// and act accordingly:
int collect_and_act(MpiContext& mpiContext) {
  int toSend = 0;
  reset_statistical_data(mpiContext);
  // Parsing received data:
  // cout << "Rank " << rank << " parsing received data" << endl;
  for (int iRank = 0; iRank < mpiContext.nprocs; iRank++) {
    if (iRank != mpiContext.rank) {
      // Prepare received bytes for deserializing using POP macro:
      unsigned char* arrayRank = mpiContext.fromRank[iRank].data();
      // Prepare deserializing (reading) position:
      int nArrayRank = 0;
      // Initialize total number of received bytes from rank iRank:
      int nBytes = mpiContext.fromRank[iRank].size();
      // Parse received bytes from rank iRank:
      while (nArrayRank < nBytes) {
        // if(nBytes) std::cout << "rank " << rank <<
        //         " looping requests from iRank = " << iRank
        //         << ", nBytes = " << nBytes << std::endl;
        int complexId,    // which complex does the request apply to
            requestType,  // type of the request
            idMolecule,   // ID of a molecule that tried association or so
            nAssociationTriesFromRank;  // number of molecules that tried
                                        // to associate from iRank
        MPIComplexInfo* m;
        // Call POP macro to extract bytes from arrayRank
        // needed for an integer; store these in complexID variable:
        POP(complexId);
        // Extract requestType from the following bytes of arrayRank:
        POP(requestType);
        // Perform action based on the type of request to the owner:
        switch (requestType) {
          case COMPLEX_REQUEST_ASSOCIATE:
            // Collect all requests for complex complexId
            // before making decision who associates, if anyone:
            POP(idMolecule);  // TODO: decide what you want to do
            POP(nAssociationTriesFromRank);
            std::cout << "Rank " << mpiContext.rank
                      << " received request COMPLEX_REQUEST_ASSOCIATE by rank "
                      << iRank << std::endl;
            // Store received request into MPIComplexInfo struct:
            m = mpiContext.myComplexes[complexId];
            m->nAssociationTries += nAssociationTriesFromRank;
            m->iRankAssociationTries.push_back(iRank);
            m->nRankAssociationTries.push_back(nAssociationTriesFromRank);
            break;
          case COMPLEX_ASSOCIATE:
            std::cout << "Rank " << mpiContext.rank
                      << " received requirement COMPLEX_ASSOCIATE by rank "
                      << iRank << " and associating." << std::endl;
            // TODO
            break;
          case COMPLEX_MOVE:
            std::cout << "Rank " << mpiContext.rank
                      << " received requirement COMPLEX_MOVE by rank " << iRank
                      << " and moving." << std::endl;
            // TODO
            break;
          case COMPLEX_ADD_RANK:
            // Update minRank or maxRank:
            std::cout << "Rank " << mpiContext.rank
                      << " received requirement COMPLEX_ADD_RANK by rank "
                      << iRank << " and adding a rank." << std::endl;
            m = mpiContext.myComplexes[complexId];
            if (m->minRank - 1 == iRank)
              m->minRank = iRank;
            else if (m->maxRank + 1 == iRank)
              m->maxRank = iRank;
            break;
          case COMPLEX_REMOVE_RANK:
            // Update minRank or maxRank:
            std::cout << "Rank " << mpiContext.rank
                      << " received requirement COMPLEX_REMOVE_RANK by rank "
                      << iRank << " and updating." << std::endl;
            m = mpiContext.myComplexes[complexId];
            if (m->minRank == iRank)
              m->minRank++;
            else if (m->maxRank == iRank)
              m->maxRank--;
            break;
          case COMPLEX_UPDATE:
            std::cout << "Rank " << mpiContext.rank
                      << " received requirement COMPLEX_UPDATE by rank "
                      << iRank << " and updating." << std::endl;
            // TODO
            break;
        }
      }
      // As everything received from rank iRank is parsed,
      // clear vector of bytes from rank iRank to this rank:
      mpiContext.fromRank[iRank].clear();
    }
  }

  // Parse and decide what to do with owned complexes:
  for (auto& complex :
       mpiContext.myComplexes) {  // for all complexes this rank owns
    if (complex.second
            ->nAssociationTries) {  // if any rank required association
      int id = complex.first;       // complex ID
      cout << endl
           << "Rank " << mpiContext.rank << " decides on complex " << id
           << ", minRank = " << complex.second->minRank
           << ", maxRank = " << complex.second->maxRank << endl;
      // Number of molecule that will be accepted to associate:
      int nAccept = rand() % complex.second->nAssociationTries;
      // cout << "Accepting molecule number " << nAccept << endl;
      //  Find what rank holds nAccept-th request if they are put in a row
      //  e.g. rank 1 has 3 requests, rank 5 has 4. Row: 1 1 1 5 5 5 5:
      int sum = 0;          // sum of requests up to currently observed rank
      int iAssocTries = 0;  // currently observed rank
      while (sum + complex.second->nRankAssociationTries[iAssocTries] < nAccept)
        sum += complex.second->nRankAssociationTries[iAssocTries++];
      int acceptRank = complex.second->iRankAssociationTries[iAssocTries];
      cout << "Rank " << mpiContext.rank
           << " accepts COMPLEX_REQUEST_ASSOCIATE from " << acceptRank << endl;

      PUSH_TO(id, acceptRank);  // e.g. 44
      PUSH_TO(COMPLEX_ASSOCIATE, acceptRank);
      std::cout << "Rank " << mpiContext.rank << " pushed to rank "
                << acceptRank << ", requiring COMPLEX_ASSOCIATE\n";

      for (int iRank = complex.second->minRank;
           iRank <= complex.second->maxRank; iRank++) {
        if (iRank != mpiContext.rank) {
          PUSH_TO(id, iRank);  // e.g. 44
          PUSH_TO(COMPLEX_MOVE, iRank);
          std::cout << "Rank " << mpiContext.rank << " pushed to rank " << iRank
                    << ", requiring COMPLEX_MOVE\n";
        } else {
          cout << "Rank " << mpiContext.rank << " moving..." << endl;
        }
      }
      toSend = 1;
    }
  }
  return toSend;
}

// Check whether a rank required further communication to another rank or ranks
int check_for_more(MpiContext& mpiContext, int& toSend) {
  int toReceive = 0;  // At least one rank required further communication
  // Check whether any rank requested further complex-related communication by
  // all ranks
  if (toSend)
    cout << "Rank " << mpiContext.rank << " requires further communication."
         << endl
         << endl;
  MPI_Allreduce(&toSend, &toReceive, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  // if(toReceive) std::cout << "Rank " << rank << " needs further
  // communication. Received communication requests " << toReceive << std::endl;
  return toReceive;
}

// Free resources for reserved by complex owners for MPIComplexInfo structs:
void free_myComplexes(MpiContext& mpiContext) {
  // Traversing an unordered map
  for (auto x : mpiContext.myComplexes) delete x.second;
}

void example_complex_scenario(MpiContext& mpiContext,
                              vector<Complex>& complexList) {
  // Scenario:
  // Rank 1-4 hold portions of complex e.g. 44 and of 66, both owned by rank 3

  // Initialize 2 complexes:
  if ((mpiContext.rank >= 1) && (mpiContext.rank <= 4)) {
    if (complexList.size() < 2) {
      Complex c;
      c.id = 44;  // e.g. 44
      c.ownerRank = 3;
      complexList.push_back(c);
      complexList.push_back(c);
    }
    Complex& c = complexList[0];
    c.id = 44;
    c.ownerRank = 3;
    // complexList.push_back(c);
    // c.id = 66;
    // complexList.push_back(c);
    complexList[1].id = 66;
  }

  if (mpiContext.rank == 3) {
    add_my_complex(mpiContext, 44, 1, 4);  // e.g. 44
    add_my_complex(mpiContext, 66, 1, 4);
  }
  // delete_my_complex(mpiContext, 66);
  // MPIComplexInfo *m = find_my_complex(mpiContext, 44); //e.g. 44

  // Ranks 1 and 2 sending to complex owner over complex e.g. 44:
  if ((mpiContext.rank == 1) || (mpiContext.rank == 2)) {
    int toRank = complexList[0].ownerRank;
    PUSH_TO(complexList[0].id, toRank);  // e.g. 44
    int requestType = COMPLEX_REQUEST_ASSOCIATE;
    PUSH_TO(requestType, toRank);
    PUSH_TO(77 + 1000 * mpiContext.rank, toRank);  // molecule ID
    PUSH_TO(3, toRank);  // how many molecules tried associating
    std::cout << "Rank " << mpiContext.rank
              << " pushed COMPLEX_REQUEST_ASSOCIATE and "
              << 77 + 1000 * mpiContext.rank << " to rank "
              << complexList[0].ownerRank << " for complex "
              << complexList[0].id << "\n";
  }

  // TODO: synchronize with owner actions:
  // add_rank_to_complex(mpiContext, 5, 3, complexList[0].id); // add rank 5 to
  // complex owner 3 as holder of complexID remove_rank_from_complex(mpiContext,
  // 5, 3, complexList[0].id); // add rank 5 to complex owner 3 as holder of
  // complexID
}

int manage_complexes(int argc, char** argv) {
  // int main(int argc, char **argv){
  cout << endl;
  // Common MPI, MpiContext, random, and complexList initialization:
  MpiContext mpiContext(MPI::COMM_WORLD.Get_size(), 100000);
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiContext.rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpiContext.nprocs);
  int rank = mpiContext.rank;
  srand(time(NULL));
  // Mock list of complexes:
  vector<Complex> complexList;

  // Initialize sending and receiving vectors before first use:
  init_send_recv_vectors(mpiContext);

  // TODO: example scenario is for demonstrating purposes:
  example_complex_scenario(mpiContext, complexList);

  // Communicate with complex owners as long as any rank requires it:
  int toReceive;  // At least one rank required further communication
  do {
    // Send and receive everything ranks have prepared for each other:
    send_and_receive(mpiContext);
    // Collect all requests to each complex by its owner and
    // make and propagate decisions:
    int toSend = collect_and_act(mpiContext);
    // Check whether another round of communication between ranks is needed:
    toReceive = check_for_more(mpiContext, toSend);
  } while (toReceive);  // continue communication while any rank requires it
  // Free resources for reserved by complex owners for MPIComplexInfo structs:
  free_myComplexes(mpiContext);

  // Common MPI and main finalization:
  MPI_Finalize();

  return 0;
}
