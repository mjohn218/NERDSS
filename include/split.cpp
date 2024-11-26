#pragma once

#include <math.h>
#ifdef mpi_
#include <mpi.h>
#endif
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "classes/mpi_functions.hpp"

#define ACTION_INCREASE_SPECIES 1
#define ACTION2 2
#define ACTION3 3
#define ACTION4 4

#define REQUIRE_ROTATE 1
#define REQUIRE2 2

#define min2(a, b) ((a) < (b) ? (a) : (b))

struct Membrane;
struct SimulVolume;
struct Molecule;
struct Complex;

// MPIComplexInfo holds information about complex needed for the complex owner:
struct MPIComplexInfo {
 public:
  // How big is the complex in terms of number of ranks:
  int minRank, maxRank;
  // Vector holding rank numbers that requested association to the complex:
  std::vector<int> iRankAssociationTries;
  // Vector holding number of association requests by these ranks:
  std::vector<int> nRankAssociationTries;
  // Total number of molecules that requested association to the rank:
  int nAssociationTries = 0;
};

typedef struct structMpiContext {  // holds MPI related data
  structMpiContext(int nRanks, int neigborRankBufferSize) {
    xOffset = 0;
    startCell = 0;
    endCell = 0;
    startGhosted = 0;
    endGhosted = 0;
    nprocs = nRanks;

    #ifdef mpi_
    // Reserve memory for sending/receiving requests regarding complexes:
    requestsSend = (MPI_Request*)malloc(sizeof(MPI_Request) * nRanks);
    requestsRecv = (MPI_Request*)malloc(sizeof(MPI_Request) * nRanks);
    #endif

    // Reserve memory for communication with neighboring ranks:
    sendBufferSize = neigborRankBufferSize;
    recvBufferSize = neigborRankBufferSize;
    MPIArrayFromLeft = (unsigned char*)malloc(recvBufferSize);
    MPIArrayToLeft = (unsigned char*)malloc(sendBufferSize);
    MPIArrayFromRight = (unsigned char*)malloc(recvBufferSize);
    MPIArrayToRight = (unsigned char*)malloc(sendBufferSize);
  }

  void init_x_domain_and_offset(int xTotal, int& tempRank) {
    // Distribute as similar number of cells among x axis as possible per rank,
    // std::cout << "Cells per x axis: " << simulVolume.numSubCells.x <<
    // std::endl; Splitting cells onto ranks is done in greedy-manner, but
    // without having two ranks that process number of cells that differ by more
    // than 1 e.g. 12 cells onto 5 ranks is divided as: 3, 3, 2, 2, 2.
    int divCells = xTotal / nprocs;
    int modCells = xTotal % nprocs;
    // Determine the domain of cells that tempRank is responsible for in gready
    // manner, e.g. divide 12 cells onto 5 ranks as: 3, 3, 2, 2, 2:
    startCell = divCells * tempRank + min2(tempRank, modCells);
    endCell = divCells * (tempRank + 1) + min2(tempRank + 1, modCells) - 1;
    // Include ghosted zones:
    startGhosted = startCell;
    endGhosted = endCell;
    if (startGhosted > 0) startGhosted--;
    if (endGhosted < xTotal - 1) endGhosted++;
    // xOffset represents x bin coordinate of most left cell belonging to
    // current rank; It is used for calculating currBin in update_memberMolLists
    // that requires substracting this offset.
    xOffset = startGhosted;
    // printf("== Rank %2d cell-domain:   startCell: %3d   endCell: %3d, diff+1:
    // %3d ==\n", tempRank, startCell, endCell, endCell - startCell + 1);
    //  determine the boundary of the rank in x axis
    int leftNonGhosted = 0;
    if (tempRank > 0) leftNonGhosted = xOffset + 1;
    xLeft = 1.0 * leftNonGhosted / xTotal;
    int rightNonGhosted = xOffset + (endCell - startCell + 1);
    if (tempRank > 0) rightNonGhosted = xOffset + 1 + (endCell - startCell + 1);
    xRight = 1.0 * rightNonGhosted / xTotal;
    if (tempRank == nprocs - 1) xRight = 1.0;
  }

  inline void check_buffer_size(int serializedSize) {
    if (serializedSize > sendBufferSize * 0.8) {
      std::cout << " --- WARNING ---" << std::endl;
      std::cout << "Serialized buffer size = " << serializedSize << std::endl;
      std::cout << "Sending buffer size = " << sendBufferSize << std::endl;
    }
    if (serializedSize > sendBufferSize) {
      std::cout << "Serialized buffer size = " << serializedSize << std::endl;
      std::cout << "Sending buffer size = " << sendBufferSize << std::endl;
      std::cout << "Please increase NEIGHBOR_BUFFER_SIZE" << std::endl;
      std::cout << "Exiting..." << std::endl;
      exit(1);
    }
  }
  void increase_recv_buffers() {
    if ((nMPIArrayFromLeft > 0.8 * recvBufferSize) ||
        (nMPIArrayFromRight > 0.8 * recvBufferSize)) {
      free(MPIArrayFromLeft);
      free(MPIArrayFromRight);
      recvBufferSize *= 1.2;
      MPIArrayFromLeft = (unsigned char*)malloc(recvBufferSize);
      MPIArrayFromRight = (unsigned char*)malloc(recvBufferSize);
    }
  }
  void increase_send_buffers() {
    if ((nMPIArrayToLeft > 0.8 * sendBufferSize) ||
        (nMPIArrayToRight > 0.8 * sendBufferSize)) {
      free(MPIArrayToLeft);
      free(MPIArrayToRight);
      sendBufferSize *= 1.2;
      MPIArrayToLeft = (unsigned char*)malloc(sendBufferSize);
      MPIArrayToRight = (unsigned char*)malloc(sendBufferSize);
    }
  }
  ~structMpiContext() {
    #ifdef mpi_
    free(requestsSend);
    free(requestsRecv);
    #endif
    free(MPIArrayFromLeft);
    free(MPIArrayToLeft);
    free(MPIArrayFromRight);
    free(MPIArrayToRight);
  }
  void print_spaces() {
    if (rank) {
      printf("%*c", rank, ' ');
      printf("%*c", rank, ' ');
    }
  }

  void print_xBins() {
    printf(
        "-----xBins of rank %d: startCell = %d, endCell = %d, startGhosted = "
        "%d, endGhosted = %d, xOffset = %d\n",
        rank, startCell, endCell, startGhosted, endGhosted, xOffset);
  }

  /*
  Function serialize serializes structMpiSynchronizationData struct
  into arrayRank of bytes.
  */
  void serialize(unsigned char* arrayRank, int& nArrayRank) {
    // Determines the domain of cells that tempRank is responsible for in gready
    // manner, e.g. divide 12 cells onto 5 ranks as: 3, 3, 2, 2, 2:
    PUSH(xOffset);
    PUSH(startCell);
    PUSH(endCell);
    PUSH(startGhosted);
    PUSH(endGhosted);
    PUSH(xLeft);
    PUSH(xRight);
  }
  /*
  Function deserialize deserializes structMpiSynchronizationData struct
  from arrayRank of bytes.
  */
  void deserialize(unsigned char* arrayRank, int& nArrayRank) {
    POP(xOffset);
    POP(startCell);
    POP(endCell);
    POP(startGhosted);
    POP(endGhosted);
    POP(xLeft);
    POP(xRight);
  }

  #ifdef mpi_
  // For non-blocking sending of actions:
  MPI_Request requestSendToLeft, requestSendToRight, requestRecvFromLeft,
      requestRecvFromRight;
  MPI_Status statusRecvFromLeft,
      statusRecvFromRight;  // for counting number of doubles received
  #endif

  int sendBufferSize, recvBufferSize;
  unsigned char* MPIArrayToRight;
  unsigned char* MPIArrayFromRight;
  unsigned char* MPIArrayToLeft;
  unsigned char* MPIArrayFromLeft;
  int nMPIArrayToRight;
  int nMPIArrayFromRight;
  int nMPIArrayToLeft;
  int nMPIArrayFromLeft;

  int xOffset;
  // Determine the domain of cells that tempRank is responsible for in gready
  // manner, e.g. divide 12 cells onto 5 ranks as: 3, 3, 2, 2, 2:
  int startCell, endCell;
  int startGhosted, endGhosted;
  double xLeft, xRight;

  int nprocs;

  int rank,
      size;  // MPI current rank and number of processing elements per width

  int simItr;
  bool checkUnimoleculeReactionPopulation;
  bool hasRankCommunicationForLargeComplex;

  #ifdef mpi_
  // MPI_Requests for sending requests to other ranks regarding complexes:
  MPI_Request *requestsSend, *requestsRecv;
  #endif

  // Vector of vectors will hold information one rank has to send to other
  // ranks, and receive from them:
  std::vector<std::vector<unsigned char> > toRank;
  std::vector<std::vector<unsigned char> > fromRank;

  // For a complex index, define information about which complex holder requires
  // what:
  std::unordered_map<int, MPIComplexInfo*> myComplexes;

  Membrane* membraneObject;
  SimulVolume* simulVolume;
  std::vector<Molecule>* moleculeList;
  std::vector<Complex>* complexList;
} MpiContext;
