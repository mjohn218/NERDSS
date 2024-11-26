
#include "classes/mpi_functions.hpp"

#include <array>
#include <chrono>
#include <iostream>
#include <string>
#include <vector>

#include "mpi.h"

void serialize_string(std::string to_serialize, unsigned char *arrayRank,
                      int &nArrayRank) {
  *((int *)(arrayRank + nArrayRank)) = to_serialize.length();
  nArrayRank += sizeof(int);
  for (int i = 0; i < to_serialize.length(); i++) {
    *((char *)(arrayRank + nArrayRank)) = to_serialize[i];
    nArrayRank++;
  }
}

void deserialize_string(std::string &to_deserialize, unsigned char *arrayRank,
                        int &nArrayRank) {
  int nString;
  nString = *((int *)(arrayRank + nArrayRank));
  nArrayRank += sizeof(int);
  to_deserialize = "";
  for (int i = 0; i < nString; i++) {
    to_deserialize += *((char *)(arrayRank + nArrayRank));
    nArrayRank++;
  }
}

void serialize_vector_strings(std::vector<std::string> to_serialize,
                              unsigned char *arrayRank, int &nArrayRank) {
  // serialize to_serialize matrix of abstract type T
  *((int *)(arrayRank + nArrayRank)) =
      to_serialize
          .size();  // store number or rows in first bytes needed for an int
  nArrayRank +=
      sizeof(int);  // set current number of bytes read to the size of int
  for (auto &it : to_serialize)  // for all rows
    serialize_string(it, arrayRank,
                     nArrayRank);  // read abstract vector T
                                   // and increase number of bytes read
                                   // Note: T has to implement serialize method
}

void deserialize_vector_strings(std::vector<std::string> &to_deserialize,
                                unsigned char *arrayRank, int &nArrayRank) {
  // Note: if -std=c++14 or -std=gnu++14 is used, typename T would not have to
  // be passed, but auto could be used
  to_deserialize.clear();
  int vector_size =
      *((int *)(arrayRank +
                nArrayRank));  // extract number of rows
                               // from first bytes of array needed for one int
  nArrayRank +=
      sizeof(int);  // set current number of bytes read to the size of int
  for (int i = 0; i < vector_size; i++) {  // for all rows
    std::string temp;                      // Create local vector of type T
    deserialize_string(
        temp, arrayRank,
        nArrayRank);  // extract one vector
                      // and increase number of bytes read
                      // Note: T has to implement deserialize method
    to_deserialize.push_back(
        temp);  // push resulting vector to the back of the matrix
  }
}
