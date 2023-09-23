/*
Copyright (c) 2015-2016 Xiaowei Zhu, Tsinghua University

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#ifndef FILESYSTEM_HPP
#define FILESYSTEM_HPP

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "type.hpp"
#include <fcntl.h>

#ifdef _WIN32
#include "unistdWindows.hpp"
#else
#include <unistd.h>
#endif


inline bool file_exists(std::string filename) {
  struct stat st;
  return stat(filename.c_str(), &st)==0;
}

inline long file_size(std::string filename) {
  struct stat st;
  assert(stat(filename.c_str(), &st)==0);
  return st.st_size;
}

/// Represent the part of the files that will be processed by a partition
/// TODO: Make it shared
class FilePartition{
public:
  /// Constructor
  /// Open file, calcuate the data part to read
  FilePartition(std::string path,
                int current_partition, 
                int total_partition,
                std::size_t edge_uint_size);
  
  /// Set the location to the start of that partition
  void ResetReadOffset();
  /// Read data into buffer
  long ReadNext(void *buffer, size_t buffer_size);
  
private:
  int fin_;
  long start_offset_;
  long bytes_to_read_; //Excluding
  long bytes_read_;
};

// FilePartition::FilePartition(std::string path,
//                              int current_partition, 
//                              int total_partition,
//                              std::size_t edge_unit_size){
//   long total_bytes = file_size(path.c_str());

//   edge_size_t total_edges = total_bytes / edge_unit_size;
//   edge_size_t read_edges = total_edges / total_partition;

//   if(current_partition == total_partition - 1){
//     read_edges += total_edges % total_partition;
//   }

//   bytes_to_read_ = edge_unit_size * read_edges;
//   start_offset_ = edge_unit_size * (total_edges / total_partition * current_partition);
  
//   bytes_read_ = 0;
//   fin_ = open(path.c_str(), O_RDONLY);
//   assert(lseek(fin_, start_offset_, SEEK_SET) == start_offset_);
// }

// void FilePartition::ResetReadOffset(){
//   assert(lseek(fin_, start_offset_, SEEK_SET) == start_offset_);
//   bytes_read_ = 0;
// }

// long FilePartition::ReadNext(void *buffer, size_t buffer_size){
//   long curr_bytes_read;

//   if (size_t(bytes_to_read_ - bytes_read_) > buffer_size) {
//     curr_bytes_read = read(fin_, buffer, buffer_size);
//   } else {
//     curr_bytes_read = read(fin_, buffer, bytes_to_read_ - bytes_read_);
//   }
//   bytes_read_ += curr_bytes_read;

//   return curr_bytes_read;
// }
#endif