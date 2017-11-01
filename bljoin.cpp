/*
 * Unix join command, but for biological sequences.
 *
 * This program is used to join together records with matching names
 * in different files (and actually in the same file too). Useful for
 * concatenating gene sequences.
 *
 * Assumes just one of each unique identifier in every file - otherwise
 * there are problems.
 *
 * Also assumes that every file is an alignment, so that it only checks
 * the size of the first record in each file. However, the ends of the
 * sequences are padded to even up the length.
 *
 */

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <seqan/seq_io.h>

#include <tclap/CmdLine.h>

#include <SeqFileInWrapper.h>

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::ifstream;
using std::istream;
using std::map;
using std::pair;
using std::string;
using std::vector;

using namespace seqan;
using namespace bltools;

int main(int argc, char * argv[]) {
  
  TCLAP::CmdLine cmd("Equivalent of `join' for sequence files", ' ', "0.0");
  TCLAP::SwitchArg no_pad_arg("n", "no-pad",
                              "Do not attempt to pad sequences to the same length",
                              cmd);
  TCLAP::SwitchArg ignore_case_arg("i", "ignore-case",
                                   "Ignore case when matching", cmd);
  TCLAP::ValueArg<int> field_arg("f", "field",
                                 "Field (1-based) of ID to join on (after splitting); will ignore sequences for which the field is empty; default is to use whole ID",
                                 false, 0, "int", cmd);
  TCLAP::ValueArg<string> delim_arg("d", "delim",
                                 "Field separator", false, " ", "int", cmd);
  TCLAP::UnlabeledMultiArg<string> files("FILE(s)", "filenames", false,
                                         "file name(s)", cmd, false);
  cmd.parse(argc, argv);
  bool ignore_case = ignore_case_arg.getValue();
  bool no_pad = no_pad_arg.getValue();
  int field = field_arg.getValue();
  string delim = delim_arg.getValue();
  vector<string> infiles = files.getValue();
  if(infiles.size() == 0) infiles.push_back("-");

  CharString id;
  CharString seq_;              // CharString more flexible than Dna5String
  string seq;
  SeqFileInWrapper seq_handle;
  map<string, string> seqs;
  unsigned long total_bases = 0;
  unsigned long seq_size = 0;

  for(string& infile: infiles) {
    
    seq_size = 0;

    try {
        seq_handle.open(infile);
    } catch(Exception const &e) {
      cerr << "Could not open " << infile << endl;
      seq_handle.close();
      return 1;
    }
    
    while(!seq_handle.atEnd()) {

      try {

        readRecord(id, seq_, seq_handle.sqh);
        //nrecs_read++; 

      } catch (Exception const &e) {

        cerr << "Error: " << e.what() << endl;
        seq_handle.close();
        return 1;

      } // End try-catch for record reading.
      
      // Check the size of the first sequence in the file
      if(seq_size == 0) {
        seq_size = length(seq_);
      }
      
      if(seq_size != length(seq_)) {
        cerr << "Warning " << id << " is not the same size as other seqs"
             << " in the same file" << endl;
      }
      
      seq = string(toCString(seq_));
      
      // Simple method: just use the whole sequence ID to join
      string join_id = string(toCString(id));

      // Fancier: use parts of the ID (Needs conversion to C++)
      /*
      if(ignore_case) {
        // Convert to upper case
        // TODO: This needs to be converted to C++ - toupper won't work
        join_id = toupper(join_id)
      }
      */
      /*
      if(field > 0) {
        // TODO: this also needs to be converted to C++
        vector<string> splitted = join_id.split(delim);
        if(splitted.size() < (field - 1)) {
          // If the field is not found in this record, skip it
          continue;
        }
        join_id = splitted[field];
      }
      */

      // TODO: Test if id is in map, and add if not
      // Also add enough sequence to fill in missed sequences
      map<string, string>::iterator it = seqs.find(join_id);
      if(it != seqs.end()) {
        // Found: do nothing
      } else {
        if(total_bases > 0 && !no_pad) {
          seqs[join_id].reserve(total_bases + seq_size * 2);
          for(unsigned long j = 0; j < total_bases; j++) {
            seqs[join_id] += "-";
          }
        } else {
          seqs[join_id] = "";
        }
      } // End test for existence of ID
      
      // Add the current sequence
      seqs[join_id] += seq;
      
    } // End single file reading loop
    
    total_bases += seq_size;
    
    if(!no_pad) {
      for(map<string, string>::iterator it = seqs.begin(); it != seqs.end(); it++) {
        if(it->second.length() < total_bases) {
          for(unsigned long j = it->second.length(); j < total_bases; j++) {
            it->second += "-";
          }
        }
      }
    }

    if(!seq_handle.close()) {
        cerr << "Problem closing " << infile << endl;
        return 1;
    }

  } // End loop over files
  
  // Write the output in fasta format
  for(pair<string, string> item: seqs) {
    cout << ">" << item.first << endl;
    cout << item.second << endl;
  }

  return 0;
}
