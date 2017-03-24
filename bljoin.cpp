/*
 * Unix join command, but for biological sequences.
 *
 */

#include <iostream>
#include <map>
#include <queue>
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
using std::queue;
using std::string;
using std::vector;

using namespace seqan;
using namespace bltools;

int main(int argc, char * argv[]) {
  
  TCLAP::CmdLine cmd("Equivalent of `join' for sequence files", ' ', "0.0");
  TCLAP::SwitchArg ignore_case_arg("i", "ignore-case",
                                   "Ignore case when matching", cmd);
  TCLAP::ValueArg<int> field_arg("f", "field"
                                 "Field (1-based) of ID to join on (after splitting); will ignore sequences for which the field is empty; default is to use whole ID",
                                 false, 0, "int", cmd);
  TCLAP::ValueArg<int> delim_arg("d", "delim"
                                 "Field separator", false, " ", "int", cmd);
  TCLAP::UnlabeledMultiArg<string> files("FILE(s)", "filenames", false,
                                         "file name(s)", cmd, false);
  cmd.parse(argc, argv);
  bool ignore_case = ignore_case_arg.getValue();
  int field = field_arg.getValue();
  string delim = delim_arg.getValue();
  vector<string> infiles = files.getValue();
  if(infiles.size() == 0) infiles.push_back("-");

  CharString id;
  CharString seq;              // CharString more flexible than Dna5String
  SeqFileInWrapper seq_handle;
  map<string, string> seqs;
  int total_bases = 0;
  int seq_size = 0;

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

        readRecord(id, seq, seq_handle.sqh);
        nrecs_read++; 

      } catch (Exception const &e) {

        cerr << "Error: " << e.what() << endl;
        seq_handle.close();
        return 1;

      } // End try-catch for record reading.
      
      if(seq_size == 0) {
        // TODO Is this how to get the sequence length?
        seq_size = seq.length();
      }
      
      string join_id = id;
      // Get the field to join on
      if(ignore_case) {
        // Convert to upper case
        // TODO: This needs to be converted to C++ - toupper won't work
        join_id = toupper(join_id)
      }
      if(field > 0) {
        // TODO: this also needs to be converted to C++
        vector<string> splitted = join_id.split(delim);
        if(splitted.size() < (field - 1)) {
          // If the field is not found in this record, skip it
          continue;
        }
        join_id = splitted[field];
      }

      // TODO: Test if id is in map, and add if not
      // Also add enough sequence to fill in missed sequences
      
      // TODO: check how to do this too:
      seqs[join_id].append(string(seq));
      
    } // End single file reading loop
    
    total_bases += seq_size;

    // TODO: check if a sequence was missed and fill in with gaps
    // Will also need to keep track of sequence size from previous
    // file; use first sequence

    if(!seq_handle.close()) {
        cerr << "Problem closing " << infile << endl;
        return 1;
    }

  } // End loop over files

  return 0;
}
