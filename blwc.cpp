/*
 * Unix wc command, but for biological sequences.
 *
 */

#include <iostream>
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
  
  TCLAP::CmdLine cmd("Equivalent of `wc' for sequence files", ' ', "0.0");
  TCLAP::SwitchArg rec_count_arg("m", "length",
                                 "Give the length of each record", cmd);
  TCLAP::SwitchArg gc_arg("g", "gc",
                          "Give the GC proportion (of file or of each record with -m)",
                          cmd);
  TCLAP::SwitchArg include_gap_arg("i", "include-gap",
                                 "Include gaps ('-') in the base count", cmd);
  TCLAP::UnlabeledMultiArg<string> files("FILE(s)", "filenames", false,
                                         "file name(s)", cmd, false);
  cmd.parse(argc, argv);
  bool include_gaps = include_gap_arg.getValue();
  bool rec_count = rec_count_arg.getValue();
  bool gc = gc_arg.getValue();
  vector<string> infiles = files.getValue();
  if(infiles.size() == 0) infiles.push_back("-");

  CharString id;
  CharString seq;              // CharString more flexible than Dna5String
  SeqFileInWrapper seq_handle;
  unsigned base_count = 0;
  unsigned gc_count = 0;

  for(string& infile: infiles) {
    base_count = 0;
    gc_count = 0;

    try {
        seq_handle.open(infile);
    } catch(Exception const &e) {
      cerr << "Could not open " << infile << endl;
      seq_handle.close();
      return 1;
    }
    
    int nrecs_read = 0;

    while(!seq_handle.atEnd()) {

      try {

        readRecord(id, seq, seq_handle.sqh);
        nrecs_read++; 

      } catch (Exception const &e) {

        cerr << "Error: " << e.what() << endl;
        seq_handle.close();
        return 1;

      } // End try-catch for record reading.
      
      if(gc || rec_count) {
        // TODO Use the size() or length() or w/e function to get the base_count
        // unless you want to avoid counting gaps
        for(char& b: seq) {
          if(include_gaps || b != '-') {
              base_count += 1;
            if(gc) {
              if(b == 'G' || b == 'C' || b == 'g' || b == 'c') {
                  gc_count += 1;
              }
            }
          }
        }
      }
      
      if(rec_count) {
        if(gc) {
          cout << infile << "\t" << id << "\t" << ((double)gc_count) / (base_count) << endl;
        } else {
          cout << infile << "\t" << id << "\t" << base_count << endl;
        }
        gc_count = 0;
        base_count = 0;
      } // End rec_count output

    } // End single file reading loop

    if(!seq_handle.close()) {
        cerr << "Problem closing " << infile << endl;
        return 1;
    }

    if(!rec_count) {
      if(gc) {
        cout << infile << "\t" << ((double)gc_count) / (base_count) << endl;
      } else {
        cout << infile << "\t" << nrecs_read << endl;
      }
    }

  } // End loop over files

  return 0;
}
