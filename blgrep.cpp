/*
 * Seqan-based program for grepping sequence files. Right now it uses
 * a regex approach to all the matching, but I would like to implement
 * a biological matching approach too.
 *
 * TODO:
 *    - Output format based on input format
 *    - Handle input formats with quality scores
 *
 */

#include <iostream>
#include <regex>
#include <string>

#include <seqan/find.h>
#include <seqan/seq_io.h>

#include <tclap/CmdLine.h>

using std::cout;
using std::cerr;
using std::endl;
using std::regex;
using std::regex_match;
using std::string;

using namespace seqan;

int main(int argc, char * argv[]) {
  
  // Regex flags
  std::regex_constants::syntax_option_type regex_flags =
    regex::extended | regex::optimize;
  std::regex_constants::match_flag_type regex_match_flags =
    std::regex_constants::match_any | std::regex_constants::match_not_null;
  
  /*
   * Command line arguments.
   */
  TCLAP::CmdLine cmd("Program to search fasta file", ' ', "0.0");
  TCLAP::SwitchArg seq_regex_arg("S", "sequence-regex",
				 "use regex for sequences instead of name; sets -i",
				 cmd);
  TCLAP::SwitchArg invert_regex_arg("v", "invert-match",
				    "Invert matching, like grep -v", cmd);
  TCLAP::SwitchArg ignore_case_arg("i", "ignore-case",
				   "Ignore case in pattern and input", cmd);
  TCLAP::SwitchArg case_sensitive_arg("I", "case-sensitive",
				      "Do not ignore case in pattern and input; only for -S",
				      cmd);
  TCLAP::UnlabeledValueArg<string> regex_string_arg("pattern", "regex pattern",
						    true, "",
						    "regex", cmd, false);
  TCLAP::UnlabeledValueArg<string> infile_name("infile", "input file", true,
					       "", "file name", cmd, false);
  cmd.parse(argc, argv);
  const char* infile = infile_name.getValue().c_str();
  bool seq_regex = seq_regex_arg.getValue();
  bool inverted = invert_regex_arg.getValue();
  bool case_sensitive = case_sensitive_arg.getValue();
  if(ignore_case_arg.getValue() || (seq_regex && !case_sensitive))
    regex_flags |= regex::icase;
  regex regex_pattern(regex_string_arg.getValue(), regex_flags);
  
  CharString id;
  CharString seq;              // CharString more flexible than Dna5String
  SeqFileIn seq_handle(infile);

  if(!open(seq_handle, infile)) {
    cerr << "Could not open " << infile << endl;
    return 1;
  }
  
  bool matched = false;
  int nmatched = 0;
  while(!atEnd(seq_handle)) {

    try {

      readRecord(id, seq, seq_handle);
      
      // Regex
      if(seq_regex) {
	matched = regex_search(toCString(seq), regex_pattern,
			       regex_match_flags);
      } else {
	matched = regex_search(toCString(id), regex_pattern,
			       regex_match_flags);
      } // End regex
      
      // Write out if matched
      if((matched && !inverted) || (!matched && inverted)) {
	nmatched++;
	cout << ">" << id << endl << seq << endl;
      }

    } catch (Exception const &e) {

      cout << "Error: " << e.what() << endl;
      return 1;

    } // End try-catch for record reading.

  } // End file reading loop
  
  if(nmatched) {
    return 0;
  } else {
    return 1;
  }

}
