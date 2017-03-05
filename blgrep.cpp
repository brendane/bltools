/*
 * Seqan-based program for grepping sequence files. Right now it uses
 * a regex approach to all the matching, but I would like to implement
 * a biological matching approach too.
 *
 * It refuses to read genbank files with a ".gb" extension unless they are
 * piped in. I can't seem to convince SeqFileIn that genbank is the format.
 *
 * TODO:
 *    - Output format based on input format, or at least some choices
 *    - Trick from blhead for forcing SeqAn to not use the file extension
 *      to guess at input sequence format.
 *    - Checks for input file issues, and make sure file handles are
 *      closed in all error conditions.
 *    - Ability to search within a range of sequence coordinates
 *    - Ability to spit out context around matches or the number of matches
 *    - Ability to split a fasta file into records using a GFF file -
 *      or maybe that should be a different tool?
 *    - Biologically-based search instead of regex; I think if I make a
 *      bl_search function that either takes the same arguments as
 *      regex_search or runs a biological pattern matching search, I
 *      could make it work. If it is a template, I might be able to do
 *      something about the problem of treating all sequences as Dna.
 *
 */

#include <iostream>
#include <regex>
#include <string>
#include <vector>

#include <seqan/seq_io.h>
#include <seqan/translation.h>

#include <tclap/CmdLine.h>

#include <SeqFileInWrapper.h>

using std::cout;
using std::cerr;
using std::endl;
using std::regex;
using std::regex_match;
using std::string;
using std::vector;

using namespace seqan;
using namespace bltools;

int main(int argc, char * argv[]) {

  /*
   * Command line arguments.
   */
  TCLAP::CmdLine cmd("Program to regex search sequence files", ' ', "0.0");
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
  TCLAP::ValueArg<string> match_type_arg("M", "match-type",
                                         "Match type: f=fwd, r=rev, c=compl., R=revcomp; ignored for names",
                                         false, "f", "string", cmd);
  TCLAP::ValueArg<int> frame_arg("F", "frame",
                                 "Frame for translation: 0=fwd frame, 1=fwd + revcomp, 2=all 3 fwd, 3=all 6",
                                 false, 0, "string", cmd);
  TCLAP::ValueArg<string> format_arg("o", "output-format",
                                     "Output format: fasta or fastq; fasta is default",
                                     false, "fasta", "fast[aq]", cmd);
  TCLAP::UnlabeledValueArg<string> regex_string_arg("PATTERN", "regex pattern",
                                                    true, "",
                                                    "regex", cmd, false);
  TCLAP::UnlabeledMultiArg<string> infile_name("FILE(s)", "input file(s) use '-' for stdin or leave blank",
                                               false, "file name(s)", cmd, false);
  cmd.parse(argc, argv);
  vector<string> infiles = infile_name.getValue();
  if(infiles.size() == 0) infiles.push_back("-"); // stdin
  bool seq_regex = seq_regex_arg.getValue();
  bool inverted = invert_regex_arg.getValue();
  string match_type = match_type_arg.getValue();
  int frame = frame_arg.getValue();
  string format = format_arg.getValue();

  // Regex setup
  std::regex_constants::syntax_option_type regex_flags =
    regex::extended | regex::optimize;
  std::regex_constants::match_flag_type regex_match_flags =
    std::regex_constants::match_any | std::regex_constants::match_not_null;
  if(ignore_case_arg.getValue() ||
     (seq_regex && !case_sensitive_arg.getValue())) {
    regex_flags |= regex::icase;
  }
  regex regex_pattern(regex_string_arg.getValue(), regex_flags);

  // Translation frame setup
  TranslationFrames tframe;
  switch(frame) {
  case 0:
    tframe = SINGLE_FRAME;
    break;
  case 1:
    tframe = WITH_REVERSE_COMPLEMENT;
    break;
  case 2:
    tframe = WITH_FRAME_SHIFTS;
    break;
  case 3:
    tframe = SIX_FRAME;
    break;
  default:
    tframe = SINGLE_FRAME;
    break;
  }

  // Output file setup
  SeqFileOut out_handle(cout, Fasta());
  if(format == "fasta") {
    setFormat(out_handle, Fasta());
  } else if(format == "fastq") {
    setFormat(out_handle, Fastq());
  } else {
    cerr << "Unrecognized output format";
    return 1;
  }

  // Loop variables
  CharString id;
  CharString seq;              // CharString more flexible than Dna5String
  CharString qual;
  SeqFileInWrapper seq_handle;

  // Loop over input files
  int nmatched = 0;
  for(string& infile: infiles) {

    try {
        seq_handle.open(infile);
    } catch(Exception const &e) {
      cerr << "Could not open " << infile << endl;
      seq_handle.close();
      return 1;
    }
 
    bool matched = false;
    while(!seq_handle.atEnd()) {

      try {

        readRecord(id, seq, qual, seq_handle.sqh);
      } catch (Exception const &e) {

        cerr << "Error: " << e.what() << endl;
        seq_handle.close();
        close(out_handle);
        return 1;

      } // End try-catch for record reading.


      // Regex
       if(seq_regex) {

         if(regex_search(match_type, regex("a"))) {
           match_type = "frcR";
         }
         if(regex_search(match_type, regex("A"))) {
           match_type = "frcRt";
         }

         // Due to some quirks of Seqan, I have to do a number of format
         // conversions, so this isn't as elegant as I would like.
         // Also this assumes DNA, not RNA, even though RNA could work
         // fine. Note that any type of sequence will work with regular
         // forward matching.
         matched = false;
         for(char& c: match_type) {
           switch (c) {
           case 'f':
             {
               CharString _seq = seq;
               matched |= regex_search(toCString(_seq), regex_pattern,
                                       regex_match_flags);
               break;
             }
           case 'r':
             {
               ModifiedString<CharString, ModReverse> rseq(seq);
               CharString _seq(rseq);
               matched |= regex_search(toCString(_seq), regex_pattern,
                                       regex_match_flags);
               break;
             }
           case 'c':
             {
               Dna5String dseq(seq);
               complement(dseq);
               CharString _seq(dseq);
               matched |= regex_search(toCString(_seq), regex_pattern,
                                       regex_match_flags);
               break;
             }
            case 'R':
             {
               Dna5String dseq(seq);
               reverseComplement(dseq);
               CharString _seq(dseq);
               matched |= regex_search(toCString(_seq), regex_pattern,
                                       regex_match_flags);
               break;
             }
           case 't':
             {
               // template<typename T> trans_search(seq, pattern) ... 
               // use with <Dna5String> for DNA or <Rna5String>...
               StringSet< String<AminoAcid> > aseqs;
               Dna5String dseq(seq);
               translate(aseqs, dseq, tframe);
               // Loop over translation frames
               for(String<AminoAcid>& _aseq: aseqs) {
                 CharString _seq(_aseq);
                 matched |= regex_search(toCString(_seq), regex_pattern,
                                        regex_match_flags);
              } // End loop over translation frames
              break;
            }

          } // End switch statement

        } // End match_type loop

      } else {

        // Simple regex on sequence IDs
        matched = regex_search(toCString(id), regex_pattern,
                               regex_match_flags);

      } // End regex

      // Write out if matched
      if((matched && !inverted) || (!matched && inverted)) {
        nmatched++;
        try {
            writeRecord(out_handle, id, seq, qual);
        } catch (Exception const &e) {
            cerr << "Error: " << e.what() << endl;
            cerr << "Error writing output" << endl;
            seq_handle.close();
            return 1;
        }
      } // End write out if matched

      
    } // End single file reading loop

    
    // Close the input handle and check for errors
    if(!seq_handle.close()) {
        cerr << "Problem closing " << infile << endl;
        close(out_handle);
        return 1;
    }

  } // End loop over files

  close(out_handle);

  if(nmatched) {
    return 0;
  } else {
    return 1;
  }

}
