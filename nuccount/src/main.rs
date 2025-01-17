// TODO: implement a nucleotide counter
//
// command line argument parsing has been provided
// you must use the PackedDna struct you previously implemented
// if there is any functionality you would like to add to PackedDna feel free to do so in the DNA
// crate
//
// If run with `nuccount --dna ACGTTT" it should print the following to stdout:
// ```
// Input: ACGTTT
//
// A: 1
// C: 1
// G: 1
// T: 3
// ```
//
// be sure to exit with informative error messages if the input is invalid

use dna::PackedDna;
use std::str::FromStr;
use structopt::StructOpt;
// These need to be imported if need to use from_iter construct function
// of PackedDNA struct
// use dna::Nuc;
// use std::iter::FromIterator;

/// Count the number of occurrences of each nucleotide in the provided DNA.
#[derive(Debug, StructOpt)]
struct Opts {
    /// The DNA sequence for which we should retrieve a nucleotide count.
    ///
    /// It is case insensitive but only nucleotides A, C, G and T are supported.
    #[structopt(short = "d", long, required = true)]
    dna: String,
}

fn main() {
    let opts = Opts::from_args();
    let dna1 = opts.dna;
    println!("Input: {}", &dna1);
    // let nu_data = vec![Nuc::A, Nuc::C, Nuc::G, Nuc::T, Nuc::T, Nuc::T, Nuc::G];
    // let c = PackedDna::from_iter(vec![]);
    // c.print_data();

    // calling the from str constructor from DNA crate to build the
    // PackedDNA struct based on input strings
    let d = PackedDna::from_str(&dna1);
    // prints the frequencies of the nucleotides present in the input string
    d.expect("REASON").print_data();
}
