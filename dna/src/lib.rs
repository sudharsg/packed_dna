//! A general-purpose genomics crate for dealing with DNA.

#![warn(missing_docs)]

use std::{convert::TryFrom, fmt::Display, str::FromStr};

// TODO: add a packed module with the PackedDna struct
//
// this struct must have the following:
// 1. A representation that is more memory efficient that simply storing a vector of `Nuc`
// 2. A FromStr implementation (should be case insensitive like the `Nuc` impl)
// 3. A `FromIterator` implementation to construct it from an iterator over `Nuc`s
// 4. A `fn get(&self, idx: usize) -> Nuc` getter for a particular nucleotide
//
// Make sure to unit test and document all elements
// Also, the internal representation of the PackedDna struct should be privately scoped

pub fn bits_char_convert(value: u8) -> char {
    match value {
        0u8 => 'A',
        1u8 => 'C',
        2u8 => 'G',
        3u8 => 'T',
        _ => 'X',
    }
}

pub fn char_bits_convert(value: char) -> u8 {
    match value.to_ascii_uppercase() {
        'A' => 0u8,
        'C' => 1u8,
        'G' => 2u8,
        'T' => 3u8,
        _ => 8u8,
    }
}

#[derive(Debug)]
pub struct PackedDna {
   // dna_data: String, 
    encode_data: u64, 
    en_size: u8
}

impl PackedDna {
    // pub fn get(&self, idx:usize) -> Nuc {
    //     let val = ((self.encode_data >> (idx*2)) & (3u64)) as u8;
    //     let c = bits_char_convert(val);
    //     return Nuc::c;
    // }
    pub fn new(dna_data: String) -> PackedDna {
        let mut en_size = 0u8;
        let mut encode_data = 0u64;
        for c in dna_data.chars(){
            let val = char_bits_convert(c);
            if val != 8u8 {
                encode_data = (encode_data << 2) | (val as u64);
                en_size += 1;
            }
        }
        PackedDna {encode_data, en_size}
    }
    pub fn iterate(&self) {
        let (mut a, mut c, mut g, mut t) = (0,0,0,0);
        for i in 0..self.en_size {
            let val = (self.encode_data >> (i*2)) & (3u64);
            if val == 0  { 
                a+= 1; 
            } else if val == 1 {
                c+=1;
            } else if val == 2 { 
                g+=1;
            } else if val == 3 { 
                t+=1;
            }
        }
        println!("A: {}", a);
        println!("C: {}", c);
        println!("G: {}", g);
        println!("T: {}", t);
    }
    
}

// impl FromIterator for PackedDna{
//     type Item = char;
//     fn next(&mut self) -> Option<char> {
//         let data = match self.decoded_rep {
//             'A' => self.dna.A,

//         }
//     }
// }

/// A nucleotide
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Nuc {
    /// Adenine
    A,
    /// Cytosine
    C,
    /// Guanine
    G,
    /// Thymine
    T,
}

/// An error that can occur when parsing a nucleotide.
#[derive(Debug, thiserror::Error)]
#[error("failed to parse nucleotide from {0}")]
pub struct ParseNucError<T: Display>(T);

impl TryFrom<char> for Nuc {
    type Error = ParseNucError<char>;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value.to_ascii_uppercase() {
            'A' => Ok(Self::A),
            'C' => Ok(Self::C),
            'G' => Ok(Self::G),
            'T' => Ok(Self::T),
            _ => Err(ParseNucError(value)),
        }
    }
}

impl FromStr for Nuc {
    type Err = ParseNucError<String>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let upper = s.to_ascii_uppercase();
        match upper.as_str() {
            "A" => Ok(Self::A),
            "C" => Ok(Self::C),
            "G" => Ok(Self::G),
            "T" => Ok(Self::T),
            _ => Err(ParseNucError(upper)),
        }
    }
}

#[cfg(test)]
mod tests {
    // TODO: fill in tests

    #[test]
    fn tryfrom_char() {
        assert!(false);
    }

    #[test]
    fn fromstr() {
        assert!(false);
    }
}
