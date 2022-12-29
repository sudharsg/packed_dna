//! A general-purpose genomics crate for dealing with DNA.

#![warn(missing_docs)]

use std::{convert::TryFrom, fmt::Display, str::FromStr};
use std::iter::FromIterator;


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


#[derive(Debug)]
pub struct PackedDna(Vec<u8>);

impl PackedDna {
    fn new() -> PackedDna {
        PackedDna(Vec::new())
    }

    fn add(&mut self, elem: u8) {
        self.0.push(elem);
    }

    fn get(&self, idx:usize) -> Nuc {
        return PackedDna::bits_enum_convert(self.0[idx]);
    }

    fn bits_char_convert(value: u8) -> char {
        match value {
            0u8 => 'A',
            1u8 => 'C',
            2u8 => 'G',
            3u8 => 'T',
            _ => 'X',
        }
    }

    fn char_bits_convert(value: char) -> u8 {
        match value.to_ascii_uppercase() {
            'A' => 0u8,
            'C' => 1u8,
            'G' => 2u8,
            'T' => 3u8,
            _ => 8u8,
        }
    }

    fn enum_bits_convert(value: Nuc) -> u8 {
        match value {
            Nuc::A => 0u8,
            Nuc::C => 1u8,
            Nuc::G => 2u8,
            Nuc::T => 3u8,
            // _ => 8u8,
        }
    }

    fn bits_enum_convert(value: u8) -> Nuc {
        match value {
            0u8 => Nuc::A,
            1u8 => Nuc::C,
            2u8 => Nuc::G,
            3u8 => Nuc::T,
            _ => Nuc::T,
        }
    }

    pub fn print_data(&self) {
        let (mut a, mut c, mut g, mut t) = (0,0,0,0);
        let data_size = self.0[0];
        let mut data = 0u8;
        // println!("data_size -- {:?}", data_size);
        for inx in 0..data_size{
            let i_index = inx as usize;
            if i_index%4 == 0 {
                data = self.0[(i_index+4)/4];
            }
            let item = (data >> ((i_index%4)*2)) & (3u8);
            let val = PackedDna::bits_char_convert(item);
            if val == 'A' { 
                a+= 1; 
            } else if val == 'C'  {
                c+=1;
            } else if val == 'G'  { 
                g+=1;
            } else if val == 'T' { 
                t+=1;
            }
        }
        println!("A: {}", a);
        println!("C: {}", c);
        println!("G: {}", g);
        println!("T: {}", t);
    }

}

impl FromIterator<Nuc> for PackedDna {
    fn from_iter<I: IntoIterator<Item=Nuc>>(iter: I) -> Self {
        let mut arr = PackedDna::new();
        let mut size = 0u8;
        let mut local_itr = 0u8;
        let mut local_data = 0u8;
        for nuc_data in iter {
            let val = PackedDna::enum_bits_convert(nuc_data);
            if local_itr == 4u8{
                arr.add(local_data);
                local_itr = 0u8;
                local_data = 0u8;
            }
            local_data = (local_data << 2) | (val as u8);
            local_itr += 1;
            // arr.add();
            size += 1;
        }
        if local_itr != 0u8 {
            arr.add(local_data);
        }
        arr.0.insert(0,size);
        arr
    }
}

impl FromStr for PackedDna {
    type Err = ParseNucError<String>;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let dna_data = s.to_ascii_uppercase();
        let mut arr = PackedDna::new();
        let mut size = 0u8;
        let mut local_itr = 0u8;
        let mut local_data = 0u8;
        for c in dna_data.chars(){
            let val = PackedDna::char_bits_convert(c);
            if val == 8u8 {
                continue;
            }
            if local_itr == 4u8{
                arr.add(local_data);
                local_itr = 0u8;
                local_data = 0u8;
            }
            local_data = (local_data << 2) | (val as u8);
            local_itr += 1;
            size +=1;
        }
        if local_itr != 0u8 {
            arr.add(local_data);
        }
        arr.0.insert(0,size);
        Ok(arr)
    }
}

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
