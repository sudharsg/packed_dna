//! A general-purpose genomics crate for dealing with DNA.

#![warn(missing_docs)]

use std::{convert::TryFrom, fmt::Display, str::FromStr, iter::FromIterator, process};
// use std::iter::FromIterator;
// use std::

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
pub struct PackedDna(Vec<u32>);

impl PackedDna {
    fn new() -> PackedDna {
        PackedDna(Vec::new())
    }

    fn add(&mut self, elem: u32) {
        self.0.push(elem);
    }

    fn get(&self, idx:usize) -> Nuc {
        let data = self.0[(idx+16)/16];
        let item = ((data >> ((idx%16)*2)) & (3u32)) as u8;
        return PackedDna::bits_enum_convert(item);
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
        if data_size == 0 {
            print!("Input DNA sequence is empty; ");
            println!("Please enter a valid sequence using {{A,C,G,T}}");
            process::exit(1);   
        }
        for inx in 0..data_size{
            let i_index = inx as usize;
            let val = self.get(i_index);
            if val == Nuc::A { 
                a+= 1; 
            } else if val == Nuc::C  {
                c+=1;
            } else if val == Nuc::G  { 
                g+=1;
            } else if val == Nuc::T { 
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
        let mut size = 0u32;
        let mut local_data = 0u32;
        arr.add(0u32); // size of the data
        for nuc_data in iter {
            let val = PackedDna::enum_bits_convert(nuc_data);
            if (size%16u32 == 0) && (size != 0u32){
                arr.add(local_data);
                local_data = 0u32;
            }
            local_data = (local_data << 2) | (val as u32);
            size += 1;
        }
        if size%16u32 != 0u32 {
            arr.add(local_data);
        }
        arr.0[0] = size;
        arr
    }
}

impl FromStr for PackedDna {
    type Err = ParseNucError<String>;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let dna_data = s.to_ascii_uppercase();
        let mut arr = PackedDna::new();
        let mut size = 0u32;
        let mut local_data = 0u32;
        arr.add(0u32); // size of the data
        for c in dna_data.chars(){
            // input validity check
            if let Err(_parse_nuc_err) = Nuc::try_from(c){
                println!("Invalid Input - {:?} ; Please correct and rerun", c);
                process::exit(1);
            }
            let val = PackedDna::char_bits_convert(c);
            if (size%16u32 == 0) && (size != 0u32){
                arr.add(local_data);
                local_data = 0u32;
            }
            local_data = (local_data << 2) | (val as u32);
            size +=1;
        }
        if size%16u32 != 0u32 {
            arr.add(local_data);
        }
        arr.0[0] = size;
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
