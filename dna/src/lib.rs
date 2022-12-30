//! A general-purpose genomics crate for dealing with DNA.

#![warn(missing_docs)]

use std::{convert::TryFrom, fmt::Display, str::FromStr, iter::FromIterator, process};


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


#[derive(Debug, PartialEq)]
pub struct PackedDna{
    data: Vec<u8>,
    data_len: u32,
}

impl PackedDna {
    fn new(data:Vec<u8>, data_len:u32) -> PackedDna {
        PackedDna {data, data_len}
    }

    fn get(&self, idx:usize) -> Nuc {
        let data = self.data[(idx)/4];
        let item = (data >> ((idx%4)*2)) & (3u8);
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
        let data_size = self.data_len;
        if data_size == 0u32 {
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
        let mut arr = Vec::<u8>::new();
        let mut size = 0u32;
        let mut local_data = 0u8;
        for nuc_data in iter {
            let val = PackedDna::enum_bits_convert(nuc_data);
            if ((size%4u32) as u8 == 0u8) && (size != 0u32){
                arr.push(local_data);
                local_data = 0u8;
            }
            local_data = (local_data << 2) | (val as u8);
            size += 1;
        }
        if size%4u32 != 0u32 {
            arr.push(local_data);
        }
        PackedDna::new(arr, size)
    }
}

impl FromStr for PackedDna {
    type Err = ParseNucError<String>;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let dna_data = s.to_ascii_uppercase();
        let mut arr = Vec::<u8>::new();
        let mut size = 0u32;
        let mut local_data = 0u8;
        let mut err_data = Vec::new();
        for c in dna_data.chars(){
            // input validity check
            if let Err(_parse_nuc_err) = Nuc::try_from(c){
                err_data.push(c);
            }
            let val = PackedDna::char_bits_convert(c);
            if ((size%4u32) as u8 == 0u8) && (size != 0u32){
                arr.push(local_data);
                local_data = 0u8;
            }
            local_data = (local_data << 2) | (val as u8);
            size +=1;
        }
        if size%4u32 != 0u32 {
            arr.push(local_data);
        }
        // error handling
        if err_data.len() != 0 {
            println!("Invalid chars in input {:?}.\nPlease remove and rerun using only {{A,C,G,T}}",err_data);
            process::exit(1);
        }
        Ok(PackedDna::new(arr, size))
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
    use super::*;
    #[test]
    fn test_dna_from_iter() {
        assert_eq!(PackedDna::from_iter(vec![Nuc::A, Nuc::C, Nuc::G, Nuc::T, Nuc::T, Nuc::T, Nuc::G])
            ,PackedDna{data:vec![0b00011011,0b111110],data_len:7});
        assert_eq!(PackedDna::from_iter(vec![Nuc::A,Nuc::G,Nuc::C,Nuc::T,Nuc::G,Nuc::C,Nuc::T,Nuc::A,
            Nuc::G,Nuc::C,Nuc::T,Nuc::G,Nuc::A,Nuc::T,Nuc::C,Nuc::G,Nuc::A,Nuc::C])
            ,PackedDna{data:vec![0b00100111,0b10011100,0b10011110,0b00110110,0b0001],data_len:18});
        assert_eq!(PackedDna::from_iter(vec![]), PackedDna{data:vec![],data_len:0});
    }

    #[test]
    fn test_dna_from_str() {
        let res = PackedDna::from_str("ACGTTT").unwrap();
        assert_eq!(res,PackedDna{data:vec![0b00011011,0b1111],data_len:6});
        let res2 = PackedDna::from_str("AGCTGCTAGCTGATCGAAGTCAAAAAgggggtgAattttttttttttttttttttttgatgatcgtgacgtagtcgtacttagcta").unwrap();
        assert_eq!(res2
            ,PackedDna{data:vec![0b00100111,0b10011100,0b10011110,0b00110110,0b00001011,0b01000000,
                0b00001010,0b10101011,0b10000011,0b11111111,0b11111111,0b11111111,0b11111111,
                0b11111111,0b11100011,0b10001101,0b10111000,0b01101100,0b10110110,0b11000111,0b11001001,0b1100],data_len:86});
        let res3 = PackedDna::from_str("").unwrap();
        assert_eq!(res3, PackedDna{data:vec![],data_len:0});
    }

    #[test]
    fn tryfrom_char() {
        assert_eq!(Nuc::try_from('C').unwrap(), Nuc::C);
        assert_eq!(Nuc::try_from('A').unwrap(), Nuc::A);
        assert_eq!(Nuc::try_from('G').unwrap(), Nuc::G);
        assert_eq!(Nuc::try_from('T').unwrap(), Nuc::T);
    }

    #[test]
    fn fromstr() {
        assert_eq!(Nuc::from_str("C").unwrap(), Nuc::C);
        assert_eq!(Nuc::from_str("A").unwrap(), Nuc::A);
        assert_eq!(Nuc::from_str("G").unwrap(), Nuc::G);
        assert_eq!(Nuc::from_str("T").unwrap(), Nuc::T);
    }
}
