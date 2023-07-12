#![no_std]

mod alaw;
mod tables;
mod ulaw;

pub use alaw::{alaw_to_linear, linear_to_alaw};
pub use ulaw::{linear_to_ulaw, ulaw_to_linear};

#[inline]
pub fn alaw_to_ulaw(alaw: u8) -> u8 {
    tables::ALAW_TO_ULAW[alaw as usize]
}

#[inline]
pub fn ulaw_to_alaw(ulaw: u8) -> u8 {
    tables::ULAW_TO_ALAW[ulaw as usize]
}

pub fn encode_alaw(samples: &[i16], output: &mut [u8]) {
    assert_eq!(samples.len(), output.len());
    assert_eq!(samples.as_ptr() as usize % 64, 0);
    assert_eq!(samples.len() % 64, 0);
    assert_eq!(output.as_ptr() as usize % 64, 0);
    assert_eq!(output.len() % 64, 0);
    output.iter_mut().zip(samples).for_each(|(o, &i)| *o = linear_to_alaw(i))
}
