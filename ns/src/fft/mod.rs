mod fft4g;

use super::{ComplexSpectrum, Signal, FFT_SIZE, FFT_SIZE_BY_2};

pub trait Fft {
    fn real_dft<'a>(&mut self, signal: &'a mut Signal) -> &'a mut ComplexSpectrum;
    fn inverse_real_dft<'a>(&mut self, spectrum: &'a mut ComplexSpectrum) -> &'a mut Signal;
}

fft::complex_dft_adapter!(pub ComplexDftAdapter, FFT_SIZE, true);

impl<F: FnMut(&mut ComplexSpectrum)> Fft for ComplexDftAdapter<F> {
    fn real_dft<'a>(&mut self, signal: &'a mut Signal) -> &'a mut ComplexSpectrum {
        self.real_dft(signal, true)
    }

    fn inverse_real_dft<'a>(&mut self, spectrum: &'a mut ComplexSpectrum) -> &'a mut Signal {
        self.inverse_real_dft(spectrum, true)
    }
}

pub struct Fft4g {
    bit_reversal_state: [usize; FFT_SIZE_BY_2],
    tables: [f32; FFT_SIZE_BY_2],
}

impl Fft4g {
    pub const fn new() -> Self {
        Self { bit_reversal_state: [0; FFT_SIZE_BY_2], tables: [0.; FFT_SIZE_BY_2] }
    }

    unsafe fn rdft(&mut self, inv: bool, time_data: &mut Signal) {
        fft4g::rdft(
            FFT_SIZE,
            if inv { -1 } else { 1 },
            time_data.as_mut_ptr(),
            self.bit_reversal_state.as_mut_ptr(),
            self.tables.as_mut_ptr(),
        )
    }
}

impl Fft for Fft4g {
    fn real_dft<'a>(&mut self, time_data: &'a mut Signal) -> &'a mut ComplexSpectrum {
        unsafe {
            self.rdft(false, time_data);
            &mut *(time_data as *mut _ as *mut _)
        }
    }

    fn inverse_real_dft<'a>(&mut self, time_data: &'a mut ComplexSpectrum) -> &'a mut Signal {
        unsafe {
            let time_data = &mut *(time_data as *mut _ as *mut _);
            self.rdft(true, time_data);

            // Scale the output
            const SCALING: f32 = 2. / FFT_SIZE as f32;

            time_data.iter_mut().for_each(|d| *d *= SCALING);
            time_data
        }
    }
}

use static_assertions::{assert_eq_align, assert_eq_size};

pub fn cast_complex_to_real(spectrum: &mut ComplexSpectrum) -> &mut Signal {
    assert_eq_size!(ComplexSpectrum, Signal);
    assert_eq_align!(ComplexSpectrum, Signal);
    unsafe { &mut *(spectrum as *mut _ as *mut _) }
}

pub fn cast_real_to_complex(signal: &mut Signal) -> &mut ComplexSpectrum {
    assert_eq_size!(Signal, ComplexSpectrum);
    assert_eq_align!(Signal, ComplexSpectrum);
    unsafe { &mut *(signal as *mut _ as *mut _) }
}
