#[derive(Copy, Clone, Debug)]
pub enum SuppressionLevel {
    Level6dB,
    Level12dB,
    Level18dB,
    Level21dB,
}

#[derive(Copy, Clone, Debug)]
pub enum SampleRate {
    Rate16kHz = 0,
    Rate32kHz = 1,
    Rate48kHz = 2,
}

impl SampleRate {
    pub const fn num_bands_minus_1(self) -> usize {
        self as _
    }
}
