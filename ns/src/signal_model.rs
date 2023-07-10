use super::{Spectrum, FFT_SIZE_BY_2_PLUS_1, LRT_FEATURE_THR};

#[derive(Copy, Clone)]
pub struct SignalModel {
    pub lrt: f32,
    pub spectral_diff: f32,
    pub spectral_flatness: f32,
    pub avg_log_lrt: Spectrum,
}

impl SignalModel {
    pub const fn new() -> Self {
        const SF_FEATURE_THR: f32 = 0.5;
        Self {
            lrt: LRT_FEATURE_THR,
            spectral_diff: SF_FEATURE_THR,
            spectral_flatness: SF_FEATURE_THR,
            avg_log_lrt: [LRT_FEATURE_THR; FFT_SIZE_BY_2_PLUS_1],
        }
    }
}
