use super::{SignalModel, BIN_SIZE_LRT, BIN_SIZE_SPEC_DIFF, BIN_SIZE_SPEC_FLAT};

pub const HISTOGRAM_SIZE: usize = 1000;
pub type Histogram = [i32; HISTOGRAM_SIZE];

#[derive(Copy, Clone)]
pub struct Histograms {
    lrt: Histogram,
    spectral_flatness: Histogram,
    spectral_diff: Histogram,
}

impl Histograms {
    pub const fn new() -> Self {
        Self {
            lrt: [0; HISTOGRAM_SIZE],
            spectral_flatness: [0; HISTOGRAM_SIZE],
            spectral_diff: [0; HISTOGRAM_SIZE],
        }
    }

    pub fn clear(&mut self) {
        self.lrt.fill(0);
        self.spectral_flatness.fill(0);
        self.spectral_diff.fill(0);
    }

    pub fn lrt(&self) -> &Histogram {
        &self.lrt
    }

    pub fn spectral_flatness(&self) -> &Histogram {
        &self.spectral_flatness
    }

    pub fn spectral_diff(&self) -> &Histogram {
        &self.spectral_diff
    }

    pub fn update(&mut self, features: &SignalModel) {
        fn update_inner(hist: &mut Histogram, val: f32, one_by: f32) {
            if val >= 0. {
                if let Some(x) = hist.get_mut((one_by * val) as usize) {
                    *x += 1;
                }
            }
        }
        update_inner(&mut self.lrt, features.lrt, 1. / BIN_SIZE_LRT);
        update_inner(
            &mut self.spectral_flatness,
            features.spectral_flatness,
            1. / BIN_SIZE_SPEC_FLAT,
        );
        update_inner(&mut self.spectral_diff, features.spectral_diff, 1. / BIN_SIZE_SPEC_DIFF);
    }
}
