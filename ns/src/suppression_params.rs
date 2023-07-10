use super::SuppressionLevel::{self};

pub struct SuppressionParams {
    pub over_subtraction_factor: f32,
    pub minimum_attenuating_gain: f32,
    pub use_attenuation_adjustment: bool,
}

impl SuppressionParams {
    const LEVEL_6_DB: Self = Self::new(1., 0.5, false);
    const LEVEL_12_DB: Self = Self::new(1., 0.25, true);
    const LEVEL_18_DB: Self = Self::new(1.1, 0.125, true);
    const LEVEL_21_DB: Self = Self::new(1.25, 0.09, true);

    const fn new(
        over_subtraction_factor: f32,
        minimum_attenuating_gain: f32,
        use_attenuation_adjustment: bool,
    ) -> Self {
        Self { over_subtraction_factor, minimum_attenuating_gain, use_attenuation_adjustment }
    }

    pub const fn from_level(level: SuppressionLevel) -> &'static Self {
        match level {
            SuppressionLevel::Level6dB => &Self::LEVEL_6_DB,
            SuppressionLevel::Level12dB => &Self::LEVEL_12_DB,
            SuppressionLevel::Level18dB => &Self::LEVEL_18_DB,
            SuppressionLevel::Level21dB => &Self::LEVEL_21_DB,
        }
    }
}
