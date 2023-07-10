use super::{
    QuantileNoiseEstimator, Spectrum, SuppressionParams, FFT_SIZE_BY_2_PLUS_1,
    ONE_BY_FFT_SIZE_BY_2_PLUS_1, ONE_BY_SHORT_STARTUP_PHASE_BLOCKS, SHORT_STARTUP_PHASE_BLOCKS,
};

// Struct for estimating the spectral characteristics of the noise in an incoming
// signal.
#[derive(Copy, Clone)]
pub struct NoiseEstimator {
    suppression_params: &'static SuppressionParams,
    white_noise_level: f32,
    pink_noise_numerator: f32,
    pink_noise_exp: f32,
    prev_noise_spectrum: Spectrum,
    conservative_noise_spectrum: Spectrum,
    parametric_noise_spectrum: Spectrum,
    noise_spectrum: Spectrum,
    quantile_noise_estimator: QuantileNoiseEstimator,
}

impl NoiseEstimator {
    pub const fn new(suppression_params: &'static SuppressionParams) -> Self {
        Self {
            suppression_params,
            white_noise_level: 0.,
            pink_noise_numerator: 0.,
            pink_noise_exp: 0.,
            prev_noise_spectrum: [0.; FFT_SIZE_BY_2_PLUS_1],
            conservative_noise_spectrum: [0.; FFT_SIZE_BY_2_PLUS_1],
            parametric_noise_spectrum: [0.; FFT_SIZE_BY_2_PLUS_1],
            noise_spectrum: [0.; FFT_SIZE_BY_2_PLUS_1],
            quantile_noise_estimator: QuantileNoiseEstimator::new(),
        }
    }

    pub fn prepare_analysis(&mut self) {
        self.prev_noise_spectrum = self.noise_spectrum;
    }

    pub fn pre_update(
        &mut self,
        num_analyzed_frames: i32,
        signal_spectrum: &Spectrum,
        signal_spectral_sum: f32,
    ) {
        self.quantile_noise_estimator.estimate(signal_spectrum, &mut self.noise_spectrum);

        if num_analyzed_frames >= SHORT_STARTUP_PHASE_BLOCKS {
            return;
        }

        // Compute simplified noise model during startup.
        const START_BAND: usize = 5;
        let mut sum_log_i_log_magn = 0.;
        let mut sum_log_i = 0.;
        let mut sum_log_i_square = 0.;
        let mut sum_log_magn = 0.;
        for i in START_BAND..FFT_SIZE_BY_2_PLUS_1 {
            let log_i = LOG_TABLE[i];
            sum_log_i += log_i;
            sum_log_i_square += log_i * log_i;
            let log_signal = signal_spectrum[i].ln();
            sum_log_magn += log_signal;
            sum_log_i_log_magn += log_i * log_signal;
        }

        // Estimate the parameter for the level of the white noise.
        self.white_noise_level += signal_spectral_sum
            * ONE_BY_FFT_SIZE_BY_2_PLUS_1
            * self.suppression_params.over_subtraction_factor;

        // Estimate pink noise parameters.
        let denom =
            sum_log_i_square * (FFT_SIZE_BY_2_PLUS_1 - START_BAND) as f32 - sum_log_i * sum_log_i;
        let mut num = sum_log_i_square * sum_log_magn - sum_log_i * sum_log_i_log_magn;
        debug_assert_ne!(denom, 0.);
        let mut pink_noise_adjustment = num / denom;

        // Constrain the estimated spectrum to be positive.
        pink_noise_adjustment = pink_noise_adjustment.max(0.);
        self.pink_noise_numerator += pink_noise_adjustment;
        num = sum_log_i * sum_log_magn
            - (FFT_SIZE_BY_2_PLUS_1 - START_BAND) as f32 * sum_log_i_log_magn;
        debug_assert_ne!(denom, 0.);
        pink_noise_adjustment = num / denom;

        // Constrain the pink noise power to be in the interval [0, 1].
        pink_noise_adjustment = pink_noise_adjustment.clamp(0., 1.);

        self.pink_noise_exp += pink_noise_adjustment;

        let one_by_num_analyzed_frames_plus_1 = 1. / (num_analyzed_frames + 1) as f32;

        // Calculate the frequency-independent parts of parametric noise estimate.
        let mut parametric_exp = 0.;
        let mut parametric_num = 0.;
        if self.pink_noise_exp > 0. {
            // Use pink noise estimate.
            parametric_num = (self.pink_noise_numerator * one_by_num_analyzed_frames_plus_1).exp();
            parametric_num *= (num_analyzed_frames + 1) as f32;
            parametric_exp = self.pink_noise_exp * one_by_num_analyzed_frames_plus_1;
        }

        for (i, (parametric_noise_spectrum, noise_spectrum)) in self
            .parametric_noise_spectrum
            .iter_mut()
            .zip(self.noise_spectrum.iter_mut())
            .enumerate()
        {
            // Estimate the background noise using the white and pink noise
            // parameters.
            *parametric_noise_spectrum = if self.pink_noise_exp == 0. {
                // Use white noise estimate.
                self.white_noise_level
            } else {
                // Use pink noise estimate.
                let use_band = if i < START_BAND { START_BAND } else { i };
                let denom = (use_band as f32).powf(parametric_exp);
                debug_assert_ne!(denom, 0.);
                parametric_num / denom
            };

            // Weight quantile noise with modeled noise.
            *noise_spectrum *= num_analyzed_frames as f32;
            let tmp = *parametric_noise_spectrum
                * (SHORT_STARTUP_PHASE_BLOCKS - num_analyzed_frames) as f32;
            *noise_spectrum += tmp * one_by_num_analyzed_frames_plus_1;
            *noise_spectrum *= ONE_BY_SHORT_STARTUP_PHASE_BLOCKS;
        }
    }

    pub fn post_update(&mut self, speech_probability: &Spectrum, signal_spectrum: &Spectrum) {
        // Time-avg parameter for noise_spectrum update.
        const NOISE_UPDATE: f32 = 0.9;

        let mut gamma = NOISE_UPDATE;
        for i in 0..FFT_SIZE_BY_2_PLUS_1 {
            let prob_speech = speech_probability[i];
            let prob_non_speech = 1. - prob_speech;

            // Temporary noise update used for speech frames if update value is less
            // than previous.
            let noise_update_tmp = gamma * self.prev_noise_spectrum[i]
                + (1. - gamma)
                    * (prob_non_speech * signal_spectrum[i]
                        + prob_speech * self.prev_noise_spectrum[i]);

            // Time-constant based on speech/noise_spectrum state.
            let gamma_old = gamma;

            // Increase gamma for frame likely to be speech.
            const PROB_RANGE: f32 = 0.2;
            gamma = if prob_speech > PROB_RANGE { 0.99 } else { NOISE_UPDATE };

            // Conservative noise_spectrum update.
            if prob_speech < PROB_RANGE {
                self.conservative_noise_spectrum[i] +=
                    0.05 * (signal_spectrum[i] - self.conservative_noise_spectrum[i]);
            }

            // Noise_spectrum update.
            self.noise_spectrum[i] = if gamma == gamma_old {
                noise_update_tmp
            } else {
                let noise_spectrum = gamma * self.prev_noise_spectrum[i]
                    + (1. - gamma)
                        * (prob_non_speech * signal_spectrum[i]
                            + prob_speech * self.prev_noise_spectrum[i]);
                // Allow for noise_spectrum update downwards: If noise_spectrum update
                // decreases the noise_spectrum, it is safe, so allow it to happen.
                noise_spectrum.min(noise_update_tmp)
            }
        }
    }

    pub fn noise_spectrum(&self) -> &Spectrum {
        &self.noise_spectrum
    }

    pub fn prev_noise_spectrum(&self) -> &Spectrum {
        &self.prev_noise_spectrum
    }

    pub fn parametric_noise_spectrum(&self) -> &Spectrum {
        &self.parametric_noise_spectrum
    }

    pub fn conservative_noise_spectrum(&self) -> &Spectrum {
        &self.conservative_noise_spectrum
    }
}

#[allow(clippy::approx_constant, clippy::excessive_precision)]
const LOG_TABLE: &Spectrum = &[
    0., 0., 0., 0., 0., 1.609438, 1.791759, 1.945910, 2.079442, 2.197225, 2.302585, 2.397895,
    2.484907, 2.564949, 2.639057, 2.708050, 2.772589, 2.833213, 2.890372, 2.944439, 2.995732,
    3.044522, 3.091043, 3.135494, 3.178054, 3.218876, 3.258097, 3.295837, 3.332205, 3.367296,
    3.401197, 3.433987, 3.465736, 3.496507, 3.526361, 3.555348, 3.583519, 3.610918, 3.637586,
    3.663562, 3.688879, 3.713572, 3.737669, 3.761200, 3.784190, 3.806663, 3.828641, 3.850147,
    3.871201, 3.891820, 3.912023, 3.931826, 3.951244, 3.970292, 3.988984, 4.007333, 4.025352,
    4.043051, 4.060443, 4.077538, 4.094345, 4.110874, 4.127134, 4.143135, 4.158883, 4.174387,
    4.189655, 4.204693, 4.219508, 4.234107, 4.248495, 4.262680, 4.276666, 4.290460, 4.304065,
    4.317488, 4.330733, 4.343805, 4.356709, 4.369448, 4.382027, 4.394449, 4.406719, 4.418841,
    4.430817, 4.442651, 4.454347, 4.465908, 4.477337, 4.488636, 4.499810, 4.510859, 4.521789,
    4.532599, 4.543295, 4.553877, 4.564348, 4.574711, 4.584968, 4.595119, 4.605170, 4.615121,
    4.624973, 4.634729, 4.644391, 4.653960, 4.663439, 4.672829, 4.682131, 4.691348, 4.700480,
    4.709530, 4.718499, 4.727388, 4.736198, 4.744932, 4.753591, 4.762174, 4.770685, 4.779124,
    4.787492, 4.795791, 4.804021, 4.812184, 4.820282, 4.828314, 4.836282, 4.844187, 4.852030,
];
