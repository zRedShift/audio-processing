use super::{
    Spectrum, SuppressionParams, FFT_SIZE_BY_2_PLUS_1, LONG_STARTUP_PHASE_BLOCKS,
    ONE_BY_SHORT_STARTUP_PHASE_BLOCKS, SHORT_STARTUP_PHASE_BLOCKS,
};

// Estimates a Wiener-filter based frequency domain noise reduction filter.
#[derive(Copy, Clone)]
pub struct WienerFilter {
    suppression_params: &'static SuppressionParams,
    spectrum_prev_process: Spectrum,
    initial_spectral_estimate: Spectrum,
    filter: Spectrum,
}

impl WienerFilter {
    pub const fn new(suppression_params: &'static SuppressionParams) -> Self {
        Self {
            suppression_params,
            spectrum_prev_process: [0.; FFT_SIZE_BY_2_PLUS_1],
            initial_spectral_estimate: [0.; FFT_SIZE_BY_2_PLUS_1],
            filter: [1.; FFT_SIZE_BY_2_PLUS_1],
        }
    }

    pub fn filter(&self) -> &Spectrum {
        &self.filter
    }

    pub fn update(
        &mut self,
        num_analyzed_frames: i32,
        noise_spectrum: &Spectrum,
        prev_noise_spectrum: &Spectrum,
        parametric_noise_spectrum: &Spectrum,
        signal_spectrum: &Spectrum,
    ) {
        for i in 0..FFT_SIZE_BY_2_PLUS_1 {
            // Previous estimate based on previous frame with gain filter.
            let prev_tsa =
                self.spectrum_prev_process[i] / (prev_noise_spectrum[i] + 0.0001) * self.filter[i];

            // Current estimate.
            let current_tsa = if signal_spectrum[i] > noise_spectrum[i] {
                signal_spectrum[i] / (noise_spectrum[i] + 0.0001) - 1.
            } else {
                0.
            };

            // Directed decision estimate is sum of two terms: current estimate and
            // previous estimate.
            let snr_prior = 0.98 * prev_tsa + (1. - 0.98) * current_tsa;
            self.filter[i] = (snr_prior
                / (self.suppression_params.over_subtraction_factor + snr_prior))
                .clamp(self.suppression_params.minimum_attenuating_gain, 1.);
        }

        if num_analyzed_frames < SHORT_STARTUP_PHASE_BLOCKS {
            for i in 0..FFT_SIZE_BY_2_PLUS_1 {
                self.initial_spectral_estimate[i] += signal_spectrum[i];
                let mut filter_initial = self.initial_spectral_estimate[i]
                    - self.suppression_params.over_subtraction_factor
                        * parametric_noise_spectrum[i];
                filter_initial /= self.initial_spectral_estimate[i] + 0.0001;
                filter_initial =
                    filter_initial.clamp(self.suppression_params.minimum_attenuating_gain, 1.);

                // Weight the two suppression filters.
                filter_initial *= (SHORT_STARTUP_PHASE_BLOCKS - num_analyzed_frames) as f32;
                self.filter[i] *= num_analyzed_frames as f32;
                self.filter[i] += filter_initial;
                self.filter[i] *= ONE_BY_SHORT_STARTUP_PHASE_BLOCKS;
            }
        }

        self.spectrum_prev_process = *signal_spectrum;
    }

    pub fn compute_overall_scaling_factor(
        &self,
        num_analyzed_frames: i32,
        prior_speech_probability: f32,
        energy_before_filtering: f32,
        energy_after_filtering: f32,
    ) -> f32 {
        if !self.suppression_params.use_attenuation_adjustment
            || num_analyzed_frames <= LONG_STARTUP_PHASE_BLOCKS
        {
            return 1.;
        }

        let gain = (energy_after_filtering / (energy_before_filtering + 1.)).sqrt();

        // Scaling for new version. Threshold in final energy gain factor calculation.
        const B_LIM: f32 = 0.5;
        let scale_factor1 = if gain > B_LIM {
            let factor = 1. + 1.3 * (gain - B_LIM);
            if gain * factor > 1. {
                1. / gain
            } else {
                factor
            }
        } else {
            1.
        };

        let scale_factor2 = if gain < B_LIM {
            // Do not reduce scale too much for pause regions: attenuation here should
            // be controlled by flooring.
            1. - 0.3 * (B_LIM - gain.max(self.suppression_params.minimum_attenuating_gain))
        } else {
            1.
        };

        // Combine both scales with speech/noise prob: note prior
        // (prior_speech_probability) is not frequency dependent.
        prior_speech_probability * scale_factor1 + (1. - prior_speech_probability) * scale_factor2
    }
}
