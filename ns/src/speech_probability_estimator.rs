use super::{SignalModelEstimator, Spectrum, FFT_SIZE_BY_2_PLUS_1, LONG_STARTUP_PHASE_BLOCKS};

#[derive(Copy, Clone)]
pub struct SpeechProbabilityEstimator {
    signal_model_estimator: SignalModelEstimator,
    prior_speech_prob: f32,
    speech_probability: Spectrum,
}

impl SpeechProbabilityEstimator {
    pub const fn new() -> Self {
        Self {
            signal_model_estimator: SignalModelEstimator::new(),
            prior_speech_prob: 0.5,
            speech_probability: [0.; FFT_SIZE_BY_2_PLUS_1],
        }
    }

    pub fn prior_probability(&self) -> f32 {
        self.prior_speech_prob
    }

    pub fn probability(&self) -> &Spectrum {
        &self.speech_probability
    }

    // Compute speech probability.
    #[allow(clippy::too_many_arguments)]
    pub fn update(
        &mut self,
        num_analyzed_frames: i32,
        prior_snr: &Spectrum,
        post_snr: &Spectrum,
        conservative_noise_spectrum: &Spectrum,
        signal_spectrum: &Spectrum,
        signal_spectral_sum: f32,
        signal_energy: f32,
    ) {
        if num_analyzed_frames < LONG_STARTUP_PHASE_BLOCKS {
            self.signal_model_estimator.adjust_normalization(num_analyzed_frames, signal_energy);
        }
        self.signal_model_estimator.update(
            prior_snr,
            post_snr,
            conservative_noise_spectrum,
            signal_spectrum,
            signal_spectral_sum,
            signal_energy,
        );

        let model = self.signal_model_estimator.model();
        let prior_model = self.signal_model_estimator.prior_model();

        // Compute indicator function: sigmoid map.
        fn sigmoid_map(first: f32, second: f32) -> f32 {
            // Width parameter in sigmoid map for prior model.
            const WIDTH_PRIOR_0: f32 = 4.;
            // Width for pause region: lower range, so increase width in tanh map.
            const WIDTH_PRIOR_1: f32 = 2. * WIDTH_PRIOR_0;

            let d = first - second;
            let tmp = (if d < 0. { WIDTH_PRIOR_1 } else { WIDTH_PRIOR_0 }) * d;
            0.5 * (tmp.tanh() + 1.)
        }

        // Average LRT feature: use larger width in tanh map for pause regions.
        let indicator0 = sigmoid_map(model.lrt, prior_model.lrt);

        // Spectral flatness feature: use larger width in tanh map for pause regions.
        let indicator1 = sigmoid_map(prior_model.flatness_threshold, model.spectral_flatness);

        // For template spectrum-difference : use larger width in tanh map for pause
        // regions.
        let indicator2 = sigmoid_map(model.spectral_diff, prior_model.template_diff_threshold);

        // Combine the indicator function with the feature weights.
        let ind_prior = prior_model.lrt_weighting * indicator0
            + prior_model.flatness_weighting * indicator1
            + prior_model.difference_weighting * indicator2;

        // Compute the prior probability.
        self.prior_speech_prob += 0.1 * (ind_prior - self.prior_speech_prob);

        // Make sure probabilities are within range: keep floor to 0.01.
        self.prior_speech_prob = self.prior_speech_prob.clamp(0.01, 1.);

        // Final speech probability: combine prior model with LR factor:.
        let gain_prior = (1. - self.prior_speech_prob) / (self.prior_speech_prob + 0.0001);

        for (&avg_log_lrt, speech_probability) in
            model.avg_log_lrt.iter().zip(self.speech_probability.iter_mut())
        {
            *speech_probability = 1. / (1. + gain_prior * (-avg_log_lrt).exp());
        }
    }
}
