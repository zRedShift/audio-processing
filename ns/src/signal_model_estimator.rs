use super::{
    Histograms, PriorSignalModel, PriorSignalModelEstimator, SignalModel, Spectrum,
    FEATURE_UPDATE_WINDOW_SIZE, LRT_FEATURE_THR, ONE_BY_FFT_SIZE_BY_2_PLUS_1,
};

#[derive(Copy, Clone)]
pub struct SignalModelEstimator {
    diff_normalization: f32,
    signal_energy_sum: f32,
    histograms: Histograms,
    histogram_analysis_counter: i32,
    prior_model_estimator: PriorSignalModelEstimator,
    features: SignalModel,
}

impl SignalModelEstimator {
    pub const fn new() -> Self {
        Self {
            diff_normalization: 0.,
            signal_energy_sum: 0.,
            histograms: Histograms::new(),
            histogram_analysis_counter: FEATURE_UPDATE_WINDOW_SIZE,
            prior_model_estimator: PriorSignalModelEstimator::new(LRT_FEATURE_THR),
            features: SignalModel::new(),
        }
    }

    pub fn prior_model(&self) -> &PriorSignalModel {
        self.prior_model_estimator.prior_model()
    }

    pub fn model(&self) -> &SignalModel {
        &self.features
    }

    pub fn adjust_normalization(&mut self, num_analyzed_frames: i32, signal_energy: f32) {
        self.diff_normalization *= num_analyzed_frames as f32;
        self.diff_normalization += signal_energy;
        self.diff_normalization /= (num_analyzed_frames + 1) as f32;
    }

    pub fn update(
        &mut self,
        prior_snr: &Spectrum,
        post_snr: &Spectrum,
        conservative_noise_spectrum: &Spectrum,
        signal_spectrum: &Spectrum,
        signal_spectral_sum: f32,
        signal_energy: f32,
    ) {
        // Compute spectral flatness on input spectrum.
        update_spectral_flatness(
            signal_spectrum,
            signal_spectral_sum,
            &mut self.features.spectral_flatness,
        );
        // Compute difference of input spectrum with learned/estimated noise spectrum.
        let spectral_diff = compute_spectral_diff(
            conservative_noise_spectrum,
            signal_spectrum,
            signal_spectral_sum,
            self.diff_normalization,
        );
        // Compute time-avg update of difference feature.
        self.features.spectral_diff += 0.3 * (spectral_diff - self.features.spectral_diff);

        self.signal_energy_sum += signal_energy;

        // Compute histograms for parameter decisions (thresholds and weights for
        // features). Parameters are extracted periodically.
        self.histogram_analysis_counter -= 1;

        if self.histogram_analysis_counter > 0 {
            self.histograms.update(&self.features);
        } else {
            // Compute model parameters.
            self.prior_model_estimator.update(&self.histograms);

            // Clear histograms for next update.
            self.histograms.clear();

            self.histogram_analysis_counter = FEATURE_UPDATE_WINDOW_SIZE;

            // Update every window:
            // Compute normalization for the spectral difference for next estimation.
            self.signal_energy_sum /= FEATURE_UPDATE_WINDOW_SIZE as f32;
            self.diff_normalization = 0.5 * (self.signal_energy_sum + self.diff_normalization);
            self.signal_energy_sum = 0.;
        }

        // Compute the LRT.
        update_spectral_lrt(prior_snr, post_snr, &mut self.features);
    }
}

// Updates the spectral flatness based on the input spectrum.
fn update_spectral_flatness(
    signal_spectrum: &Spectrum,
    signal_spectral_sum: f32,
    spectral_flatness: &mut f32,
) {
    // Compute log of ratio of the geometric to arithmetic mean (handle the log(0)
    // separately).
    const AVERAGING: f32 = 0.3;
    let [first, tail @ ..] = signal_spectrum;
    if tail.iter().any(|&x| x == 0.) {
        *spectral_flatness -= AVERAGING * *spectral_flatness;
        return;
    }

    let avg_spect_flatness_num =
        tail.iter().map(|x| x.ln()).sum::<f32>() * ONE_BY_FFT_SIZE_BY_2_PLUS_1;
    let avg_spect_flatness_denom = (signal_spectral_sum - *first) * ONE_BY_FFT_SIZE_BY_2_PLUS_1;

    let spectral_tmp = avg_spect_flatness_num.exp() / avg_spect_flatness_denom;

    // Time-avg update of spectral flatness feature.
    *spectral_flatness = AVERAGING * (spectral_tmp - *spectral_flatness);
}

// Computes the difference measure between input spectrum and a template/learned
// noise spectrum.
fn compute_spectral_diff(
    conservative_noise_spectrum: &Spectrum,
    signal_spectrum: &Spectrum,
    signal_spectral_sum: f32,
    diff_normalization: f32,
) -> f32 {
    // spectral_diff = var(signal_spectrum) - cov(signal_spectrum, magnAvgPause)^2
    // / var(magnAvgPause)

    // Compute average quantities.
    // Conservative smooth noise spectrum from pause frames.
    let noise_average =
        conservative_noise_spectrum.iter().sum::<f32>() * ONE_BY_FFT_SIZE_BY_2_PLUS_1;
    let signal_average = signal_spectral_sum * ONE_BY_FFT_SIZE_BY_2_PLUS_1;

    // Compute variance and covariance quantities.
    let (mut covariance, mut noise_variance, mut signal_variance) = (0., 0., 0.);
    for (&signal, &noise) in signal_spectrum.iter().zip(conservative_noise_spectrum.iter()) {
        let signal_diff = signal - signal_average;
        let noise_diff = noise - noise_average;
        covariance += signal_diff * noise_diff;
        noise_variance += noise_diff * noise_diff;
        signal_variance += signal_diff * signal_diff;
    }
    covariance *= ONE_BY_FFT_SIZE_BY_2_PLUS_1;
    noise_variance *= ONE_BY_FFT_SIZE_BY_2_PLUS_1;
    signal_variance *= ONE_BY_FFT_SIZE_BY_2_PLUS_1;

    // Update of average magnitude spectrum.
    let spectral_diff = signal_variance - (covariance * covariance) / (noise_variance + 0.0001);
    // Normalize.
    spectral_diff / (diff_normalization + 0.0001)
}

// Updates the log LRT measures.
fn update_spectral_lrt(prior_snr: &Spectrum, post_snr: &Spectrum, features: &mut SignalModel) {
    for ((&prior, &post), avg_log_lrt) in
        prior_snr.iter().zip(post_snr.iter()).zip(features.avg_log_lrt.iter_mut())
    {
        let tmp1 = 1. + 2. * prior;
        let tmp2 = 2. * prior / (tmp1 + 0.0001);
        let bessel_tmp = (post + 1.) * tmp2;
        *avg_log_lrt += 0.5 * (bessel_tmp - tmp1.ln() - *avg_log_lrt);
    }

    features.lrt = features.avg_log_lrt.iter().sum::<f32>() * ONE_BY_FFT_SIZE_BY_2_PLUS_1;
}
