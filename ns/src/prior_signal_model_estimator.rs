use super::{
    Histogram, Histograms, PriorSignalModel, BIN_SIZE_LRT, BIN_SIZE_SPEC_DIFF, BIN_SIZE_SPEC_FLAT,
    FEATURE_UPDATE_WINDOW_SIZE,
};

#[derive(Copy, Clone)]
pub struct PriorSignalModelEstimator {
    prior_model: PriorSignalModel,
}

impl PriorSignalModelEstimator {
    pub const fn new(lrt_initial_value: f32) -> Self {
        let prior_model = PriorSignalModel::new(lrt_initial_value);
        Self { prior_model }
    }

    pub fn prior_model(&self) -> &PriorSignalModel {
        &self.prior_model
    }

    pub fn update(&mut self, histograms: &Histograms) {
        let low_lrt_fluctuations;
        (low_lrt_fluctuations, self.prior_model.lrt) = update_lrt(histograms.lrt());

        // For spectral flatness and spectral difference: compute the main peaks of
        // the histograms.
        let (spectral_flatness_peak_position, spectral_flatness_peak_weight) =
            find_first_of_two_larger_peaks(BIN_SIZE_SPEC_FLAT, histograms.spectral_flatness());
        let (spectral_diff_peak_position, spectral_diff_peak_weight) =
            find_first_of_two_larger_peaks(BIN_SIZE_SPEC_DIFF, histograms.spectral_diff());

        // Reject if weight of peaks is not large enough, or peak value too small.
        // Peak limit for spectral flatness (varies between 0 and 1).
        let use_spec_flat =
            !(spectral_flatness_peak_weight < 150 || spectral_flatness_peak_position < 0.6) as i32;

        // Reject if weight of peaks is not large enough or if fluctuation of the LRT
        // feature are very low, indicating a noise state.
        let use_spec_diff = !(spectral_diff_peak_weight < 150 || low_lrt_fluctuations) as i32;

        // Update the model.
        self.prior_model.template_diff_threshold =
            (1.2 * spectral_diff_peak_position).clamp(0.16, 1.);

        let one_by_feature_sum = 1. / (1. + use_spec_flat as f32 + use_spec_diff as f32);
        self.prior_model.lrt_weighting = one_by_feature_sum;

        self.prior_model.flatness_weighting = if use_spec_flat == 1 {
            self.prior_model.flatness_threshold =
                (0.9 * spectral_flatness_peak_position).clamp(0.1, 0.95);
            one_by_feature_sum
        } else {
            0.
        };

        self.prior_model.difference_weighting =
            if use_spec_diff == 1 { one_by_feature_sum } else { 0. }
    }
}

fn update_lrt(lrt_histogram: &Histogram) -> (bool, f32) {
    let (mut average, mut average_compl, mut average_squared, mut count) = (0., 0., 0., 0);

    for (i, &lrt) in lrt_histogram[..10].iter().enumerate() {
        average += lrt as f32 * bin_mid(i, BIN_SIZE_LRT);
        count += lrt;
    }

    if count > 0 {
        average /= count as f32;
    }

    for (i, &lrt) in lrt_histogram.iter().enumerate() {
        let bin_mid = bin_mid(i, BIN_SIZE_LRT);
        average_squared += lrt as f32 * bin_mid * bin_mid;
        average_compl += lrt as f32 * bin_mid;
    }

    const ONE_FEATURE_UPDATE_WINDOW_SIZE: f32 = 1. / FEATURE_UPDATE_WINDOW_SIZE as f32;
    average_squared *= ONE_FEATURE_UPDATE_WINDOW_SIZE;
    average_compl *= ONE_FEATURE_UPDATE_WINDOW_SIZE;

    // Fluctuation limit of LRT feature.
    let low_lrt_fluctuations = average_squared - average * average_compl < 0.05;

    // Get threshold for LRT feature.
    const MAX_LRT: f32 = 1.;
    const MIN_LRT: f32 = 0.2;

    (
        low_lrt_fluctuations,
        if low_lrt_fluctuations {
            // Very low fluctuation, so likely noise.
            MAX_LRT
        } else {
            (1.2 * average).clamp(MIN_LRT, MAX_LRT)
        },
    )
}

// Identifies the first of the two largest peaks in the histogram.
fn find_first_of_two_larger_peaks(bin_size: f32, spectral_flatness: &Histogram) -> (f32, i32) {
    let (mut peak_value, mut secondary_peak_value) = (0, 0);
    let (mut peak_position, mut secondary_peak_position) = (0., 0.);
    let (mut peak_weight, mut secondary_peak_weight) = (0, 0);

    // Identify the two largest peaks.
    for (i, &spectral_flatness) in spectral_flatness.iter().enumerate() {
        let bin_mid = bin_mid(i, bin_size);
        if spectral_flatness > peak_value {
            // Found new "first" peak candidate.
            secondary_peak_value = peak_value;
            secondary_peak_weight = peak_weight;
            secondary_peak_position = peak_position;

            peak_value = spectral_flatness;
            peak_weight = spectral_flatness;
            peak_position = bin_mid;
        } else if spectral_flatness > secondary_peak_value {
            // Found new "second" peak candidate.
            secondary_peak_value = spectral_flatness;
            secondary_peak_weight = spectral_flatness;
            secondary_peak_position = bin_mid;
        }
    }

    // Merge the peaks if they are close.
    if ((secondary_peak_position - peak_position).abs() < 2. * bin_size)
        && (2 * secondary_peak_weight > peak_weight)
    {
        peak_weight += secondary_peak_weight;
        peak_position = 0.5 * (peak_position + secondary_peak_position);
    }

    (peak_position, peak_weight)
}

fn bin_mid(i: usize, bin_size: f32) -> f32 {
    (i as f32 + 0.5) * bin_size
}
