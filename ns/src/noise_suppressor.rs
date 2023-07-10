use super::{
    ComplexSpectrum, Fft, Frame, NoiseEstimator, Signal, Spectrum, SpeechProbabilityEstimator,
    SplitAudioFrames, SuppressionLevel, SuppressionParams, WienerFilter, FFT_SIZE, FFT_SIZE_BY_2,
    FFT_SIZE_BY_2_PLUS_1, NS_FRAME_SIZE, ONE_BY_FFT_SIZE_BY_2_PLUS_1, OVERLAP_SIZE,
};

type Memory = [f32; OVERLAP_SIZE];

const OVERLAP_REMAINDER: usize = NS_FRAME_SIZE - OVERLAP_SIZE;

// Struct for suppressing noise in a signal.
pub struct NoiseSuppressor<F, const B_M_1: usize = 0, const CH: usize = 1> {
    suppression_params: &'static SuppressionParams,
    num_analyzed_frames: i32,
    fft: F,
    channels: [ChannelState<B_M_1>; CH],
}

impl<F: Fft, const B_M_1: usize, const CH: usize> NoiseSuppressor<F, B_M_1, CH> {
    const ASSERT: () = if B_M_1 > 2 || CH < 1 {
        panic!("1 <= BANDS <= 3 && CHANNELS > 0")
    };
    const BANDS: usize = B_M_1 + 1;
    pub const fn new(fft: F, suppression_level: SuppressionLevel) -> Self {
        matches!(Self::ASSERT, ());
        let suppression_params = SuppressionParams::from_level(suppression_level);
        Self {
            suppression_params,
            num_analyzed_frames: -1,
            fft,
            channels: [ChannelState::<B_M_1>::new(suppression_params); CH],
        }
    }

    pub fn analyze<const BANDS: usize>(&mut self, audio: &SplitAudioFrames<BANDS, CH>) {
        assert_eq!(BANDS, Self::BANDS);
        self.analyze_with_scratch(audio, &mut [0.; FFT_SIZE])
    }

    pub fn analyze_with_scratch<const BANDS: usize>(
        &mut self,
        audio: &SplitAudioFrames<BANDS, CH>,
        extended_frame: &mut Signal,
    ) {
        assert_eq!(BANDS, Self::BANDS);
        for ch in self.channels.iter_mut() {
            ch.noise_estimator.prepare_analysis();
        }

        // Check for zero frames.
        let zero_frame = !self.channels.iter().zip(audio.iter()).any(|(ch, f)| {
            compute_energy_of_extended_frame_subcomponents(&f[0], &ch.analyze_analysis_memory) > 0.
        });

        if zero_frame {
            // We want to avoid updating statistics in this case:
            // Updating feature statistics when we have zeros only will cause
            // thresholds to move towards zero signal situations. This in turn has the
            // effect that once the signal is "turned on" (non-zero values) everything
            // will be treated as speech and there is no noise suppression effect.
            // Depending on the duration of the inactive signal it takes a
            // considerable amount of time for the system to learn what is noise and
            // what is speech.
            return;
        }

        self.num_analyzed_frames = self.num_analyzed_frames.checked_add(1).unwrap_or_default();

        for (ch, frame) in self.channels.iter_mut().zip(audio.iter()) {
            // Form an extended frame and apply analysis filter bank windowing.
            form_extended_frame(&frame[0], &mut ch.analyze_analysis_memory, extended_frame);
            apply_filter_bank_window(extended_frame);

            // Compute energies and the magnitude spectrum.
            let complex_spectrum = self.fft.real_dft(extended_frame);
            let signal_energy = compute_signal_energy(complex_spectrum);
            let signal_spectrum = &*compute_magnitude_spectrum_in_place(complex_spectrum);
            let signal_spectral_sum = signal_spectrum.iter().sum();

            // Estimate the noise spectra and the probability estimates of speech
            // presence.
            ch.noise_estimator.pre_update(
                self.num_analyzed_frames,
                signal_spectrum,
                signal_spectral_sum,
            );

            let post_snr = compute_post_snr(signal_spectrum, ch.noise_estimator.noise_spectrum());

            let prior_snr = compute_prior_snr(
                ch.wiener_filter.filter(),
                &ch.prev_analysis_signal_spectrum,
                ch.noise_estimator.prev_noise_spectrum(),
                &post_snr,
            );

            ch.speech_probability_estimator.update(
                self.num_analyzed_frames,
                &prior_snr,
                &post_snr,
                ch.noise_estimator.conservative_noise_spectrum(),
                signal_spectrum,
                signal_spectral_sum,
                signal_energy,
            );

            ch.noise_estimator
                .post_update(ch.speech_probability_estimator.probability(), signal_spectrum);

            // Store the magnitude spectrum to make it available for the process
            // method.
            ch.prev_analysis_signal_spectrum.copy_from_slice(signal_spectrum);
        }
    }

    pub fn process<const BANDS: usize>(&mut self, audio: &mut SplitAudioFrames<BANDS, CH>) {
        self.process_with_scratch(audio, &mut [[0.; FFT_SIZE]; CH], &mut [0.; FFT_SIZE_BY_2_PLUS_1])
    }

    pub fn process_with_scratch<const BANDS: usize>(
        &mut self,
        audio: &mut SplitAudioFrames<BANDS, CH>,
        extended_frames: &mut [Signal; CH],
        signal_spectrum: &mut Spectrum,
    ) {
        assert_eq!(B_M_1 + 1, BANDS);
        let mut energies_before_filtering = [0.; CH];
        let mut upper_band_gain = f32::MAX;
        let mut gain_adjustment = f32::MAX;

        // Compute the suppression filters for all channels.
        for (((extended_frame, ch), frame), energy_before_filtering) in extended_frames
            .iter_mut()
            .zip(&mut self.channels)
            .zip(&*audio)
            .zip(&mut energies_before_filtering)
        {
            form_extended_frame(&frame[0], &mut ch.process_analysis_memory, extended_frame);
            apply_filter_bank_window(extended_frame);
            *energy_before_filtering = compute_energy_of_extended_frame(extended_frame);

            // Perform filter bank analysis and compute the magnitude spectrum.
            let complex_spectrum = self.fft.real_dft(extended_frame);
            compute_magnitude_spectrum(complex_spectrum, signal_spectrum);

            // Compute the frequency domain gain filter for noise attenuation.
            ch.wiener_filter.update(
                self.num_analyzed_frames,
                ch.noise_estimator.noise_spectrum(),
                ch.noise_estimator.prev_noise_spectrum(),
                ch.noise_estimator.parametric_noise_spectrum(),
                signal_spectrum,
            );
            if BANDS > 1 {
                // Compute the time-domain gain for attenuating the noise in the upper
                // bands.
                let curr_gain = compute_upper_bands_gain(
                    self.suppression_params.minimum_attenuating_gain,
                    ch.wiener_filter.filter(),
                    ch.speech_probability_estimator.probability(),
                    &ch.prev_analysis_signal_spectrum,
                    signal_spectrum,
                );
                upper_band_gain = upper_band_gain.min(curr_gain);
            }
        }

        // Aggregate the Wiener filters for all channels.
        let [filter_first, filter @ .., filter_last] = if CH == 1 {
            self.channels[0].wiener_filter.filter()
        } else {
            self.aggregate_wiener_filters(signal_spectrum);
            signal_spectrum
        };

        for extended_frame in extended_frames.iter_mut() {
            let complex_spectrum = super::cast_real_to_complex(extended_frame);
            let [spectrum_first, spectrum @ ..] = complex_spectrum;

            // Apply the filter to the lower band.
            spectrum_first.re *= filter_first;
            spectrum_first.im *= filter_last;
            for (spectrum, &filter) in spectrum.iter_mut().zip(filter.iter()) {
                *spectrum *= filter;
            }

            // Perform filter bank synthesis
            self.fft.inverse_real_dft(complex_spectrum);
        }

        for ((extended_frame, ch), energy_before_filtering) in
            extended_frames.iter_mut().zip(&mut self.channels).zip(energies_before_filtering)
        {
            let energy_after_filtering = compute_energy_of_extended_frame(extended_frame);

            // Apply synthesis window.
            apply_filter_bank_window(extended_frame);

            // Compute the adjustment of the noise attenuation filter based on the
            // effect of the attenuation.
            let curr_adjustment = ch.wiener_filter.compute_overall_scaling_factor(
                self.num_analyzed_frames,
                ch.speech_probability_estimator.prior_probability(),
                energy_before_filtering,
                energy_after_filtering,
            );

            // Select the noise attenuating gain to apply to the upper band.
            gain_adjustment = gain_adjustment.min(curr_adjustment);
        }

        for ((extended_frame, ch), frame) in
            extended_frames.iter_mut().zip(&mut self.channels).zip(audio)
        {
            let (first, rest) = frame.split_first_mut().unwrap();
            let higher_bands: &mut [Frame; B_M_1] = rest.try_into().unwrap();

            // Select and apply adjustment of the noise attenuation filter based on the
            // effect of the attenuation.
            extended_frame.iter_mut().for_each(|f| *f *= gain_adjustment);

            // Use overlap-and-add to form the output frame of the lowest band.
            overlap_and_add(extended_frame, &mut ch.process_synthesis_memory, first);

            // Process the upper bands.
            for (memory, frame) in ch.process_delay_memory.iter_mut().zip(higher_bands) {
                // Delay the upper bands to match the delay of the filterbank applied to
                // the lowest band.
                delay_signal_in_place(frame, memory, extended_frame);

                // Apply the time-domain noise-attenuating gain.
                frame.iter_mut().for_each(|f| *f *= upper_band_gain);
            }

            // Limit the output the allowed range.
            frame.iter_mut().flatten().for_each(|f| *f = f.clamp(-32768., 32767.));
        }
    }

    // Aggregates the Wiener filters into a single filter to use.
    fn aggregate_wiener_filters(&self, filter: &mut Spectrum) {
        filter.copy_from_slice(self.channels[0].wiener_filter.filter());

        for ChannelState { wiener_filter, .. } in &self.channels[1..] {
            for (filter, &filter_ch) in filter.iter_mut().zip(wiener_filter.filter()) {
                *filter = filter.min(filter_ch);
            }
        }
    }
}

#[derive(Copy, Clone)]
struct ChannelState<const B_M_1: usize> {
    speech_probability_estimator: SpeechProbabilityEstimator,
    wiener_filter: WienerFilter,
    noise_estimator: NoiseEstimator,
    prev_analysis_signal_spectrum: Spectrum,
    analyze_analysis_memory: Memory,
    process_analysis_memory: Memory,
    process_synthesis_memory: Memory,
    process_delay_memory: [Memory; B_M_1],
}

impl<const B_M_1: usize> ChannelState<B_M_1> {
    const fn new(suppression_params: &'static SuppressionParams) -> Self {
        Self {
            speech_probability_estimator: SpeechProbabilityEstimator::new(),
            wiener_filter: WienerFilter::new(suppression_params),
            noise_estimator: NoiseEstimator::new(suppression_params),
            prev_analysis_signal_spectrum: [1.; FFT_SIZE_BY_2_PLUS_1],
            analyze_analysis_memory: [0.; OVERLAP_SIZE],
            process_analysis_memory: [0.; OVERLAP_SIZE],
            process_synthesis_memory: [0.; OVERLAP_SIZE],
            process_delay_memory: [[0.; OVERLAP_SIZE]; B_M_1],
        }
    }
}

#[allow(clippy::approx_constant, clippy::excessive_precision)]
const BLOCKS_160_W_256_FIRST_HALF: &Memory = &[
    0.00000000, 0.01636173, 0.03271908, 0.04906767, 0.06540313, 0.08172107, 0.09801714, 0.11428696,
    0.13052619, 0.14673047, 0.16289547, 0.17901686, 0.19509032, 0.21111155, 0.22707626, 0.24298018,
    0.25881905, 0.27458862, 0.29028468, 0.30590302, 0.32143947, 0.33688985, 0.35225005, 0.36751594,
    0.38268343, 0.39774847, 0.41270703, 0.42755509, 0.44228869, 0.45690388, 0.47139674, 0.48576339,
    0.50000000, 0.51410274, 0.52806785, 0.54189158, 0.55557023, 0.56910015, 0.58247770, 0.59569930,
    0.60876143, 0.62166057, 0.63439328, 0.64695615, 0.65934582, 0.67155895, 0.68359230, 0.69544264,
    0.70710678, 0.71858162, 0.72986407, 0.74095113, 0.75183981, 0.76252720, 0.77301045, 0.78328675,
    0.79335334, 0.80320753, 0.81284668, 0.82226822, 0.83146961, 0.84044840, 0.84920218, 0.85772861,
    0.86602540, 0.87409034, 0.88192126, 0.88951608, 0.89687274, 0.90398929, 0.91086382, 0.91749450,
    0.92387953, 0.93001722, 0.93590593, 0.94154407, 0.94693013, 0.95206268, 0.95694034, 0.96156180,
    0.96592583, 0.97003125, 0.97387698, 0.97746197, 0.98078528, 0.98384601, 0.98664333, 0.98917651,
    0.99144486, 0.99344778, 0.99518473, 0.99665524, 0.99785892, 0.99879546, 0.99946459, 0.99986614,
];

// Applies the filterbank window to a buffer.
fn apply_filter_bank_window(signal: &mut Signal) {
    for (signal, &filter) in signal.iter_mut().zip(BLOCKS_160_W_256_FIRST_HALF.iter()) {
        *signal *= filter;
    }

    const BLOCKS_WITHOUT_ZERO: &[f32; OVERLAP_SIZE - 1] = {
        let [_, rest @ ..] = BLOCKS_160_W_256_FIRST_HALF;
        rest
    };

    for (signal, &filter) in signal.iter_mut().rev().zip(BLOCKS_WITHOUT_ZERO.iter()) {
        *signal *= filter;
    }
}

// Extends a frame with previous data.
fn form_extended_frame(frame: &Frame, old_data: &mut Memory, extended_frame: &mut Signal) {
    let (left, right) = extended_frame.split_at_mut(OVERLAP_SIZE);
    left.copy_from_slice(old_data);
    right.copy_from_slice(frame);
    old_data.copy_from_slice(&extended_frame[NS_FRAME_SIZE..]);
}

// Uses overlap-and-add to produce an output frame.
fn overlap_and_add(extended_frame: &Signal, overlap_memory: &mut Memory, output_frame: &mut Frame) {
    let (left, right) = extended_frame.split_at(OVERLAP_SIZE);
    let (mid, right) = right.split_at(OVERLAP_REMAINDER);
    let (output_left, output_right) = output_frame.split_at_mut(OVERLAP_SIZE);
    output_left
        .iter_mut()
        .zip(overlap_memory.iter())
        .zip(left.iter())
        .for_each(|((out, &overlap), &ext)| *out = overlap + ext);
    output_right.copy_from_slice(mid);
    overlap_memory.copy_from_slice(right);
}

// Produces a delayed frame in place with scratch.
// fn delay_signal(frame: &Frame, delay_buffer: &mut Memory, delayed_frame: &mut Frame) {
//     let (delayed_left, delayed_right) = delayed_frame.split_at_mut(OVERLAP_SIZE);
//     let (frame_left, frame_right) = frame.split_at(OVERLAP_REMAINDER);
//     delayed_left.copy_from_slice(delay_buffer);
//     delayed_right.copy_from_slice(frame_left);
//     delay_buffer.copy_from_slice(frame_right);
// }
fn delay_signal_in_place(frame: &mut Frame, delay_buffer: &mut Memory, scratch: &mut Signal) {
    let scratch = &mut scratch[..OVERLAP_SIZE];
    scratch.copy_from_slice(delay_buffer);
    delay_buffer.copy_from_slice(&frame[OVERLAP_REMAINDER..]);
    let (frame_left, frame_right) = frame.split_at_mut(OVERLAP_SIZE);
    frame_right.copy_from_slice(&frame_left[..OVERLAP_REMAINDER]);
    frame_left.copy_from_slice(scratch);
}

// Computes the energy of an extended frame.
fn compute_energy_of_extended_frame(x: &Signal) -> f32 {
    x.iter().map(|&x| x * x).sum()
}

// Computes the energy of an extended frame based on its subcomponents.
fn compute_energy_of_extended_frame_subcomponents(frame: &Frame, old_data: &Memory) -> f32 {
    old_data.iter().map(|&x| x * x).sum::<f32>() + frame.iter().map(|&x| x * x).sum::<f32>()
}

// Computes the magnitude spectrum based on an FFT output.
fn compute_magnitude_spectrum(fft: &ComplexSpectrum, signal_spectrum: &mut Spectrum) {
    let [first, spectrum @ .., last] = signal_spectrum;
    let [fft_first, fft @ ..] = fft;
    *first = fft_first.re.abs() + 1.;
    *last = fft_first.im.abs() + 1.;
    for (&c, spectrum) in fft.iter().zip(spectrum.iter_mut()) {
        *spectrum = c.re.hypot(c.im) + 1.;
    }
}

fn compute_magnitude_spectrum_in_place(fft: &mut ComplexSpectrum) -> &mut Spectrum {
    let fft = super::cast_complex_to_real(fft);
    fft[0] = fft[0].abs() + 1.;
    fft[FFT_SIZE_BY_2] = fft[1].abs() + 1.;
    for i in 1..FFT_SIZE_BY_2 {
        fft[i] = fft[i * 2].hypot(fft[i * 2 + 1]) + 1.;
    }
    (&mut fft[..FFT_SIZE_BY_2_PLUS_1]).try_into().unwrap()
}

fn compute_signal_energy(fft: &ComplexSpectrum) -> f32 {
    fft.iter().map(|c| c.norm_sqr()).sum::<f32>() * ONE_BY_FFT_SIZE_BY_2_PLUS_1
}

// Compute post SNR.
fn compute_post_snr(signal_spectrum: &Spectrum, noise_spectrum: &Spectrum) -> Spectrum {
    core::array::from_fn(|i| {
        let (signal, noise) = (signal_spectrum[i], noise_spectrum[i]);
        if signal > noise {
            signal / (noise + 0.0001) - 1.
        } else {
            0.
        }
    })
}

// Compute prior SNR.
fn compute_prior_snr(
    filter: &Spectrum,
    prev_signal_spectrum: &Spectrum,
    prev_noise_spectrum: &Spectrum,
    post_snr: &Spectrum,
) -> Spectrum {
    core::array::from_fn(|i| {
        let prev_estimate = prev_signal_spectrum[i] / (prev_noise_spectrum[i] + 0.0001) * filter[i];
        0.98 * prev_estimate + (1. - 0.98) * post_snr[i]
    })
}

// Computes the attenuating gain for the noise suppression of the upper bands.
fn compute_upper_bands_gain(
    minimum_attenuating_gain: f32,
    filter: &Spectrum,
    speech_probability: &Spectrum,
    prev_analysis_signal_spectrum: &Spectrum,
    signal_spectrum: &Spectrum,
) -> f32 {
    // Average speech prob and filter gain for the end of the lowest band.
    const NUM_AVG_BINS: usize = 32;
    const UPPER_BANDS_START: usize = FFT_SIZE_BY_2 - NUM_AVG_BINS;
    const ONE_BY_NUM_AVG_BINS: f32 = 1. / NUM_AVG_BINS as f32;

    let mut avg_prob_speech: f32 =
        speech_probability[UPPER_BANDS_START..FFT_SIZE_BY_2].iter().sum::<f32>();
    let mut avg_filter_gain: f32 = filter[UPPER_BANDS_START..FFT_SIZE_BY_2].iter().sum();
    avg_prob_speech *= ONE_BY_NUM_AVG_BINS;
    avg_filter_gain *= ONE_BY_NUM_AVG_BINS;

    // If the speech was suppressed by a component between Analyze and Process, an
    // example being by an AEC, it should not be considered speech for the purpose
    // of high band suppression. To that end, the speech probability is scaled
    // accordingly.
    let sum_analysis_spectrum: f32 = prev_analysis_signal_spectrum.iter().sum();
    let sum_processing_spectrum: f32 = signal_spectrum.iter().sum();

    // The magnitude spectrum computation enforces the spectrum to be strictly
    // positive.
    avg_prob_speech *= sum_processing_spectrum / sum_analysis_spectrum;

    // Compute gain based on speech probability.
    let gain = 0.5 * (1. + (2. * avg_prob_speech - 1.).tanh());

    // Combine gain with low band gain.
    let gain = if avg_prob_speech >= 0.5 {
        0.25 * gain + 0.75 * avg_filter_gain
    } else {
        0.5 * gain + 0.5 * avg_filter_gain
    };

    // Make sure gain is within flooring range.
    gain.clamp(minimum_attenuating_gain, 1.)
}
