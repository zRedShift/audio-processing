use super::{
    DctModulation, Filter, FilterCoefficients, Frame, FullBandFrame, SplitBandFrames, State,
    States, MAX_IN_SHIFT, MEMORY_SIZE, SPLIT_BAND_SIZE, STRIDE, SUBSAMPLING, UPSAMPLING_SCALING,
};

const FILTER_COEFFS: &FilterCoefficients<f32> = &[
    [-0.00047749, -0.00496888, 0.16547118, 0.00425496],
    [-0.00173287, -0.01585778, 0.14989004, 0.00994113],
    [-0.00304815, -0.02536082, 0.12154542, 0.01157993],
    [-0.00346946, -0.02587886, 0.04760441, 0.00607594],
    [-0.00154717, -0.01136076, 0.01387458, 0.00186353],
    [0.00186353, 0.01387458, -0.01136076, -0.00154717],
    [0.00607594, 0.04760441, -0.02587886, -0.00346946],
    [0.00983212, 0.08543175, -0.02982767, -0.00383509],
    [0.00994113, 0.14989004, -0.01585778, -0.00173287],
    [0.00425496, 0.16547118, -0.00496888, -0.00047749],
];

#[allow(clippy::excessive_precision)]
const DCT_MODULATION: &DctModulation<f32> = &[
    [2., 2., 2.],
    [1.73205077, 0., -1.73205077],
    [1., -2., 1.],
    [-1., 2., -1.],
    [-1.73205077, 0., 1.73205077],
    [-2., -2., -2.],
    [-1.73205077, 0., 1.73205077],
    [-1., 2., -1.],
    [1., -2., 1.],
    [1.73205077, 0., -1.73205077],
];

// Filters the input signal |in| with the filter |filter| using a shift by
// |in_shift|, taking into account the previous state.
fn filter_core(
    filter: &Filter<f32>,
    input: &Frame<f32>,
    in_shift: usize,
    output: &mut Frame<f32>,
    state: &mut State<f32>,
) {
    debug_assert!(in_shift <= MAX_IN_SHIFT);
    let (rest, third) = output.split_at_mut(MEMORY_SIZE + 1);
    let (first, second) = rest.split_at_mut(in_shift);
    for (k, out) in first.iter_mut().enumerate() {
        let mut j = MEMORY_SIZE + k - in_shift;
        for &f in filter {
            *out += state[j] * f;
            j = j.wrapping_sub(STRIDE);
        }
    }

    let mut shift = 0;
    for out in second {
        let mut j = shift;
        let loop_limit = super::loop_limit(shift);
        let (pre, post) = filter.split_at(loop_limit);
        for &f in pre {
            *out += input[j] * f;
            j = j.wrapping_sub(STRIDE);
        }
        j = MEMORY_SIZE + shift - loop_limit * STRIDE;
        for &f in post {
            *out += state[j] * f;
            j = j.wrapping_sub(STRIDE);
        }
        shift += 1;
    }

    for out in third {
        let mut j = shift;
        for &f in filter {
            *out += input[j] * f;
            j = j.wrapping_sub(STRIDE);
        }
        shift += 1;
    }

    state.copy_from_slice(&input[SPLIT_BAND_SIZE - MEMORY_SIZE..]);
}

impl super::Num for f32 {
    const ZERO: Self = 0.;

    fn analysis(
        state: &mut States<Self>,
        input: &FullBandFrame<Self>,
        output: &mut SplitBandFrames<Self>,
    ) {
        for downsampling_index in 0..SUBSAMPLING {
            // Downsample to form the filter input.
            let input_subsampled: Frame<_> = core::array::from_fn(|k| {
                input[(SUBSAMPLING - 1) - downsampling_index + SUBSAMPLING * k]
            });

            for in_shift in 0..STRIDE {
                // Choose filter, skip zero filters.
                let Some(filter_index) = super::filter_index(downsampling_index + in_shift * SUBSAMPLING) else {
                    continue;
                };

                // Filter.
                let mut output_subsampled = [0.; SPLIT_BAND_SIZE];
                filter_core(
                    &FILTER_COEFFS[filter_index],
                    &input_subsampled,
                    in_shift,
                    &mut output_subsampled,
                    &mut state[filter_index],
                );

                // Band and modulate the output.
                for (output, &modulation) in output.iter_mut().zip(&DCT_MODULATION[filter_index]) {
                    for (output, &output_subsampled) in output.iter_mut().zip(&output_subsampled) {
                        *output += modulation * output_subsampled;
                    }
                }
            }
        }
    }

    fn synthesis(
        state: &mut States<Self>,
        input: &SplitBandFrames<Self>,
        output: &mut FullBandFrame<Self>,
    ) {
        for upsampling_index in 0..SUBSAMPLING {
            for in_shift in 0..STRIDE {
                // Choose filter, skip zero filters.
                let Some(filter_index) = super::filter_index(upsampling_index + in_shift * SUBSAMPLING) else {
                    continue;
                };

                // Prepare filter input by modulating the banded input.
                let mut input_subsampled = [0.; SPLIT_BAND_SIZE];
                for (input, &modulation) in input.iter().zip(&DCT_MODULATION[filter_index]) {
                    for (input_subsampled, &input) in input_subsampled.iter_mut().zip(input) {
                        *input_subsampled += modulation * input;
                    }
                }

                // Filter.
                let mut output_subsampled = [0.; SPLIT_BAND_SIZE];
                filter_core(
                    &FILTER_COEFFS[filter_index],
                    &input_subsampled,
                    in_shift,
                    &mut output_subsampled,
                    &mut state[filter_index],
                );

                // Upsample.
                for (k, output_subsampled) in output_subsampled.into_iter().enumerate() {
                    const UPSAMPLING: f32 = UPSAMPLING_SCALING as f32;
                    output[upsampling_index + SUBSAMPLING * k] += UPSAMPLING * output_subsampled;
                }
            }
        }
    }
}
