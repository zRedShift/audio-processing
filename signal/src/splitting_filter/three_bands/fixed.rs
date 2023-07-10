use super::{
    DctModulation, Filter, FilterCoefficients, Frame, FullBandFrame, SplitBandFrames, State,
    States, MAX_IN_SHIFT, MEMORY_SIZE, SPLIT_BAND_SIZE, STRIDE, SUBSAMPLING, UPSAMPLING_SCALING,
};

type Q14 = i16;
type Q15 = i16;
type Q17 = i16;
type Q32 = i32;

fn mul_i32(x: i16, y: i16) -> i32 {
    x as i32 * y as i32
}

fn q32_as_q15(x: Q32) -> Q15 {
    (x >> 17) as Q15
}

fn mul_q15_q17(x: Q15, y: Q17) -> Q32 {
    mul_i32(x, y)
}

fn mul_q14_q15(x: Q14, y: Q15) -> Q15 {
    (mul_i32(x, y) >> 15) as Q15
}

fn q14_to_q15(x: Q14) -> Q15 {
    ((x as i32) << 1).clamp(i16::MIN as i32, i16::MAX as i32) as i16
}

// (FILTER_COEFFS[i][j] * (1i32 << 17) as f32).round()
const FILTER_COEFFS: &FilterCoefficients<Q17> = &[
    [-63, -651, 21689, 558],
    [-227, -2079, 19646, 1303],
    [-400, -3324, 15931, 1518],
    [-455, -3392, 6240, 796],
    [-203, -1489, 1819, 244],
    [244, 1819, -1489, -203],
    [796, 6240, -3392, -455],
    [1289, 11198, -3910, -503],
    [1303, 19646, -2079, -227],
    [558, 21689, -651, -63],
];

// (DCT_MODULATION[i][j] * (1i32 << 14) as f32).round() clamped
const DCT_MODULATION: &DctModulation<Q14> = &[
    [32767, 32767, 32767],
    [28378, 0, -28378],
    [16384, -32768, 16384],
    [-16384, 32767, -16384],
    [-28378, 0, 28378],
    [-32768, -32768, -32768],
    [-28378, 0, 28378],
    [-16384, 32767, -16384],
    [16384, -32768, 16384],
    [28378, 0, -28378],
];

// Filters the input signal |in| with the filter |filter| using a shift by
// |in_shift|, taking into account the previous state.
fn filter_core(
    filter: &Filter<Q17>,
    input: &Frame<Q15>,
    in_shift: usize,
    output: &mut Frame<Q15>,
    state: &mut State<Q15>,
) {
    debug_assert!(in_shift <= MAX_IN_SHIFT);
    let (rest, third) = output.split_at_mut(MEMORY_SIZE + 1);
    let (first, second) = rest.split_at_mut(in_shift);
    for (k, out) in first.iter_mut().enumerate() {
        let mut j = MEMORY_SIZE + k - in_shift;
        let mut sum = 0;
        for &f in filter {
            sum += mul_q15_q17(state[j], f);
            j = j.wrapping_sub(STRIDE);
        }
        *out = q32_as_q15(sum);
    }

    let mut shift = 0;
    for out in second {
        let mut j = shift;
        let mut sum = 0;
        let loop_limit = super::loop_limit(shift);
        let (pre, post) = filter.split_at(loop_limit);
        for &f in pre {
            sum += mul_q15_q17(input[j], f);
            j = j.wrapping_sub(STRIDE);
        }
        j = MEMORY_SIZE + shift - loop_limit * STRIDE;
        for &f in post {
            sum += mul_q15_q17(state[j], f);
            j = j.wrapping_sub(STRIDE);
        }
        *out = q32_as_q15(sum);
        shift += 1;
    }

    for out in third {
        let mut j = shift;
        let mut sum = 0;
        for &f in filter {
            sum += mul_q15_q17(input[j], f);
            j = j.wrapping_sub(STRIDE);
        }
        *out = q32_as_q15(sum);
        shift += 1;
    }

    state.copy_from_slice(&input[SPLIT_BAND_SIZE - MEMORY_SIZE..]);
}

impl super::Num for Q15 {
    const ZERO: Self = 0;

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
                let mut output_subsampled = [0; SPLIT_BAND_SIZE];
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
                        *output = output.wrapping_add(mul_q14_q15(modulation, output_subsampled));
                    }
                }
            }
        }
        for output in output.iter_mut() {
            for output in output {
                *output = q14_to_q15(*output);
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
                let mut input_subsampled = [0i16; SPLIT_BAND_SIZE];
                for (input, &modulation) in input.iter().zip(&DCT_MODULATION[filter_index]) {
                    for (input_subsampled, &input) in input_subsampled.iter_mut().zip(input) {
                        *input_subsampled =
                            input_subsampled.wrapping_add(mul_q14_q15(modulation, input));
                    }
                }

                // Filter.
                let mut output_subsampled = [0; SPLIT_BAND_SIZE];
                filter_core(
                    &FILTER_COEFFS[filter_index],
                    &input_subsampled,
                    in_shift,
                    &mut output_subsampled,
                    &mut state[filter_index],
                );

                // Upsample.
                for (k, output_subsampled) in output_subsampled.into_iter().enumerate() {
                    const UPSAMPLING: i16 = UPSAMPLING_SCALING as i16;
                    let output = &mut output[upsampling_index + SUBSAMPLING * k];
                    *output = output.wrapping_add(UPSAMPLING * output_subsampled);
                }
            }
        }
        for output in output.iter_mut() {
            *output = q14_to_q15(*output);
        }
    }
}
