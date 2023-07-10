use super::{
    Filter, FilterCoefficients, Frame, FullBandFrame, SplitBandFrames, State, States, NUM_BANDS,
    SPLIT_BAND_SIZE,
};

type Q10 = i32;
type Q15 = i16;
type Q16 = u16;

const Q10_MULT: Q10 = 1 << 10;

// c + the 32 most significant bits of a * b
fn scale_diff(a: Q16, b: Q10, c: Q10) -> Q10 {
    c + (b >> 16) * a as Q10 + (((b as u32 & 0x0000FFFF) * a as u32) >> 16) as Q10
}

fn q10_to_q15<const SHIFT: u8, const ADD: Q10>(x: Q10) -> Q15 {
    ((x + ADD) >> SHIFT).clamp(i16::MIN as i32, i16::MAX as i32) as i16
}

// QMF filter coefficients in Q16.
const FILTER_COEFFS: &FilterCoefficients<Q16> = &[[6418, 36982, 57261], [21333, 49062, 63010]];

// Allpass filter used by the analysis and synthesis parts of the QMF filter.
fn filter_core(
    in_data: &mut Frame<Q10>,
    &filter_coefficients: &Filter<Q16>,
    filter_state: &mut State<Q10>,
) -> Frame<Q10> {
    // The procedure is to filter the input with three first order all pass filters
    // (cascade operations).
    //
    //         a_3 + q^-1    a_2 + q^-1    a_1 + q^-1
    // y[n] =  -----------   -----------   -----------   x[n]
    //         1 + a_3q^-1   1 + a_2q^-1   1 + a_1q^-1
    //
    // The input vector |filter_coefficients| includes these three filter coefficients.
    // The filter state contains the in_data state, in_data[-1], followed by
    // the out_data state, out_data[-1]. This is repeated for each cascade.
    // The first cascade filter will filter the |in_data| and store the output in
    // |out_data|. The second will the take the |out_data| as input and make an
    // intermediate storage in |in_data|, to save memory. The third, and final, cascade
    // filter operation takes the |in_data| (which is the output from the previous cascade
    // filter) and store the output in |out_data|.
    // Note that the input vector values are changed during the process.

    // First all-pass cascade; filter from in_data to out_data.

    // Let y_i[n] indicate the output of cascade filter i (with filter coefficient a_i) at
    // vector position n. Then the final output will be y[n] = y_3[n]

    // First loop, use the states stored in memory.
    // "diff" should be safe from wrap around since max values are 2^25
    let mut prev_in = filter_state[0];
    let mut prev_out = filter_state[1];

    let mut out_data: Frame<_> = core::array::from_fn(|k| {
        let curr_in = in_data[k];
        // diff = (x[n] - y_1[n-1])
        let diff = curr_in.saturating_sub(prev_out);
        // y_1[n] =  x[n-1] + a_1 * (x[n] - y_1[n-1])
        prev_out = scale_diff(filter_coefficients[0], diff, prev_in);
        prev_in = curr_in;
        prev_out
    });

    // Update states.
    filter_state[0] = prev_in; // x[N-1], becomes x[-1] next time
    filter_state[1] = prev_out; // y_1[N-1], becomes y_1[-1] next time

    // Second all-pass cascade; filter from out_data to in_data.
    prev_in = filter_state[3];
    prev_out = filter_state[2];

    for (curr_in, &curr_out) in in_data.iter_mut().zip(&out_data) {
        // diff = (y_1[n] - y_2[n-1])
        let diff = curr_out.saturating_sub(prev_in);
        // y_2[0] =  y_1[-1] + a_2 * (y_1[0] - y_2[-1])
        prev_in = scale_diff(filter_coefficients[1], diff, prev_out);
        prev_out = curr_out;
        *curr_in = prev_in;
    }

    filter_state[2] = prev_out; // y_1[N-1], becomes y_1[-1] next time
    filter_state[3] = prev_in; // y_2[N-1], becomes y_2[-1] next time

    // Third all-pass cascade; filter from in_data to out_data.
    prev_in = filter_state[4];
    prev_out = filter_state[5];

    for (curr_out, &curr_in) in out_data.iter_mut().zip(&*in_data) {
        // diff = (y_2[n] - y[n-1])
        let diff = curr_in.saturating_sub(prev_out);
        // y[n] =  y_2[n-1] + a_3 * (y_2[n] - y[n-1])
        prev_out = scale_diff(filter_coefficients[2], diff, prev_in);
        prev_in = curr_in;
        *curr_out = prev_out;
    }

    filter_state[4] = prev_in; // y_2[N-1], becomes y_2[-1] next time
    filter_state[5] = prev_out; // y[N-1], becomes y[-1] next time

    out_data
}

impl super::Num for Q15 {
    type StateType = Q10;
    const ZERO: Self::StateType = 0;

    fn analysis(
        state: &mut States<Self::StateType>,
        input: &FullBandFrame<Self>,
        [low_band, high_band]: &mut SplitBandFrames<Self>,
    ) {
        let [mut half_in1, mut half_in2] = [[0; SPLIT_BAND_SIZE]; NUM_BANDS];

        // Split even and odd samples. Also shift them to Q10.
        for ((half_in1, half_in2), &[even, odd]) in
            half_in1.iter_mut().zip(&mut half_in2).zip(super::cast_to_chunked(input))
        {
            *half_in2 = even as Q10 * Q10_MULT;
            *half_in1 = odd as Q10 * Q10_MULT;
        }

        // All pass filter even and odd samples, independently.
        let filter1 = filter_core(&mut half_in1, &FILTER_COEFFS[0], &mut state[0]);
        let filter2 = filter_core(&mut half_in2, &FILTER_COEFFS[1], &mut state[1]);

        // Take the sum and difference of filtered version of odd and even
        // branches to get upper & lower band.
        for ((low_band, filter1), (high_band, filter2)) in
            low_band.iter_mut().zip(filter1).zip(high_band.iter_mut().zip(filter2))
        {
            const SHIFT: u8 = 11;
            const ADD: Q10 = 1 << 10;
            *low_band = q10_to_q15::<SHIFT, ADD>(filter1 + filter2);
            *high_band = q10_to_q15::<SHIFT, ADD>(filter1 - filter2);
        }
    }

    fn synthesis(
        state: &mut States<Self::StateType>,
        [low_band, high_band]: &SplitBandFrames<Self>,
        output: &mut FullBandFrame<Self>,
    ) {
        let [mut half_in1, mut half_in2] = [[0; SPLIT_BAND_SIZE]; NUM_BANDS];

        // Obtain the sum and difference channels out of upper and lower-band channels.
        // Also shift to Q10 domain.
        for ((&low_band, half_in1), (&high_band, half_in2)) in
            low_band.iter().zip(&mut half_in1).zip(high_band.iter().zip(&mut half_in2))
        {
            *half_in1 = (low_band as Q10 + high_band as Q10) * Q10_MULT;
            *half_in2 = (low_band as Q10 - high_band as Q10) * Q10_MULT;
        }

        // all-pass filter the sum and difference channels
        let filter1 = filter_core(&mut half_in1, &FILTER_COEFFS[1], &mut state[0]);
        let filter2 = filter_core(&mut half_in2, &FILTER_COEFFS[0], &mut state[1]);

        // The filtered signals are even and odd samples of the output. Combine
        // them. The signals are Q10 should shift them back to Q0 and take care of
        // saturation.
        for ((filter1, filter2), [even, odd]) in
            filter1.into_iter().zip(filter2).zip(super::cast_to_chunked_mut(output))
        {
            const SHIFT: u8 = 10;
            const ADD: Q10 = 1 << 9;
            *even = q10_to_q15::<SHIFT, ADD>(filter2);
            *odd = q10_to_q15::<SHIFT, ADD>(filter1);
        }
    }
}
