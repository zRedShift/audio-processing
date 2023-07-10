const SPARSITY: usize = 4;
const STRIDE_LOG2: usize = 2;
const STRIDE: usize = 1 << STRIDE_LOG2;
const NUM_ZERO_FILTERS: usize = 2;
const FILTER_SIZE: usize = 4;
const MEMORY_SIZE: usize = FILTER_SIZE * STRIDE - 1;
const NUM_BANDS: usize = 3;
const FULL_BAND_SIZE: usize = 480;
const SPLIT_BAND_SIZE: usize = FULL_BAND_SIZE / NUM_BANDS;
const NUM_NON_ZERO_FILTERS: usize = SPARSITY * NUM_BANDS - NUM_ZERO_FILTERS;
const SUBSAMPLING: usize = NUM_BANDS;
const DCT_SIZE: usize = NUM_BANDS;
const MAX_IN_SHIFT: usize = STRIDE - 1;
const UPSAMPLING_SCALING: usize = SUBSAMPLING;

use static_assertions::const_assert_eq;
const_assert_eq!(MEMORY_SIZE, 15);
const_assert_eq!(NUM_BANDS * SPLIT_BAND_SIZE, FULL_BAND_SIZE);

type Filter<T> = [T; FILTER_SIZE];
type Frame<T> = [T; SPLIT_BAND_SIZE];
type State<T> = [T; MEMORY_SIZE];
type States<T> = [State<T>; NUM_NON_ZERO_FILTERS];
type FilterCoefficients<T> = [Filter<T>; NUM_NON_ZERO_FILTERS];
type DctModulation<T> = [[T; DCT_SIZE]; NUM_NON_ZERO_FILTERS];
type FullBandFrame<T> = [T; FULL_BAND_SIZE];
type SplitBandFrames<T> = [Frame<T>; NUM_BANDS];

mod fixed;
mod float;

// Factors to take into account when choosing |kFilterSize|:
//   1. Higher |kFilterSize|, means faster transition, which ensures less
//      aliasing. This is especially important when there is non-linear
//      processing between the splitting and merging.
//   2. The delay that this filter bank introduces is
//      |kNumBands| * |kSparsity| * |kFilterSize| / 2, so it increases linearly
//      with |kFilterSize|.
//   3. The computation complexity also increases linearly with |kFilterSize|.

// The Matlab code to generate these |kFilterCoeffs| is:
//
// N = kNumBands * kSparsity * kFilterSize - 1;
// h = fir1(N, 1 / (2 * kNumBands), kaiser(N + 1, 3.5));
// reshape(h, kNumBands * kSparsity, kFilterSize);
//
// The code below uses the values of kFilterSize, kNumBands and kSparsity
// specified in the header.

// Because the total bandwidth of the lower and higher band is double the middle
// one (because of the spectrum parity), the low-pass prototype is half the
// bandwidth of 1 / (2 * |kNumBands|) and is then shifted with cosine modulation
// to the right places.
// A Kaiser window is used because of its flexibility and the alpha is set to
// 3.5, since that sets a stop band attenuation of 40dB ensuring a fast
// transition.

// An implementation of a 3-band FIR filter-bank with DCT modulation, similar to
// the proposed in "Multirate Signal Processing for Communication Systems" by
// Fredric J Harris.
// The low-pass filter prototype has these characteristics:
// * Pass-band ripple = 0.3dB
// * Pass-band frequency = 0.147 (7kHz at 48kHz)
// * Stop-band attenuation = 40dB
// * Stop-band frequency = 0.192 (9.2kHz at 48kHz)
// * Delay = 24 samples (500us at 48kHz)
// * Linear phase
// This filter bank does not satisfy perfect reconstruction. The SNR after
// analysis and synthesis (with no processing in between) is approximately 9.5dB
// depending on the input signal after compensating for the delay.
pub struct ThreeBandFilterBank<T> {
    state_analysis: States<T>,
    state_synthesis: States<T>,
}

fn loop_limit(shift: usize) -> usize {
    FILTER_SIZE.min(1 + (shift >> STRIDE_LOG2))
}

fn filter_index(index: usize) -> Option<usize> {
    match index {
        0..=2 => Some(index),
        4..=8 => Some(index - 1),
        10..=11 => Some(index - 2),
        _ => None,
    }
}

pub trait Num: Copy {
    const ZERO: Self;

    fn analysis(
        state: &mut States<Self>,
        input: &FullBandFrame<Self>,
        output: &mut SplitBandFrames<Self>,
    );

    fn synthesis(
        state: &mut States<Self>,
        input: &SplitBandFrames<Self>,
        output: &mut FullBandFrame<Self>,
    );
}

impl<T: Num> ThreeBandFilterBank<T> {
    // Because the low-pass filter prototype has half bandwidth it is possible to
    // use a DCT to shift it in both directions at the same time, to the center
    // frequencies [1 / 12, 3 / 12, 5 / 12].
    pub const fn new() -> Self {
        Self {
            state_analysis: [[T::ZERO; MEMORY_SIZE]; NUM_NON_ZERO_FILTERS],
            state_synthesis: [[T::ZERO; MEMORY_SIZE]; NUM_NON_ZERO_FILTERS],
        }
    }

    pub fn reset(&mut self) {
        for state in &mut self.state_analysis {
            state.fill(T::ZERO);
        }
        for state in &mut self.state_synthesis {
            state.fill(T::ZERO);
        }
    }

    // The analysis can be separated in these steps:
    //   1. Serial to parallel downsampling by a factor of |kNumBands|.
    //   2. Filtering of |kSparsity| different delayed signals with polyphase
    //      decomposition of the low-pass prototype filter and upsampled by a factor
    //      of |kSparsity|.
    //   3. Modulating with cosines and accumulating to get the desired band.
    pub fn analysis(&mut self, input: &FullBandFrame<T>, output: &mut SplitBandFrames<T>) {
        for state in output.iter_mut() {
            state.fill(T::ZERO);
        }
        T::analysis(&mut self.state_analysis, input, output)
    }

    // The synthesis can be separated in these steps:
    //   1. Modulating with cosines.
    //   2. Filtering each one with a polyphase decomposition of the low-pass
    //      prototype filter upsampled by a factor of |kSparsity| and accumulating
    //      |kSparsity| signals with different delays.
    //   3. Parallel to serial upsampling by a factor of |kNumBands|.
    pub fn synthesis(&mut self, input: &SplitBandFrames<T>, output: &mut FullBandFrame<T>) {
        output.fill(T::ZERO);
        T::synthesis(&mut self.state_synthesis, input, output)
    }
}
