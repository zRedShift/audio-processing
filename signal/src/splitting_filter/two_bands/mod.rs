const STATE_SIZE: usize = 6;
const FILTER_SIZE: usize = 3;
const NUM_BANDS: usize = 2;
const FULL_BAND_SIZE: usize = 320;
const SPLIT_BAND_SIZE: usize = FULL_BAND_SIZE / NUM_BANDS;
const NUM_FILTERS: usize = 2;

use static_assertions::const_assert_eq;
const_assert_eq!(NUM_BANDS * SPLIT_BAND_SIZE, FULL_BAND_SIZE);

type Filter<T> = [T; FILTER_SIZE];
type Frame<T> = [T; SPLIT_BAND_SIZE];
type State<T> = [T; STATE_SIZE];
type States<T> = [State<T>; NUM_FILTERS];
type FilterCoefficients<T> = [Filter<T>; NUM_FILTERS];
type FullBandFrame<T> = [T; FULL_BAND_SIZE];
type SplitBandFrames<T> = [Frame<T>; NUM_BANDS];
type FullBandChunkedFrame<T> = [[T; NUM_BANDS]; SPLIT_BAND_SIZE];

mod fixed;
mod float;

pub struct TwoBandFilterBank<T: Num> {
    state_analysis: States<T::StateType>,
    state_synthesis: States<T::StateType>,
}

fn cast_to_chunked<T>(full_band: &FullBandFrame<T>) -> &FullBandChunkedFrame<T> {
    // static_assertions::assert_eq_size!(FullBandFrame<T>, FullBandChunkedFrame<T>);
    // static_assertions::assert_eq_align!(FullBandFrame<T>, FullBandChunkedFrame<T>);
    unsafe { &*(full_band as *const _ as *const _) }
}

fn cast_to_chunked_mut<T>(full_band: &mut FullBandFrame<T>) -> &mut FullBandChunkedFrame<T> {
    // static_assertions::assert_eq_size!(FullBandFrame<T>, FullBandChunkedFrame<T>);
    // static_assertions::assert_eq_align!(FullBandFrame<T>, FullBandChunkedFrame<T>);
    unsafe { &mut *(full_band as *mut _ as *mut _) }
}

pub trait Num: Copy {
    type StateType: Copy;
    const ZERO: Self::StateType;

    fn analysis(
        state: &mut States<Self::StateType>,
        input: &FullBandFrame<Self>,
        output: &mut SplitBandFrames<Self>,
    );

    fn synthesis(
        state: &mut States<Self::StateType>,
        input: &SplitBandFrames<Self>,
        output: &mut FullBandFrame<Self>,
    );
}

impl<T: Num> TwoBandFilterBank<T> {
    pub const fn new() -> Self {
        Self {
            state_analysis: [[T::ZERO; STATE_SIZE]; NUM_FILTERS],
            state_synthesis: [[T::ZERO; STATE_SIZE]; NUM_FILTERS],
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

    // Splits a 0-2*F Hz signal into two sub bands: 0-F Hz and F-2*F Hz. The
    // current version has F = 8000, therefore, a super-wideband audio signal is
    // split to lower-band 0-8 kHz and upper-band 8-16 kHz.
    pub fn analysis(&mut self, input: &FullBandFrame<T>, output: &mut SplitBandFrames<T>) {
        T::analysis(&mut self.state_analysis, input, output)
    }

    // Combines the two sub bands (0-F and F-2*F Hz) into a signal of 0-2*F
    // Hz, (current version has F = 8000 Hz). So the filter combines lower-band
    // (0-8 kHz) and upper-band (8-16 kHz) channels to obtain super-wideband 0-16
    // kHz audio.
    pub fn synthesis(&mut self, input: &SplitBandFrames<T>, output: &mut FullBandFrame<T>) {
        T::synthesis(&mut self.state_synthesis, input, output)
    }
}
