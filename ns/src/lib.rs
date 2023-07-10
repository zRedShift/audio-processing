mod config;
mod fft;
mod histograms;
mod noise_estimator;
mod noise_suppressor;
mod prior_signal_model;
mod prior_signal_model_estimator;
mod quantile_noise_estimator;
mod signal_model;
mod signal_model_estimator;
mod speech_probability_estimator;
mod suppression_params;
mod wiener_filter;

pub use config::{SampleRate, SuppressionLevel};
use fft::{cast_complex_to_real, cast_real_to_complex};
pub use fft::{ComplexDftAdapter, Fft, Fft4g};
use histograms::{Histogram, Histograms};
use noise_estimator::NoiseEstimator;
pub use noise_suppressor::NoiseSuppressor;
use prior_signal_model::PriorSignalModel;
use prior_signal_model_estimator::PriorSignalModelEstimator;
use quantile_noise_estimator::QuantileNoiseEstimator;
use signal_model::SignalModel;
use signal_model_estimator::SignalModelEstimator;
use speech_probability_estimator::SpeechProbabilityEstimator;
use suppression_params::SuppressionParams;
use wiener_filter::WienerFilter;

const FFT_SIZE: usize = 256;
const FFT_SIZE_BY_2: usize = FFT_SIZE / 2;
const FFT_SIZE_BY_2_PLUS_1: usize = FFT_SIZE_BY_2 + 1;
const NS_FRAME_SIZE: usize = 160;
const OVERLAP_SIZE: usize = FFT_SIZE - NS_FRAME_SIZE;

const SHORT_STARTUP_PHASE_BLOCKS: i32 = 50;
const LONG_STARTUP_PHASE_BLOCKS: i32 = 200;
const FEATURE_UPDATE_WINDOW_SIZE: i32 = 500;

const LRT_FEATURE_THR: f32 = 0.5;
const BIN_SIZE_LRT: f32 = 0.1;
const BIN_SIZE_SPEC_FLAT: f32 = 0.05;
const BIN_SIZE_SPEC_DIFF: f32 = 0.1;

const ONE_BY_FFT_SIZE_BY_2_PLUS_1: f32 = 1. / FFT_SIZE_BY_2_PLUS_1 as f32;
const ONE_BY_SHORT_STARTUP_PHASE_BLOCKS: f32 = 1. / SHORT_STARTUP_PHASE_BLOCKS as f32;

pub type Frame = [f32; NS_FRAME_SIZE];
pub type SplitAudioFrames<const BANDS: usize = 1, const CH: usize = 1> = [[Frame; BANDS]; CH];
type Signal = [f32; FFT_SIZE];
type Spectrum = [f32; FFT_SIZE_BY_2_PLUS_1];
type ComplexSpectrum = [num_complex::Complex32; FFT_SIZE_BY_2];
