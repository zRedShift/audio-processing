use ns::{ComplexDftAdapter, NoiseSuppressor, SuppressionLevel};
use rustfft::FftPlanner;
use signal::ThreeBandFilterBank;
use std::fs::File;
use std::io::{Read, Write};

fn main() {
    let mut f = File::open("~/Downloads/car_10dB_e.pcm").unwrap();
    let mut arr = vec![0i16; f.metadata().unwrap().len() as usize / 2];
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(128);
    let mut scratch = [Default::default(); 128];
    let fft = ComplexDftAdapter::new(move |s| fft.process_with_scratch(s, &mut scratch));
    let bytes = bytemuck::cast_slice_mut(&mut arr);
    f.read_exact(bytes).unwrap();
    let arr_chunked: &mut [[i16; 480]] = bytemuck::cast_slice_mut(&mut arr);
    let split_band_frames = &mut [[[0; 160]; 3]];
    let mut filter = ThreeBandFilterBank::new();
    let now = std::time::Instant::now();
    for frame in arr_chunked.iter_mut() {
        filter.analysis(frame, &mut split_band_frames[0]);
        filter.synthesis(&split_band_frames[0], frame);
    }
    // let mut v: Vec<f32> = arr.iter().map(|&x| x as _).collect();
    // let v_chunked: &mut [[f32; 480]] = bytemuck::cast_slice_mut(&mut v);
    // let mut suppressor = NoiseSuppressor::<_, 2>::new(fft, SuppressionLevel::Level21dB);
    // let now = std::time::Instant::now();
    // let mut extended_frames = [[0.; 256]; 1];
    // let mut signal_spectrum = [0.; 129];
    // let mut filter = ThreeBandFilterBank::new();
    // let split_band_frames = &mut [[[0.; 160]; 3]];
    // for frame in v_chunked.iter_mut() {
    //     filter.analysis(frame, &mut split_band_frames[0]);
    //     // suppressor.analyze_with_scratch(split_band_frames, &mut extended_frames[0]);
    //     // suppressor.process_with_scratch(
    //     //     split_band_frames,
    //     //     &mut extended_frames,
    //     //     &mut signal_spectrum,
    //     // );
    //     filter.synthesis(&split_band_frames[0], frame);
    // }
    let elapsed = now.elapsed();
    println!("{} {}", elapsed.as_nanos(), arr_chunked.len());
    // v.iter().zip(arr.iter_mut()).for_each(|(&v, a)| *a = v as _);
    File::create("~/Downloads/car_10dB_48000_resynthed_s16.pcm")
        .unwrap()
        .write_all(bytemuck::cast_slice(&arr))
        .unwrap();
    // let mut planner = FftPlanner::new();
    // let mut fft = NsFft::new();
    // let mut rfft = Fft::new(&mut planner);
    // let mut orig_signal = [0.; 256];
    // for (i, s) in orig_signal.iter_mut().enumerate() {
    //     *s = (i as f32 * core::f32::consts::PI / 64.).cos();
    //     *s += (i as f32 * core::f32::consts::PI / 128.).cos();
    //     *s += (i as f32 * core::f32::consts::PI / 256.).cos();
    //     *s += (i as f32 * core::f32::consts::PI / 32.).cos();
    //     *s += (i as f32 * core::f32::consts::PI / 16.).cos();
    //     *s += (i as f32 * core::f32::consts::PI / 8.).cos();
    //     *s += (i as f32 * core::f32::consts::PI / 4.).cos();
    //     *s += (i as f32 * core::f32::consts::PI / 2.).cos();
    //     *s /= 8.;
    // }
    // let mut signal = orig_signal;
    // println!("{signal:?}");
    // let spectrum = fft.fft(&mut signal);
    // println!("{spectrum:?}");
    // fft.ifft(spectrum);
    // println!("{signal:?}");
    // let count = signal
    //     .iter()
    //     .zip(orig_signal.iter())
    //     .map(|(&s, &o)| (s - o) * (s - o))
    //     .sum::<f32>();
    // println!("fft4g: {}", (count / 256.).sqrt());
    // let mut signal = orig_signal;
    // let spectrum = rfft.real_dft(&mut signal);
    // println!("\n\n\n\n{spectrum:?}");
    // rfft.inverse_real_dft(spectrum);
    // println!("{signal:?}");
    // let count = signal
    //     .iter()
    //     .zip(orig_signal.iter())
    //     .map(|(&s, &o)| (s - o) * (s - o))
    //     .sum::<f32>();
    // println!("rustfft: {}", (count / 256.).sqrt());
}

/*
use audio_processing::NsFft;
use realfft::{num_complex::Complex32 as Complex, RealToComplex, RealToComplexEven, ComplexToRealEven, ComplexToReal};
use rustfft::FftPlanner;

fn main() {
    let mut fft = NsFft::new();
    let mut orig_signal = [0.; 256];
    let mut real = [0.; 129];
    let mut imag = [0.; 129];
    for (i, s) in orig_signal.iter_mut().enumerate() {
        *s = (i as f32 * core::f32::consts::PI / 64.).cos();
        *s += (i as f32 * core::f32::consts::PI / 128.).cos();
        *s += (i as f32 * core::f32::consts::PI / 256.).cos();
        *s += (i as f32 * core::f32::consts::PI / 32.).cos();
        *s += (i as f32 * core::f32::consts::PI / 16.).cos();
        *s += (i as f32 * core::f32::consts::PI / 8.).cos();
        *s += (i as f32 * core::f32::consts::PI / 4.).cos();
        *s += (i as f32 * core::f32::consts::PI / 2.).cos();
        *s /= 8.;
    }
    println!("{orig_signal:?}");
    let mut signal = orig_signal;
    fft.fft(&mut signal, &mut real, &mut imag);

    println!("\n\n{signal:?}");
    // println!("{real:?}");
    // println!("{imag:?}");
    // println!("{fft:?}");

    let mut plan = FftPlanner::<f32>::new();
    let fft__ = plan.plan_fft_forward(128);
    let mut signal_: [Complex; 128] = unsafe { core::mem::transmute_copy(&orig_signal) };
    let mut scratch = [Complex::new(0., 0.); 128];
    fft__.process_outofplace_with_scratch(&mut signal_, &mut scratch, &mut []);
    signal = unsafe { core::mem::transmute_copy(&scratch) };
    println!("\n{signal:?}");

    signal = orig_signal;
    let fft_ = RealToComplexEven::new(256, &mut plan);
    let ifft_ = ComplexToRealEven::new(256, &mut plan);
    fn compute_twiddle(index: usize) -> Complex {
        let constant = -2f64 * core::f64::consts::PI / 256f64;
        let angle = constant * index as f64;
        Complex {
            re: angle.cos() as _,
            im: angle.sin() as _,
        }
    }
    let mut twiddles = [0.; 126];
    for i in 0..63 {
        let c = compute_twiddle(i + 1) * 0.5;
        twiddles[i * 2] = c.re;
        twiddles[i * 2 + 1] = c.im;
    }
    println!("\n\n\n{twiddles:?}");
    let twiddles: Vec<Complex> = (1..64).map(|i| compute_twiddle(i) * 0.5).collect();
    println!("\n\n\n{twiddles:?}");
    let mut fft_res = core::array::from_fn::<_, 129, _>(|i| Complex::new(real[i], -imag[i]));
    let or = fft_res;
    println!("\n{fft_res:?}");
    fft_.process_with_scratch(&mut signal, &mut fft_res, &mut [])
        .unwrap();

    // println!("{signal:?}");
    println!("\n{fft_res:?}");
    ifft_.process_with_scratch(&mut fft_res, &mut signal, &mut []).unwrap();
    println!("\n\n\n\n{signal:?}");
    let mut arr = [0.; 129];
    let mut arr2 = [0.; 129];
    let mut arr3 = [0.; 129];
    let mut arr4 = [Complex::new(0., 0.); 129];
    let mut count = 0.;
    for i in 0..129 {
        let diff = fft_res[i] - or[i];
        arr4[i] = diff;
        count += diff.norm_sqr();
        arr[i] = fft_res[i].norm();
        arr2[i] = or[i].norm();
        arr3[i] = diff.norm();
    }
    println!("{}", (count / 129.).sqrt());
    // println!("{}", fft_res[1]);
    // println!("{}", or[1]);
    println!("\n{arr:?}");
    println!("\n{arr2:?}");
    println!("\n{arr3:?}");
    println!("\n{arr4:?}");
}

 */
