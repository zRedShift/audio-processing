use crate::{ComplexSpectrum, Signal, FFT_SIZE, FFT_SIZE_BY_2};
use num_complex::Complex32 as C;

pub struct ComplexDftAdapter<F>(F);

impl<F: FnMut(&mut ComplexSpectrum)> ComplexDftAdapter<F> {
    pub fn new(f: F) -> Self {
        Self(f)
    }
}

impl<F: FnMut(&mut ComplexSpectrum)> super::Fft for ComplexDftAdapter<F> {
    fn real_dft<'a>(&mut self, signal: &'a mut Signal) -> &'a mut ComplexSpectrum {
        let spectrum = super::cast_real_to_complex(signal);
        self.0(spectrum);
        let [left, right] = split_spectrum(spectrum);
        let [first, left @ ..] = left;
        let [mid, right @ ..] = right;
        *first = twiddle_first(*first) * 2.;
        twiddle(TWIDDLES_128, left, right);
        *mid = twiddle_mid(*mid);
        signal.iter_mut().for_each(|s| *s *= 0.5);
        super::cast_real_to_complex(signal)
    }

    fn inverse_real_dft<'a>(&mut self, spectrum: &'a mut ComplexSpectrum) -> &'a mut Signal {
        let [left, right] = split_spectrum(spectrum);
        let [first, left @ ..] = left;
        let [mid, right @ ..] = right;
        *first = twiddle_first(*first);
        twiddle_inv(TWIDDLES_128, left, right);
        *mid = twiddle_mid(*mid);
        self.0(spectrum);
        let signal = super::cast_complex_to_real(spectrum);
        const MULT: f32 = 1. / FFT_SIZE as f32;
        signal.iter_mut().for_each(|s| *s *= MULT);
        signal
    }
}

type SplitSpectrum = [[C; FFT_SIZE_BY_2 / 2]; 2];

fn split_spectrum(spectrum: &mut ComplexSpectrum) -> &mut SplitSpectrum {
    static_assertions::assert_eq_size!(ComplexSpectrum, SplitSpectrum);
    static_assertions::assert_eq_align!(ComplexSpectrum, SplitSpectrum);
    unsafe { &mut *(spectrum as *mut _ as *mut _) }
}

fn twiddle_first(C { re, im }: C) -> C {
    C::new(re + im, re - im)
}

fn twiddle_mid(C { re, im }: C) -> C {
    C::new(re + re, -im - im)
}

fn twiddle<const N: usize>(twiddles: &[C; N], left: &mut [C; N], right: &mut [C; N]) {
    for ((&twiddle, left), right) in
        twiddles.iter().zip(left.iter_mut()).zip(right.iter_mut().rev())
    {
        let sum = *left + *right;
        let diff = *left - *right;
        let twiddled_re_sum = sum * twiddle.re;
        let twiddled_im_sum = sum * twiddle.im;
        let twiddled_re_diff = diff * twiddle.re;
        let twiddled_im_diff = diff * twiddle.im;
        let output_twiddled_real = twiddled_re_sum.im + twiddled_im_diff.re;
        let output_twiddled_im = twiddled_im_sum.im - twiddled_re_diff.re;
        *left = C { re: sum.re + output_twiddled_real, im: diff.im + output_twiddled_im };
        *right = C { re: sum.re - output_twiddled_real, im: output_twiddled_im - diff.im };
    }
}

fn twiddle_inv<const N: usize>(twiddles: &[C; N], left: &mut [C; N], right: &mut [C; N]) {
    for ((&twiddle, left), right) in
        twiddles.iter().zip(left.iter_mut()).zip(right.iter_mut().rev())
    {
        let sum = *left + *right;
        let diff = *left - *right;
        let twiddled_re_sum = sum * twiddle.re;
        let twiddled_im_sum = sum * twiddle.im;
        let twiddled_re_diff = diff * twiddle.re;
        let twiddled_im_diff = diff * twiddle.im;
        let output_twiddled_real = twiddled_re_sum.im - twiddled_im_diff.re;
        let output_twiddled_im = twiddled_im_sum.im + twiddled_re_diff.re;
        *left = C { re: sum.re + output_twiddled_real, im: output_twiddled_im - diff.im };
        *right = C { re: sum.re - output_twiddled_real, im: diff.im + output_twiddled_im };
    }
}

#[rustfmt::skip]
const TWIDDLES_128: &[C; 63] = &[
    C::new(0.9996988, -0.024541229), C::new(0.99879545, -0.049067676), C::new(0.99729043, -0.07356457),
    C::new(0.9951847, -0.09801714),  C::new(0.99247956, -0.12241068),  C::new(0.9891765, -0.14673047),
    C::new(0.98527765, -0.17096189), C::new(0.98078525, -0.19509032),  C::new(0.9757021, -0.21910124),
    C::new(0.97003126, -0.24298018), C::new(0.96377605, -0.26671275),  C::new(0.95694035, -0.29028466),
    C::new(0.94952816, -0.31368175), C::new(0.94154406, -0.33688986),  C::new(0.9329928, -0.35989505),
    C::new(0.9238795, -0.38268343),  C::new(0.9142098, -0.4052413),    C::new(0.9039893, -0.42755508),
    C::new(0.8932243, -0.44961134),  C::new(0.8819213, -0.47139674),   C::new(0.87008697, -0.4928982),
    C::new(0.8577286, -0.51410276),  C::new(0.8448536, -0.53499764),   C::new(0.8314696, -0.55557024),
    C::new(0.8175848, -0.57580817),  C::new(0.8032075, -0.5956993),    C::new(0.7883464, -0.6152316),
    C::new(0.77301043, -0.6343933),  C::new(0.7572088, -0.65317285),   C::new(0.7409511, -0.671559),
    C::new(0.7242471, -0.68954057),  C::new(0.70710677, -0.70710677),  C::new(0.68954057, -0.7242471),
    C::new(0.671559, -0.7409511),    C::new(0.65317285, -0.7572088),   C::new(0.6343933, -0.77301043),
    C::new(0.6152316, -0.7883464),   C::new(0.5956993, -0.8032075),    C::new(0.57580817, -0.8175848),
    C::new(0.55557024, -0.8314696),  C::new(0.53499764, -0.8448536),   C::new(0.51410276, -0.8577286),
    C::new(0.4928982, -0.87008697),  C::new(0.47139674, -0.8819213),   C::new(0.44961134, -0.8932243),
    C::new(0.42755508, -0.9039893),  C::new(0.4052413, -0.9142098),    C::new(0.38268343, -0.9238795),
    C::new(0.35989505, -0.9329928),  C::new(0.33688986, -0.94154406),  C::new(0.31368175, -0.94952816),
    C::new(0.29028466, -0.95694035), C::new(0.26671275, -0.96377605),  C::new(0.24298018, -0.97003126),
    C::new(0.21910124, -0.9757021),  C::new(0.19509032, -0.98078525),  C::new(0.17096189, -0.98527765),
    C::new(0.14673047, -0.9891765),  C::new(0.12241068, -0.99247956),  C::new(0.09801714, -0.9951847),
    C::new(0.07356457, -0.99729043), C::new(0.049067676, -0.99879545), C::new(0.024541229, -0.9996988),
];
