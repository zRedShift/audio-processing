use const_soft_float::soft_f64::SoftF64;
use num_complex::Complex32 as C;

pub struct ComplexDftAdapter<F, H> {
    fft: F,
    helper: H,
}

#[macro_export]
macro_rules! complex_dft_adapter {
    ($vis:vis $name:ident, $x:expr) => {
        $vis type $name<F> = $crate::ComplexDftAdapter<F, $crate::FftHelper<{$x / 2}, $x, {$x / 4 - 1}>>;
    };
}

impl<
        F: FnMut(&mut [C; COMPLEX]),
        const COMPLEX: usize,
        const REAL: usize,
        const TWIDDLE: usize,
    > ComplexDftAdapter<F, FftHelper<COMPLEX, REAL, TWIDDLE>>
{
    const MULT: f32 = 1. / REAL as f32;

    pub fn new(fft: F) -> Self {
        let helper = FftHelper::<COMPLEX, REAL, TWIDDLE>;
        helper.assert();
        Self { fft, helper }
    }

    pub fn real_dft<'a>(&mut self, signal: &'a mut [f32; REAL]) -> &'a mut [C; COMPLEX] {
        let spectrum = self.helper.cast_real_to_complex(signal);
        (self.fft)(spectrum);
        self.helper.twiddle(spectrum);
        signal.iter_mut().for_each(|s| *s *= 0.5);
        self.helper.cast_real_to_complex(signal)
    }

    pub fn inverse_real_dft<'a>(&mut self, spectrum: &'a mut [C; COMPLEX]) -> &'a mut [f32; REAL] {
        self.helper.twiddle_inv(spectrum);
        (self.fft)(spectrum);
        let signal = self.helper.cast_complex_to_real(spectrum);
        signal.iter_mut().for_each(|s| *s *= Self::MULT);
        signal
    }
}

pub struct FftHelper<const N: usize, const M: usize, const K: usize>;

impl<const N: usize, const M: usize, const K: usize> FftHelper<N, M, K> {
    const OK: () = assert!(
        N * 2 == M && N >= 4 && (K + 1) * 2 == N,
        "M has to be double the size of N, K must be (N / 2) - 1"
    );

    fn assert(&self) {
        #[allow(clippy::let_unit_value)]
        let () = Self::OK;
    }

    const TWIDDLE_CONST: f64 = -2f64 * core::f64::consts::PI / M as f64;
    const TWIDDLES: &'static [C; K] = &{
        let mut arr = [C { re: 0., im: 0. }; K];
        let mut i = 0;
        while i < K {
            arr[i] = Self::compute_twiddle(i + 1);
            i += 1;
        }
        arr
    };

    fn twiddle(&self, spectrum: &mut [C; N]) {
        let SplitSpectrum { first, left, mid, right } = Self::split_spectrum(spectrum);
        *first = Self::twiddle_first(*first) * 2.;
        for ((&twiddle, left), right) in Self::TWIDDLES.iter().zip(left).zip(right.iter_mut().rev())
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
        *mid = Self::twiddle_mid(*mid);
    }

    fn twiddle_inv(&self, spectrum: &mut [C; N]) {
        let SplitSpectrum { first, left, mid, right } = Self::split_spectrum(spectrum);
        *first = Self::twiddle_first(*first);
        for ((&twiddle, left), right) in Self::TWIDDLES.iter().zip(left).zip(right.iter_mut().rev())
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
        *mid = Self::twiddle_mid(*mid);
    }

    fn twiddle_first(C { re, im }: C) -> C {
        C::new(re + im, re - im)
    }

    fn twiddle_mid(C { re, im }: C) -> C {
        C::new(re + re, -im - im)
    }

    fn cast_complex_to_real<'a>(&self, spectrum: &'a mut [C; N]) -> &'a mut [f32; M] {
        unsafe { &mut *(spectrum as *mut _ as *mut _) }
    }

    fn cast_real_to_complex<'a>(&self, signal: &'a mut [f32; M]) -> &'a mut [C; N] {
        unsafe { &mut *(signal as *mut _ as *mut _) }
    }

    fn split_spectrum(spectrum: &mut [C; N]) -> &mut SplitSpectrum<K> {
        unsafe { &mut *(spectrum as *mut _ as *mut _) }
    }

    const fn compute_twiddle(index: usize) -> C {
        let angle = SoftF64(Self::TWIDDLE_CONST).mul(SoftF64(index as f64));
        C { re: angle.cos().0 as _, im: angle.sin().0 as _ }
    }
}

#[repr(C)]
struct SplitSpectrum<const N: usize> {
    first: C,
    left: [C; N],
    mid: C,
    right: [C; N],
}

#[cfg(test)]
mod test {
    use super::*;

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

    #[test]
    fn test_const_soft_float() {
        assert_eq!(FftHelper::<128, 256, 63>::TWIDDLES, TWIDDLES_128);
    }
}
