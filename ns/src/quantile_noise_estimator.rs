use super::{Spectrum, FFT_SIZE_BY_2_PLUS_1, LONG_STARTUP_PHASE_BLOCKS};

const SIMULT: usize = 3;

#[derive(Copy, Clone)]
pub struct QuantileNoiseEstimator {
    density: [Spectrum; SIMULT],
    log_quantile: [Spectrum; SIMULT],
    quantile: Spectrum,
    counter: [i32; SIMULT],
    num_updates: i32,
}

impl QuantileNoiseEstimator {
    pub const fn new() -> Self {
        let mut counter = [0; SIMULT];
        let mut i = 0;
        while i < SIMULT {
            counter[i] = LONG_STARTUP_PHASE_BLOCKS * (i as i32 + 1) / SIMULT as i32;
            i += 1;
        }

        Self {
            density: [[0.3; FFT_SIZE_BY_2_PLUS_1]; SIMULT],
            log_quantile: [[8.; FFT_SIZE_BY_2_PLUS_1]; SIMULT],
            quantile: [0.; FFT_SIZE_BY_2_PLUS_1],
            counter,
            num_updates: 1,
        }
    }

    pub fn estimate(&mut self, signal_spectrum: &Spectrum, noise_spectrum: &mut Spectrum) {
        signal_spectrum
            .iter()
            .zip(noise_spectrum.iter_mut())
            .for_each(|(&signal_spectrum, noise_spectrum)| *noise_spectrum = signal_spectrum.ln());

        let mut quantile_index_to_return = SIMULT;

        // Loop over simultaneous estimates.
        for (s, ((counter, log_quantile), density)) in self
            .counter
            .iter_mut()
            .zip(self.log_quantile.iter_mut())
            .zip(self.density.iter_mut())
            .enumerate()
        {
            let one_by_counter_plus_1 = 1. / (*counter as f32 + 1.);

            for ((log_quantile, density), &log_spectrum) in
                log_quantile.iter_mut().zip(density.iter_mut()).zip(noise_spectrum.iter())
            {
                // Update log quantile estimate.
                let delta = if *density > 1. { 40. / *density } else { 40. };

                let multiplier = delta * one_by_counter_plus_1;
                if log_spectrum > *log_quantile {
                    *log_quantile += 0.25 * multiplier;
                } else {
                    *log_quantile -= 0.75 * multiplier;
                }

                // Update density estimate.
                const WIDTH: f32 = 0.01;
                const ONE_BY_WIDTH_PLUS_2: f32 = 1. / (2. * WIDTH);
                if (log_spectrum - *log_quantile).abs() < WIDTH {
                    *density =
                        (*counter as f32 * *density + ONE_BY_WIDTH_PLUS_2) * one_by_counter_plus_1;
                }
            }

            if *counter >= LONG_STARTUP_PHASE_BLOCKS {
                *counter = 0;
                if self.num_updates >= LONG_STARTUP_PHASE_BLOCKS {
                    quantile_index_to_return = s;
                }
            }

            *counter += 1;
        }

        // Sequentially update the noise during startup.
        if self.num_updates < LONG_STARTUP_PHASE_BLOCKS {
            // Use the last "s" to get noise during startup that differ from zero.
            quantile_index_to_return = SIMULT - 1;
            self.num_updates += 1;
        }

        if let Some(log_quantile) = self.log_quantile.get(quantile_index_to_return) {
            log_quantile
                .iter()
                .zip(self.quantile.iter_mut())
                .for_each(|(&log_quantile, quantile)| *quantile = log_quantile.exp());
        }

        *noise_spectrum = self.quantile;
    }
}
