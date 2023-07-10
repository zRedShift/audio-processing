#[derive(Copy, Clone)]
pub struct PriorSignalModel {
    pub lrt: f32,
    pub flatness_threshold: f32,
    pub template_diff_threshold: f32,
    pub lrt_weighting: f32,
    pub flatness_weighting: f32,
    pub difference_weighting: f32,
}

impl PriorSignalModel {
    pub const fn new(lrt_initial_value: f32) -> Self {
        Self {
            lrt: lrt_initial_value,
            flatness_threshold: 0.5,
            template_diff_threshold: 0.5,
            lrt_weighting: 1.,
            flatness_weighting: 0.,
            difference_weighting: 0.,
        }
    }
}
