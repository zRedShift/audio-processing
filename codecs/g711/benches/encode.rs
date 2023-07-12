use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use g711::linear_to_ulaw;
use rand::Rng;

fn bench_sketch(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let mut arr = [0i16; 160];
    let mut out = [0u8; 160];
    rng.fill(&mut arr[..]);
    let mut group = c.benchmark_group("linear_to_ulaw");
    for i in 0..1000 {
        group.bench_function(BenchmarkId::from_parameter(i), |b| {
            b.iter(|| {
                let arr = black_box(&mut arr);
                let out = black_box(&mut out);
                for i in 0..160 {
                    out[i] = linear_to_ulaw(arr[i]);
                }
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_sketch);
criterion_main!(benches);
