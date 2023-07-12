#![cfg_attr(feature = "nightly", feature(asm_experimental_arch))]
#![no_std]

#[inline]
pub fn top_bit(bits: u32) -> i32 {
    31 - leading_zeros(bits) as i32
}

#[inline]
pub fn leading_zeros(bits: u32) -> u32 {
    let lz;
    #[cfg(all(feature = "nightly", target_arch = "xtensa", target_feature = "nsa"))]
    unsafe {
        core::arch::asm!("nsau {0}, {1}", out(reg) lz, in(reg) b, options(preserves_flags))
    }
    #[cfg(not(all(feature = "nightly", target_arch = "xtensa", target_feature = "nsa")))]
    {
        lz = bits.leading_zeros()
    }
    lz
}
