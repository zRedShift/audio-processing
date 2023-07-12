// A-law is basically as follows:
//
//      Linear Input Code        Compressed Code
//      -----------------        ---------------
//      0000000wxyza             000wxyz
//      0000001wxyza             001wxyz
//      000001wxyzab             010wxyz
//      00001wxyzabc             011wxyz
//      0001wxyzabcd             100wxyz
//      001wxyzabcde             101wxyz
//      01wxyzabcdef             110wxyz
//      1wxyzabcdefg             111wxyz
//
// For further information see John C. Bellamy's Digital Telephony, 1982,
// John Wiley & Sons, pps 98-111 and 472-476.
const ALAW_AMI_MASK: u8 = 0x55;

#[inline]
pub fn linear_to_alaw(linear: i16) -> u8 {
    let (linear, mask) = if linear >= 0 {
        // Sign (bit 7) bit = 1
        (linear as i32, ALAW_AMI_MASK | 0x80)
    } else {
        // Sign (bit 7) bit = 0
        (-(linear as i32) - 1, ALAW_AMI_MASK)
    };

    // Convert the scaled magnitude to segment number.
    let seg = util::top_bit(linear as u32 | 0xFF) - 7;
    if seg >= 8 {
        if linear >= 0 {
            // Out of range. Return maximum value
            return 0x7F ^ mask;
        }
        // We must be just a tiny step below zero
        return mask;
    }

    // Combine the sign, segment, and quantization bits.
    ((seg << 4) | ((linear >> (if seg != 0 { seg + 3 } else { 4 })) & 0x0F)) as u8 ^ mask
}

#[inline]
pub fn alaw_to_linear(mut alaw: u8) -> i16 {
    alaw ^= ALAW_AMI_MASK;
    let i = (alaw as i32 & 0x0F) << 4;
    let seg = (alaw as i32 & 0x70) >> 4;
    let i = if seg != 0 { (i + 0x108) << (seg - 1) } else { i + 8 };
    (if alaw & 0x80 != 0 { i } else { -i }) as i16
}
