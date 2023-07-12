// Mu-law is basically as follows:
//
//      Biased Linear Input Code        Compressed Code
//      ------------------------        ---------------
//      00000001wxyza                   000wxyz
//      0000001wxyzab                   001wxyz
//      000001wxyzabc                   010wxyz
//      00001wxyzabcd                   011wxyz
//      0001wxyzabcde                   100wxyz
//      001wxyzabcdef                   101wxyz
//      01wxyzabcdefg                   110wxyz
//      1wxyzabcdefgh                   111wxyz
//
// Each biased linear code has a leading 1 which identifies the segment
// number. The value of the segment number is equal to 7 minus the number
// of leading 0's. The quantization interval is directly available as the
// four bits wxyz.  * The trailing bits (a - h) are ignored.
//
// Ordinarily the complement of the resulting code word is used for
// transmission, and so the code word is complemented before it is returned.
//
// For further information see John C. Bellamy's Digital Telephony, 1982,
// John Wiley & Sons, pps 98-111 and 472-476.
const ULAW_BIAS: i32 = 0x84;

#[inline]
pub fn linear_to_ulaw(linear: i16) -> u8 {
    // Get the sign and the magnitude of the value
    let (linear, mask) = if linear < 0 {
        (ULAW_BIAS - linear as i32 - 1, 0x7F)
    } else {
        (ULAW_BIAS + linear as i32, 0xFF)
    };

    // Combine the sign, segment, quantization bits,
    // and complement the code word.
    let seg = util::top_bit(linear as u32 | 0xFF) - 7;

    let u_val = if seg >= 8 { 0x7F } else { (seg << 4) | ((linear >> (seg + 3)) & 0xF) };
    (u_val ^ mask) as u8
}

#[inline]
pub fn ulaw_to_linear(mut ulaw: u8) -> i16 {
    // Complement to obtain normal u-law value
    ulaw = !ulaw;
    // Extract and bias the quantization bits. Then
    // shift up by the segment number and subtract out the bias.
    let t = (((ulaw & 0x0F) << 3) as i32 + ULAW_BIAS) << ((ulaw & 0x70) >> 4);
    (if ulaw & 0x80 != 0 { ULAW_BIAS - t } else { t - ULAW_BIAS }) as i16
}
