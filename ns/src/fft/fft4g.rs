#![allow(unused_assignments)]
use core::ffi::{c_double, c_float, c_int};

#[allow(non_camel_case_types)]
pub type size_t = usize;

unsafe fn makewt(nw: size_t, ip: *mut size_t, w: *mut c_float) {
    let mut j: size_t = 0;
    let mut nwh: size_t = 0;
    let mut delta: c_float = 0.;
    let mut x: c_float = 0.;
    let mut y: c_float = 0.;
    *ip.offset(0 as c_int as isize) = nw;
    *ip.offset(1 as c_int as isize) = 1 as c_int as size_t;
    if nw > 2 as c_int as size_t {
        nwh = nw >> 1 as c_int;
        delta = 1.0f32.atan() / nwh as c_float;
        *w.offset(0 as c_int as isize) = 1 as c_int as c_float;
        *w.offset(1 as c_int as isize) = 0 as c_int as c_float;
        *w.add(nwh) = ((delta * nwh as c_float) as c_double).cos() as c_float;
        *w.add(nwh.wrapping_add(1 as c_int as size_t)) = *w.add(nwh);
        if nwh > 2 as c_int as size_t {
            j = 2 as c_int as size_t;
            while j < nwh {
                x = ((delta * j as c_float) as c_double).cos() as c_float;
                y = ((delta * j as c_float) as c_double).sin() as c_float;
                *w.add(j) = x;
                *w.add(j.wrapping_add(1 as c_int as size_t)) = y;
                *w.add(nw.wrapping_sub(j)) = y;
                *w.add(nw.wrapping_sub(j).wrapping_add(1 as c_int as size_t)) = x;
                j = (j as size_t).wrapping_add(2 as c_int as size_t) as size_t as size_t;
            }
            bitrv2(nw, ip.offset(2 as c_int as isize), w);
        }
    }
}

unsafe fn makect(nc: size_t, ip: *mut size_t, c: *mut c_float) {
    let mut j: size_t = 0;
    let mut nch: size_t = 0;
    let mut delta: c_float = 0.;
    *ip.offset(1 as c_int as isize) = nc;
    if nc > 1 as c_int as size_t {
        nch = nc >> 1 as c_int;
        delta = 1.0f32.atan() / nch as c_float;
        *c.offset(0 as c_int as isize) = ((delta * nch as c_float) as c_double).cos() as c_float;
        *c.add(nch) = 0.5f32 * *c.offset(0 as c_int as isize);
        j = 1 as c_int as size_t;
        while j < nch {
            *c.add(j) = 0.5f32 * ((delta * j as c_float) as c_double).cos() as c_float;
            *c.add(nc.wrapping_sub(j)) =
                0.5f32 * ((delta * j as c_float) as c_double).sin() as c_float;
            j = j.wrapping_add(1);
        }
    }
}

unsafe fn bitrv2(n: size_t, ip: *mut size_t, a: *mut c_float) {
    let mut j: size_t = 0;
    let mut j1: size_t = 0;
    let mut k: size_t = 0;
    let mut k1: size_t = 0;
    let mut l: size_t = 0;
    let mut m: size_t = 0;
    let mut m2: size_t = 0;
    let mut xr: c_float = 0.;
    let mut xi: c_float = 0.;
    let mut yr: c_float = 0.;
    let mut yi: c_float = 0.;
    *ip.offset(0 as c_int as isize) = 0 as c_int as size_t;
    l = n;
    m = 1 as c_int as size_t;
    while (m << 3 as c_int) < l {
        l >>= 1 as c_int;
        j = 0 as c_int as size_t;
        while j < m {
            *ip.add(m.wrapping_add(j)) = (*ip.add(j)).wrapping_add(l);
            j = j.wrapping_add(1);
        }
        m <<= 1 as c_int;
    }
    m2 = (2 as c_int as size_t).wrapping_mul(m);
    if m << 3 as c_int == l {
        k = 0 as c_int as size_t;
        while k < m {
            j = 0 as c_int as size_t;
            while j < k {
                j1 = (2 as c_int as size_t).wrapping_mul(j).wrapping_add(*ip.add(k));
                k1 = (2 as c_int as size_t).wrapping_mul(k).wrapping_add(*ip.add(j));
                xr = *a.add(j1);
                xi = *a.add(j1.wrapping_add(1 as c_int as size_t));
                yr = *a.add(k1);
                yi = *a.add(k1.wrapping_add(1 as c_int as size_t));
                *a.add(j1) = yr;
                *a.add(j1.wrapping_add(1 as c_int as size_t)) = yi;
                *a.add(k1) = xr;
                *a.add(k1.wrapping_add(1 as c_int as size_t)) = xi;
                j1 = (j1 as size_t).wrapping_add(m2) as size_t as size_t;
                k1 = (k1 as size_t).wrapping_add((2 as c_int as size_t).wrapping_mul(m2)) as size_t
                    as size_t;
                xr = *a.add(j1);
                xi = *a.add(j1.wrapping_add(1 as c_int as size_t));
                yr = *a.add(k1);
                yi = *a.add(k1.wrapping_add(1 as c_int as size_t));
                *a.add(j1) = yr;
                *a.add(j1.wrapping_add(1 as c_int as size_t)) = yi;
                *a.add(k1) = xr;
                *a.add(k1.wrapping_add(1 as c_int as size_t)) = xi;
                j1 = (j1 as size_t).wrapping_add(m2) as size_t as size_t;
                k1 = (k1 as size_t).wrapping_sub(m2) as size_t as size_t;
                xr = *a.add(j1);
                xi = *a.add(j1.wrapping_add(1 as c_int as size_t));
                yr = *a.add(k1);
                yi = *a.add(k1.wrapping_add(1 as c_int as size_t));
                *a.add(j1) = yr;
                *a.add(j1.wrapping_add(1 as c_int as size_t)) = yi;
                *a.add(k1) = xr;
                *a.add(k1.wrapping_add(1 as c_int as size_t)) = xi;
                j1 = (j1 as size_t).wrapping_add(m2) as size_t as size_t;
                k1 = (k1 as size_t).wrapping_add((2 as c_int as size_t).wrapping_mul(m2)) as size_t
                    as size_t;
                xr = *a.add(j1);
                xi = *a.add(j1.wrapping_add(1 as c_int as size_t));
                yr = *a.add(k1);
                yi = *a.add(k1.wrapping_add(1 as c_int as size_t));
                *a.add(j1) = yr;
                *a.add(j1.wrapping_add(1 as c_int as size_t)) = yi;
                *a.add(k1) = xr;
                *a.add(k1.wrapping_add(1 as c_int as size_t)) = xi;
                j = j.wrapping_add(1);
            }
            j1 = (2 as c_int as size_t).wrapping_mul(k).wrapping_add(m2).wrapping_add(*ip.add(k));
            k1 = j1.wrapping_add(m2);
            xr = *a.add(j1);
            xi = *a.add(j1.wrapping_add(1 as c_int as size_t));
            yr = *a.add(k1);
            yi = *a.add(k1.wrapping_add(1 as c_int as size_t));
            *a.add(j1) = yr;
            *a.add(j1.wrapping_add(1 as c_int as size_t)) = yi;
            *a.add(k1) = xr;
            *a.add(k1.wrapping_add(1 as c_int as size_t)) = xi;
            k = k.wrapping_add(1);
        }
    } else {
        k = 1 as c_int as size_t;
        while k < m {
            j = 0 as c_int as size_t;
            while j < k {
                j1 = (2 as c_int as size_t).wrapping_mul(j).wrapping_add(*ip.add(k));
                k1 = (2 as c_int as size_t).wrapping_mul(k).wrapping_add(*ip.add(j));
                xr = *a.add(j1);
                xi = *a.add(j1.wrapping_add(1 as c_int as size_t));
                yr = *a.add(k1);
                yi = *a.add(k1.wrapping_add(1 as c_int as size_t));
                *a.add(j1) = yr;
                *a.add(j1.wrapping_add(1 as c_int as size_t)) = yi;
                *a.add(k1) = xr;
                *a.add(k1.wrapping_add(1 as c_int as size_t)) = xi;
                j1 = (j1 as size_t).wrapping_add(m2) as size_t as size_t;
                k1 = (k1 as size_t).wrapping_add(m2) as size_t as size_t;
                xr = *a.add(j1);
                xi = *a.add(j1.wrapping_add(1 as c_int as size_t));
                yr = *a.add(k1);
                yi = *a.add(k1.wrapping_add(1 as c_int as size_t));
                *a.add(j1) = yr;
                *a.add(j1.wrapping_add(1 as c_int as size_t)) = yi;
                *a.add(k1) = xr;
                *a.add(k1.wrapping_add(1 as c_int as size_t)) = xi;
                j = j.wrapping_add(1);
            }
            k = k.wrapping_add(1);
        }
    };
}

unsafe fn cftfsub(n: size_t, a: *mut c_float, w: *mut c_float) {
    let mut j: size_t = 0;
    let mut j1: size_t = 0;
    let mut j2: size_t = 0;
    let mut j3: size_t = 0;
    let mut l: size_t = 0;
    let mut x0r: c_float = 0.;
    let mut x0i: c_float = 0.;
    let mut x1r: c_float = 0.;
    let mut x1i: c_float = 0.;
    let mut x2r: c_float = 0.;
    let mut x2i: c_float = 0.;
    let mut x3r: c_float = 0.;
    let mut x3i: c_float = 0.;
    l = 2 as c_int as size_t;
    if n > 8 as c_int as size_t {
        cft1st(n, a, w);
        l = 8 as c_int as size_t;
        while (l << 2 as c_int) < n {
            cftmdl(n, l, a, w);
            l <<= 2 as c_int;
        }
    }
    if l << 2 as c_int == n {
        j = 0 as c_int as size_t;
        while j < l {
            j1 = j.wrapping_add(l);
            j2 = j1.wrapping_add(l);
            j3 = j2.wrapping_add(l);
            x0r = *a.add(j) + *a.add(j1);
            x0i = *a.add(j.wrapping_add(1 as c_int as size_t))
                + *a.add(j1.wrapping_add(1 as c_int as size_t));
            x1r = *a.add(j) - *a.add(j1);
            x1i = *a.add(j.wrapping_add(1 as c_int as size_t))
                - *a.add(j1.wrapping_add(1 as c_int as size_t));
            x2r = *a.add(j2) + *a.add(j3);
            x2i = *a.add(j2.wrapping_add(1 as c_int as size_t))
                + *a.add(j3.wrapping_add(1 as c_int as size_t));
            x3r = *a.add(j2) - *a.add(j3);
            x3i = *a.add(j2.wrapping_add(1 as c_int as size_t))
                - *a.add(j3.wrapping_add(1 as c_int as size_t));
            *a.add(j) = x0r + x2r;
            *a.add(j.wrapping_add(1 as c_int as size_t)) = x0i + x2i;
            *a.add(j2) = x0r - x2r;
            *a.add(j2.wrapping_add(1 as c_int as size_t)) = x0i - x2i;
            *a.add(j1) = x1r - x3i;
            *a.add(j1.wrapping_add(1 as c_int as size_t)) = x1i + x3r;
            *a.add(j3) = x1r + x3i;
            *a.add(j3.wrapping_add(1 as c_int as size_t)) = x1i - x3r;
            j = (j as size_t).wrapping_add(2 as c_int as size_t) as size_t as size_t;
        }
    } else {
        j = 0 as c_int as size_t;
        while j < l {
            j1 = j.wrapping_add(l);
            x0r = *a.add(j) - *a.add(j1);
            x0i = *a.add(j.wrapping_add(1 as c_int as size_t))
                - *a.add(j1.wrapping_add(1 as c_int as size_t));
            *a.add(j) += *a.add(j1);
            *a.add(j.wrapping_add(1 as c_int as size_t)) +=
                *a.add(j1.wrapping_add(1 as c_int as size_t));
            *a.add(j1) = x0r;
            *a.add(j1.wrapping_add(1 as c_int as size_t)) = x0i;
            j = (j as size_t).wrapping_add(2 as c_int as size_t) as size_t as size_t;
        }
    };
}

unsafe fn cftbsub(n: size_t, a: *mut c_float, w: *mut c_float) {
    let mut j: size_t = 0;
    let mut j1: size_t = 0;
    let mut j2: size_t = 0;
    let mut j3: size_t = 0;
    let mut l: size_t = 0;
    let mut x0r: c_float = 0.;
    let mut x0i: c_float = 0.;
    let mut x1r: c_float = 0.;
    let mut x1i: c_float = 0.;
    let mut x2r: c_float = 0.;
    let mut x2i: c_float = 0.;
    let mut x3r: c_float = 0.;
    let mut x3i: c_float = 0.;
    l = 2 as c_int as size_t;
    if n > 8 as c_int as size_t {
        cft1st(n, a, w);
        l = 8 as c_int as size_t;
        while (l << 2 as c_int) < n {
            cftmdl(n, l, a, w);
            l <<= 2 as c_int;
        }
    }
    if l << 2 as c_int == n {
        j = 0 as c_int as size_t;
        while j < l {
            j1 = j.wrapping_add(l);
            j2 = j1.wrapping_add(l);
            j3 = j2.wrapping_add(l);
            x0r = *a.add(j) + *a.add(j1);
            x0i = -*a.add(j.wrapping_add(1 as c_int as size_t))
                - *a.add(j1.wrapping_add(1 as c_int as size_t));
            x1r = *a.add(j) - *a.add(j1);
            x1i = -*a.add(j.wrapping_add(1 as c_int as size_t))
                + *a.add(j1.wrapping_add(1 as c_int as size_t));
            x2r = *a.add(j2) + *a.add(j3);
            x2i = *a.add(j2.wrapping_add(1 as c_int as size_t))
                + *a.add(j3.wrapping_add(1 as c_int as size_t));
            x3r = *a.add(j2) - *a.add(j3);
            x3i = *a.add(j2.wrapping_add(1 as c_int as size_t))
                - *a.add(j3.wrapping_add(1 as c_int as size_t));
            *a.add(j) = x0r + x2r;
            *a.add(j.wrapping_add(1 as c_int as size_t)) = x0i - x2i;
            *a.add(j2) = x0r - x2r;
            *a.add(j2.wrapping_add(1 as c_int as size_t)) = x0i + x2i;
            *a.add(j1) = x1r - x3i;
            *a.add(j1.wrapping_add(1 as c_int as size_t)) = x1i - x3r;
            *a.add(j3) = x1r + x3i;
            *a.add(j3.wrapping_add(1 as c_int as size_t)) = x1i + x3r;
            j = (j as size_t).wrapping_add(2 as c_int as size_t) as size_t as size_t;
        }
    } else {
        j = 0 as c_int as size_t;
        while j < l {
            j1 = j.wrapping_add(l);
            x0r = *a.add(j) - *a.add(j1);
            x0i = -*a.add(j.wrapping_add(1 as c_int as size_t))
                + *a.add(j1.wrapping_add(1 as c_int as size_t));
            *a.add(j) += *a.add(j1);
            *a.add(j.wrapping_add(1 as c_int as size_t)) = -*a
                .add(j.wrapping_add(1 as c_int as size_t))
                - *a.add(j1.wrapping_add(1 as c_int as size_t));
            *a.add(j1) = x0r;
            *a.add(j1.wrapping_add(1 as c_int as size_t)) = x0i;
            j = (j as size_t).wrapping_add(2 as c_int as size_t) as size_t as size_t;
        }
    };
}

unsafe fn cft1st(n: size_t, a: *mut c_float, w: *mut c_float) {
    let mut j: size_t = 0;
    let mut k1: size_t = 0;
    let mut k2: size_t = 0;
    let mut wk1r: c_float = 0.;
    let mut wk1i: c_float = 0.;
    let mut wk2r: c_float = 0.;
    let mut wk2i: c_float = 0.;
    let mut wk3r: c_float = 0.;
    let mut wk3i: c_float = 0.;
    let mut x0r: c_float = 0.;
    let mut x0i: c_float = 0.;
    let mut x1r: c_float = 0.;
    let mut x1i: c_float = 0.;
    let mut x2r: c_float = 0.;
    let mut x2i: c_float = 0.;
    let mut x3r: c_float = 0.;
    let mut x3i: c_float = 0.;
    x0r = *a.offset(0 as c_int as isize) + *a.offset(2 as c_int as isize);
    x0i = *a.offset(1 as c_int as isize) + *a.offset(3 as c_int as isize);
    x1r = *a.offset(0 as c_int as isize) - *a.offset(2 as c_int as isize);
    x1i = *a.offset(1 as c_int as isize) - *a.offset(3 as c_int as isize);
    x2r = *a.offset(4 as c_int as isize) + *a.offset(6 as c_int as isize);
    x2i = *a.offset(5 as c_int as isize) + *a.offset(7 as c_int as isize);
    x3r = *a.offset(4 as c_int as isize) - *a.offset(6 as c_int as isize);
    x3i = *a.offset(5 as c_int as isize) - *a.offset(7 as c_int as isize);
    *a.offset(0 as c_int as isize) = x0r + x2r;
    *a.offset(1 as c_int as isize) = x0i + x2i;
    *a.offset(4 as c_int as isize) = x0r - x2r;
    *a.offset(5 as c_int as isize) = x0i - x2i;
    *a.offset(2 as c_int as isize) = x1r - x3i;
    *a.offset(3 as c_int as isize) = x1i + x3r;
    *a.offset(6 as c_int as isize) = x1r + x3i;
    *a.offset(7 as c_int as isize) = x1i - x3r;
    wk1r = *w.offset(2 as c_int as isize);
    x0r = *a.offset(8 as c_int as isize) + *a.offset(10 as c_int as isize);
    x0i = *a.offset(9 as c_int as isize) + *a.offset(11 as c_int as isize);
    x1r = *a.offset(8 as c_int as isize) - *a.offset(10 as c_int as isize);
    x1i = *a.offset(9 as c_int as isize) - *a.offset(11 as c_int as isize);
    x2r = *a.offset(12 as c_int as isize) + *a.offset(14 as c_int as isize);
    x2i = *a.offset(13 as c_int as isize) + *a.offset(15 as c_int as isize);
    x3r = *a.offset(12 as c_int as isize) - *a.offset(14 as c_int as isize);
    x3i = *a.offset(13 as c_int as isize) - *a.offset(15 as c_int as isize);
    *a.offset(8 as c_int as isize) = x0r + x2r;
    *a.offset(9 as c_int as isize) = x0i + x2i;
    *a.offset(12 as c_int as isize) = x2i - x0i;
    *a.offset(13 as c_int as isize) = x0r - x2r;
    x0r = x1r - x3i;
    x0i = x1i + x3r;
    *a.offset(10 as c_int as isize) = wk1r * (x0r - x0i);
    *a.offset(11 as c_int as isize) = wk1r * (x0r + x0i);
    x0r = x3i + x1r;
    x0i = x3r - x1i;
    *a.offset(14 as c_int as isize) = wk1r * (x0i - x0r);
    *a.offset(15 as c_int as isize) = wk1r * (x0i + x0r);
    k1 = 0 as c_int as size_t;
    j = 16 as c_int as size_t;
    while j < n {
        k1 = (k1 as size_t).wrapping_add(2 as c_int as size_t) as size_t as size_t;
        k2 = (2 as c_int as size_t).wrapping_mul(k1);
        wk2r = *w.add(k1);
        wk2i = *w.add(k1.wrapping_add(1 as c_int as size_t));
        wk1r = *w.add(k2);
        wk1i = *w.add(k2.wrapping_add(1 as c_int as size_t));
        wk3r = wk1r - 2 as c_int as c_float * wk2i * wk1i;
        wk3i = 2 as c_int as c_float * wk2i * wk1r - wk1i;
        x0r = *a.add(j) + *a.add(j.wrapping_add(2 as c_int as size_t));
        x0i = *a.add(j.wrapping_add(1 as c_int as size_t))
            + *a.add(j.wrapping_add(3 as c_int as size_t));
        x1r = *a.add(j) - *a.add(j.wrapping_add(2 as c_int as size_t));
        x1i = *a.add(j.wrapping_add(1 as c_int as size_t))
            - *a.add(j.wrapping_add(3 as c_int as size_t));
        x2r = *a.add(j.wrapping_add(4 as c_int as size_t))
            + *a.add(j.wrapping_add(6 as c_int as size_t));
        x2i = *a.add(j.wrapping_add(5 as c_int as size_t))
            + *a.add(j.wrapping_add(7 as c_int as size_t));
        x3r = *a.add(j.wrapping_add(4 as c_int as size_t))
            - *a.add(j.wrapping_add(6 as c_int as size_t));
        x3i = *a.add(j.wrapping_add(5 as c_int as size_t))
            - *a.add(j.wrapping_add(7 as c_int as size_t));
        *a.add(j) = x0r + x2r;
        *a.add(j.wrapping_add(1 as c_int as size_t)) = x0i + x2i;
        x0r -= x2r;
        x0i -= x2i;
        *a.add(j.wrapping_add(4 as c_int as size_t)) = wk2r * x0r - wk2i * x0i;
        *a.add(j.wrapping_add(5 as c_int as size_t)) = wk2r * x0i + wk2i * x0r;
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        *a.add(j.wrapping_add(2 as c_int as size_t)) = wk1r * x0r - wk1i * x0i;
        *a.add(j.wrapping_add(3 as c_int as size_t)) = wk1r * x0i + wk1i * x0r;
        x0r = x1r + x3i;
        x0i = x1i - x3r;
        *a.add(j.wrapping_add(6 as c_int as size_t)) = wk3r * x0r - wk3i * x0i;
        *a.add(j.wrapping_add(7 as c_int as size_t)) = wk3r * x0i + wk3i * x0r;
        wk1r = *w.add(k2.wrapping_add(2 as c_int as size_t));
        wk1i = *w.add(k2.wrapping_add(3 as c_int as size_t));
        wk3r = wk1r - 2 as c_int as c_float * wk2r * wk1i;
        wk3i = 2 as c_int as c_float * wk2r * wk1r - wk1i;
        x0r = *a.add(j.wrapping_add(8 as c_int as size_t))
            + *a.add(j.wrapping_add(10 as c_int as size_t));
        x0i = *a.add(j.wrapping_add(9 as c_int as size_t))
            + *a.add(j.wrapping_add(11 as c_int as size_t));
        x1r = *a.add(j.wrapping_add(8 as c_int as size_t))
            - *a.add(j.wrapping_add(10 as c_int as size_t));
        x1i = *a.add(j.wrapping_add(9 as c_int as size_t))
            - *a.add(j.wrapping_add(11 as c_int as size_t));
        x2r = *a.add(j.wrapping_add(12 as c_int as size_t))
            + *a.add(j.wrapping_add(14 as c_int as size_t));
        x2i = *a.add(j.wrapping_add(13 as c_int as size_t))
            + *a.add(j.wrapping_add(15 as c_int as size_t));
        x3r = *a.add(j.wrapping_add(12 as c_int as size_t))
            - *a.add(j.wrapping_add(14 as c_int as size_t));
        x3i = *a.add(j.wrapping_add(13 as c_int as size_t))
            - *a.add(j.wrapping_add(15 as c_int as size_t));
        *a.add(j.wrapping_add(8 as c_int as size_t)) = x0r + x2r;
        *a.add(j.wrapping_add(9 as c_int as size_t)) = x0i + x2i;
        x0r -= x2r;
        x0i -= x2i;
        *a.add(j.wrapping_add(12 as c_int as size_t)) = -wk2i * x0r - wk2r * x0i;
        *a.add(j.wrapping_add(13 as c_int as size_t)) = -wk2i * x0i + wk2r * x0r;
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        *a.add(j.wrapping_add(10 as c_int as size_t)) = wk1r * x0r - wk1i * x0i;
        *a.add(j.wrapping_add(11 as c_int as size_t)) = wk1r * x0i + wk1i * x0r;
        x0r = x1r + x3i;
        x0i = x1i - x3r;
        *a.add(j.wrapping_add(14 as c_int as size_t)) = wk3r * x0r - wk3i * x0i;
        *a.add(j.wrapping_add(15 as c_int as size_t)) = wk3r * x0i + wk3i * x0r;
        j = (j as size_t).wrapping_add(16 as c_int as size_t) as size_t as size_t;
    }
}

unsafe fn cftmdl(n: size_t, l: size_t, a: *mut c_float, w: *mut c_float) {
    let mut j: size_t = 0;
    let mut j1: size_t = 0;
    let mut j2: size_t = 0;
    let mut j3: size_t = 0;
    let mut k: size_t = 0;
    let mut k1: size_t = 0;
    let mut k2: size_t = 0;
    let mut m: size_t = 0;
    let mut m2: size_t = 0;
    let mut wk1r: c_float = 0.;
    let mut wk1i: c_float = 0.;
    let mut wk2r: c_float = 0.;
    let mut wk2i: c_float = 0.;
    let mut wk3r: c_float = 0.;
    let mut wk3i: c_float = 0.;
    let mut x0r: c_float = 0.;
    let mut x0i: c_float = 0.;
    let mut x1r: c_float = 0.;
    let mut x1i: c_float = 0.;
    let mut x2r: c_float = 0.;
    let mut x2i: c_float = 0.;
    let mut x3r: c_float = 0.;
    let mut x3i: c_float = 0.;
    m = l << 2 as c_int;
    j = 0 as c_int as size_t;
    while j < l {
        j1 = j.wrapping_add(l);
        j2 = j1.wrapping_add(l);
        j3 = j2.wrapping_add(l);
        x0r = *a.add(j) + *a.add(j1);
        x0i = *a.add(j.wrapping_add(1 as c_int as size_t))
            + *a.add(j1.wrapping_add(1 as c_int as size_t));
        x1r = *a.add(j) - *a.add(j1);
        x1i = *a.add(j.wrapping_add(1 as c_int as size_t))
            - *a.add(j1.wrapping_add(1 as c_int as size_t));
        x2r = *a.add(j2) + *a.add(j3);
        x2i = *a.add(j2.wrapping_add(1 as c_int as size_t))
            + *a.add(j3.wrapping_add(1 as c_int as size_t));
        x3r = *a.add(j2) - *a.add(j3);
        x3i = *a.add(j2.wrapping_add(1 as c_int as size_t))
            - *a.add(j3.wrapping_add(1 as c_int as size_t));
        *a.add(j) = x0r + x2r;
        *a.add(j.wrapping_add(1 as c_int as size_t)) = x0i + x2i;
        *a.add(j2) = x0r - x2r;
        *a.add(j2.wrapping_add(1 as c_int as size_t)) = x0i - x2i;
        *a.add(j1) = x1r - x3i;
        *a.add(j1.wrapping_add(1 as c_int as size_t)) = x1i + x3r;
        *a.add(j3) = x1r + x3i;
        *a.add(j3.wrapping_add(1 as c_int as size_t)) = x1i - x3r;
        j = (j as size_t).wrapping_add(2 as c_int as size_t) as size_t as size_t;
    }
    wk1r = *w.offset(2 as c_int as isize);
    j = m;
    while j < l.wrapping_add(m) {
        j1 = j.wrapping_add(l);
        j2 = j1.wrapping_add(l);
        j3 = j2.wrapping_add(l);
        x0r = *a.add(j) + *a.add(j1);
        x0i = *a.add(j.wrapping_add(1 as c_int as size_t))
            + *a.add(j1.wrapping_add(1 as c_int as size_t));
        x1r = *a.add(j) - *a.add(j1);
        x1i = *a.add(j.wrapping_add(1 as c_int as size_t))
            - *a.add(j1.wrapping_add(1 as c_int as size_t));
        x2r = *a.add(j2) + *a.add(j3);
        x2i = *a.add(j2.wrapping_add(1 as c_int as size_t))
            + *a.add(j3.wrapping_add(1 as c_int as size_t));
        x3r = *a.add(j2) - *a.add(j3);
        x3i = *a.add(j2.wrapping_add(1 as c_int as size_t))
            - *a.add(j3.wrapping_add(1 as c_int as size_t));
        *a.add(j) = x0r + x2r;
        *a.add(j.wrapping_add(1 as c_int as size_t)) = x0i + x2i;
        *a.add(j2) = x2i - x0i;
        *a.add(j2.wrapping_add(1 as c_int as size_t)) = x0r - x2r;
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        *a.add(j1) = wk1r * (x0r - x0i);
        *a.add(j1.wrapping_add(1 as c_int as size_t)) = wk1r * (x0r + x0i);
        x0r = x3i + x1r;
        x0i = x3r - x1i;
        *a.add(j3) = wk1r * (x0i - x0r);
        *a.add(j3.wrapping_add(1 as c_int as size_t)) = wk1r * (x0i + x0r);
        j = (j as size_t).wrapping_add(2 as c_int as size_t) as size_t as size_t;
    }
    k1 = 0 as c_int as size_t;
    m2 = (2 as c_int as size_t).wrapping_mul(m);
    k = m2;
    while k < n {
        k1 = (k1 as size_t).wrapping_add(2 as c_int as size_t) as size_t as size_t;
        k2 = (2 as c_int as size_t).wrapping_mul(k1);
        wk2r = *w.add(k1);
        wk2i = *w.add(k1.wrapping_add(1 as c_int as size_t));
        wk1r = *w.add(k2);
        wk1i = *w.add(k2.wrapping_add(1 as c_int as size_t));
        wk3r = wk1r - 2 as c_int as c_float * wk2i * wk1i;
        wk3i = 2 as c_int as c_float * wk2i * wk1r - wk1i;
        j = k;
        while j < l.wrapping_add(k) {
            j1 = j.wrapping_add(l);
            j2 = j1.wrapping_add(l);
            j3 = j2.wrapping_add(l);
            x0r = *a.add(j) + *a.add(j1);
            x0i = *a.add(j.wrapping_add(1 as c_int as size_t))
                + *a.add(j1.wrapping_add(1 as c_int as size_t));
            x1r = *a.add(j) - *a.add(j1);
            x1i = *a.add(j.wrapping_add(1 as c_int as size_t))
                - *a.add(j1.wrapping_add(1 as c_int as size_t));
            x2r = *a.add(j2) + *a.add(j3);
            x2i = *a.add(j2.wrapping_add(1 as c_int as size_t))
                + *a.add(j3.wrapping_add(1 as c_int as size_t));
            x3r = *a.add(j2) - *a.add(j3);
            x3i = *a.add(j2.wrapping_add(1 as c_int as size_t))
                - *a.add(j3.wrapping_add(1 as c_int as size_t));
            *a.add(j) = x0r + x2r;
            *a.add(j.wrapping_add(1 as c_int as size_t)) = x0i + x2i;
            x0r -= x2r;
            x0i -= x2i;
            *a.add(j2) = wk2r * x0r - wk2i * x0i;
            *a.add(j2.wrapping_add(1 as c_int as size_t)) = wk2r * x0i + wk2i * x0r;
            x0r = x1r - x3i;
            x0i = x1i + x3r;
            *a.add(j1) = wk1r * x0r - wk1i * x0i;
            *a.add(j1.wrapping_add(1 as c_int as size_t)) = wk1r * x0i + wk1i * x0r;
            x0r = x1r + x3i;
            x0i = x1i - x3r;
            *a.add(j3) = wk3r * x0r - wk3i * x0i;
            *a.add(j3.wrapping_add(1 as c_int as size_t)) = wk3r * x0i + wk3i * x0r;
            j = (j as size_t).wrapping_add(2 as c_int as size_t) as size_t as size_t;
        }
        wk1r = *w.add(k2.wrapping_add(2 as c_int as size_t));
        wk1i = *w.add(k2.wrapping_add(3 as c_int as size_t));
        wk3r = wk1r - 2 as c_int as c_float * wk2r * wk1i;
        wk3i = 2 as c_int as c_float * wk2r * wk1r - wk1i;
        j = k.wrapping_add(m);
        while j < l.wrapping_add(k.wrapping_add(m)) {
            j1 = j.wrapping_add(l);
            j2 = j1.wrapping_add(l);
            j3 = j2.wrapping_add(l);
            x0r = *a.add(j) + *a.add(j1);
            x0i = *a.add(j.wrapping_add(1 as c_int as size_t))
                + *a.add(j1.wrapping_add(1 as c_int as size_t));
            x1r = *a.add(j) - *a.add(j1);
            x1i = *a.add(j.wrapping_add(1 as c_int as size_t))
                - *a.add(j1.wrapping_add(1 as c_int as size_t));
            x2r = *a.add(j2) + *a.add(j3);
            x2i = *a.add(j2.wrapping_add(1 as c_int as size_t))
                + *a.add(j3.wrapping_add(1 as c_int as size_t));
            x3r = *a.add(j2) - *a.add(j3);
            x3i = *a.add(j2.wrapping_add(1 as c_int as size_t))
                - *a.add(j3.wrapping_add(1 as c_int as size_t));
            *a.add(j) = x0r + x2r;
            *a.add(j.wrapping_add(1 as c_int as size_t)) = x0i + x2i;
            x0r -= x2r;
            x0i -= x2i;
            *a.add(j2) = -wk2i * x0r - wk2r * x0i;
            *a.add(j2.wrapping_add(1 as c_int as size_t)) = -wk2i * x0i + wk2r * x0r;
            x0r = x1r - x3i;
            x0i = x1i + x3r;
            *a.add(j1) = wk1r * x0r - wk1i * x0i;
            *a.add(j1.wrapping_add(1 as c_int as size_t)) = wk1r * x0i + wk1i * x0r;
            x0r = x1r + x3i;
            x0i = x1i - x3r;
            *a.add(j3) = wk3r * x0r - wk3i * x0i;
            *a.add(j3.wrapping_add(1 as c_int as size_t)) = wk3r * x0i + wk3i * x0r;
            j = (j as size_t).wrapping_add(2 as c_int as size_t) as size_t as size_t;
        }
        k = (k as size_t).wrapping_add(m2) as size_t as size_t;
    }
}

unsafe fn rftfsub(n: size_t, a: *mut c_float, nc: size_t, c: *mut c_float) {
    let mut j: size_t = 0;
    let mut k: size_t = 0;
    let mut kk: size_t = 0;
    let mut ks: size_t = 0;
    let mut m: size_t = 0;
    let mut wkr: c_float = 0.;
    let mut wki: c_float = 0.;
    let mut xr: c_float = 0.;
    let mut xi: c_float = 0.;
    let mut yr: c_float = 0.;
    let mut yi: c_float = 0.;
    m = n >> 1 as c_int;
    ks = (2 as c_int as size_t).wrapping_mul(nc).wrapping_div(m);
    kk = 0 as c_int as size_t;
    j = 2 as c_int as size_t;
    while j < m {
        k = n.wrapping_sub(j);
        kk = (kk as size_t).wrapping_add(ks) as size_t as size_t;
        wkr = 0.5f32 - *c.add(nc.wrapping_sub(kk));
        wki = *c.add(kk);
        xr = *a.add(j) - *a.add(k);
        xi = *a.add(j.wrapping_add(1 as c_int as size_t))
            + *a.add(k.wrapping_add(1 as c_int as size_t));
        yr = wkr * xr - wki * xi;
        yi = wkr * xi + wki * xr;
        *a.add(j) -= yr;
        *a.add(j.wrapping_add(1 as c_int as size_t)) -= yi;
        *a.add(k) += yr;
        *a.add(k.wrapping_add(1 as c_int as size_t)) -= yi;
        j = (j as size_t).wrapping_add(2 as c_int as size_t) as size_t as size_t;
    }
}

unsafe fn rftbsub(n: size_t, a: *mut c_float, nc: size_t, c: *mut c_float) {
    let mut j: size_t = 0;
    let mut k: size_t = 0;
    let mut kk: size_t = 0;
    let mut ks: size_t = 0;
    let mut m: size_t = 0;
    let mut wkr: c_float = 0.;
    let mut wki: c_float = 0.;
    let mut xr: c_float = 0.;
    let mut xi: c_float = 0.;
    let mut yr: c_float = 0.;
    let mut yi: c_float = 0.;
    *a.offset(1 as c_int as isize) = -*a.offset(1 as c_int as isize);
    m = n >> 1 as c_int;
    ks = (2 as c_int as size_t).wrapping_mul(nc).wrapping_div(m);
    kk = 0 as c_int as size_t;
    j = 2 as c_int as size_t;
    while j < m {
        k = n.wrapping_sub(j);
        kk = (kk as size_t).wrapping_add(ks) as size_t as size_t;
        wkr = 0.5f32 - *c.add(nc.wrapping_sub(kk));
        wki = *c.add(kk);
        xr = *a.add(j) - *a.add(k);
        xi = *a.add(j.wrapping_add(1 as c_int as size_t))
            + *a.add(k.wrapping_add(1 as c_int as size_t));
        yr = wkr * xr + wki * xi;
        yi = wkr * xi - wki * xr;
        *a.add(j) -= yr;
        *a.add(j.wrapping_add(1 as c_int as size_t)) =
            yi - *a.add(j.wrapping_add(1 as c_int as size_t));
        *a.add(k) += yr;
        *a.add(k.wrapping_add(1 as c_int as size_t)) =
            yi - *a.add(k.wrapping_add(1 as c_int as size_t));
        j = (j as size_t).wrapping_add(2 as c_int as size_t) as size_t as size_t;
    }
    *a.add(m.wrapping_add(1 as c_int as size_t)) = -*a.add(m.wrapping_add(1 as c_int as size_t));
}

pub unsafe fn rdft(n: size_t, isgn: c_int, a: *mut c_float, ip: *mut size_t, w: *mut c_float) {
    let mut nw: size_t = 0;
    let mut nc: size_t = 0;
    let mut xi: c_float = 0.;
    nw = *ip.offset(0 as c_int as isize);
    if n > nw << 2 as c_int {
        nw = n >> 2 as c_int;
        makewt(nw, ip, w);
    }
    nc = *ip.offset(1 as c_int as isize);
    if n > nc << 2 as c_int {
        nc = n >> 2 as c_int;
        makect(nc, ip, w.add(nw));
    }
    if isgn >= 0 as c_int {
        match n {
            5.. => {
                bitrv2(n, ip.offset(2 as c_int as isize), a);
                cftfsub(n, a, w);
                rftfsub(n, a, nc, w.add(nw));
            }
            4 => cftfsub(n, a, w),
            _ => {}
        }
        xi = *a.offset(0 as c_int as isize) - *a.offset(1 as c_int as isize);
        *a.offset(0 as c_int as isize) += *a.offset(1 as c_int as isize);
        *a.offset(1 as c_int as isize) = xi;
    } else {
        *a.offset(1 as c_int as isize) =
            0.5f32 * (*a.offset(0 as c_int as isize) - *a.offset(1 as c_int as isize));
        *a.offset(0 as c_int as isize) -= *a.offset(1 as c_int as isize);
        match n {
            5.. => {
                rftbsub(n, a, nc, w.add(nw));
                bitrv2(n, ip.offset(2 as c_int as isize), a);
                cftbsub(n, a, w);
            }
            4 => cftfsub(n, a, w),
            _ => {}
        }
    };
}
