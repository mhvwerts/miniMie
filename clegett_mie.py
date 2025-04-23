import numpy as np
from scipy.special import jv, yv
import warnings

def mie_ab(m: complex, x: float) -> np.ndarray:
    z = m * x
    nmax = np.int(np.round(2 + x + 4 * x**(1/3)))
    nmx = np.round(np.amax([nmax, np.abs(z)]) + 16).astype(int)
    n = np.arange(1, nmax + 1)
    nu = n + 0.5
    sx = np.sqrt(0.5 * np.pi * x)
    px = sx * jv(nu, x)
    p1x = np.append(np.sin(x), px[0:nmax-1])
    chx = -sx * yv(nu, x)
    ch1x = np.append(np.cos(x), chx[0:nmax-1])
    gsx = px - chx * 1j
    gs1x = p1x - ch1x * 1j
    dnx = np.zeros(nmx + 1) * 1j
    for j in np.arange(nmx, 1, -1):
        dnx[j - 1] = j / z - 1 / (dnx[j] + j / z)
    dn = dnx[n]
    da = dn / m + n / x
    db = m * dn + n / x

    an = (da * px - p1x) / (da * gsx - gs1x)
    bn = (db * px - p1x) / (db * gsx - gs1x)
    return np.array([an, bn])

def mie_abcd(m: complex, x: float) -> np.ndarray:
    nmax = np.int(np.round(2 + x + 4 * x**(1/3)))
    n = np.arange(1, nmax + 1).astype(int)
    nu = n + 0.5
    z = m * x
    m2 = m * m
    sqx = np.sqrt(0.5 * np.pi / x)
    sqz = np.sqrt(0.5 * np.pi / z)
    bx = jv(nu, x) * sqx
    bz = jv(nu, z) * sqz
    yx = yv(nu, x) * sqx
    hx = bx + yx * 1j
    b1x = np.append(np.sin(x) / x, bx[0:nmax - 1])
    b1z = np.append(np.sin(z) / z, bz[0:nmax - 1])
    y1x = np.append(-np.cos(x) / x, yx[0:nmax - 1])
    h1x = b1x + y1x * 1j
    ax = x * b1x - n * bx
    az = z * b1z - n * bz
    ahx = x * h1x - n * hx

    an = (m2 * bz * ax - bx * az) / (m2 * bz * ahx - hx * az)
    bn = (bz * ax - bx * az) / (bz * ahx - hx * az)
    cn = (bx * ahx - hx * ax) / (bz * ahx - hx * az)
    dn = m * (bx * ahx - hx * ax) / (m2 * bz * ahx - hx * az)
    return np.array([an, bn, cn, dn])

def mie(m: complex, x: float) -> np.ndarray:
    if m.real < 1:
        raise ValueError('invalid m. m.real must be greater than 1')
    elif m.imag < 0:
        raise ValueError('invalid m. m.imag must be greater than 0')
    if x < 0:
        raise ValueError('x must be >= 0')
    elif np.isclose(x, 0):
        return np.array([m.real, m.imag, 0, 0, 0, 0, 0, 0, 1.5])
    elif x > 0:
        nmax = np.int(np.round(2 + x + 4 * x**(1/3)))
        n = np.arange(1, nmax + 1)
        an, bn, _, _ = mie_abcd(m, x)
        tmp = np.zeros((4, nmax))
        tmp[0, 0:(nmax - 1)] = an[1:nmax].real
        tmp[1, 0:(nmax - 1)] = an[1:nmax].imag
        tmp[2, 0:(nmax - 1)] = bn[1:nmax].real
        tmp[3, 0:(nmax - 1)] = bn[1:nmax].imag
        qext = 2 * np.sum(((2 * n + 1) * (an.real + bn.real))) / x**2
        qsca = 2 * np.sum((2 * n + 1) * (np.abs(an)**2 + np.abs(bn)**2)) / x**2
        qabs = qext - qsca
        qb = np.abs(np.sum(((an - bn) * (2 * n + 1)) * (-1)**n))**2 / x**2
        asy1 = ((n * (n + 2) / (n + 1))
            * (an.real * tmp[0, :]
                + an.imag * tmp[1, :]
                + bn.real * tmp[2, :]
                + bn.imag * tmp[3,:]))
        asy2 = (((2 * n + 1) / n / (n + 1))
            * (an.real * bn.real + an.imag * bn.imag))
        asy = 4 / x**2 * np.sum(asy1 + asy2) / qsca
        qratio = qb / qsca
        return np.array([m.real, m.imag, x, qext, qsca, qabs, qb, asy, qratio])

def mie_pt(μ: float, nmax: int) -> np.ndarray:
    if (1 < μ) or (-1 > μ):
        raise ValueError('Invalid μ, must be -1 <= μ <= 1.')
    if nmax < 3:
        raise ValueError('Invalid nmax, must be nmax > 3')
    π = np.zeros((nmax + 1))
    τ = np.zeros((nmax + 1))
    π[1] = 1
    τ[1] = μ
    π[2] = 3 * μ
    τ[2] = 3 * np.cos(2 * np.arccos(μ))
    for n in range(3, nmax + 1):
        π1 = ((2 * n - 1) / (n - 1)) * π[n - 1] * μ
        π2 = (n / (n - 1)) * π[n - 2]
        π[n] = π1 - π2
        τ1 = n * μ * π[n]
        τ2 = (n + 1) * π[n - 1]
        τ[n] = τ1-τ2
    return np.array([π[1:], τ[1:]])

def mie_s12(m: complex, x: float, u: float) -> np.ndarray:
    nmax = np.int(np.round(2 + x + 4 * x**(1/3)))
    abcd = mie_abcd(m, x)
    an = abcd[0, :]
    bn = abcd[1, :]
    pt = mie_pt(u, nmax)
    pin = pt[0,:]
    tin = pt[1,:]
    n = np.array(range(1, nmax + 1))
    n2 = (2*n+1) / (n*(n+1))
    pin = n2 * pin
    tin = n2 * tin
    S1 = (an @ pin.transpose()) + (bn @ tin.transpose())
    S2 = (an @ tin.transpose()) + (bn @ pin.transpose())
    return np.array([S1, S2])

class InputDomainWarning(Warning):
    pass

def miecoated_ab1(m1: complex, m2: complex, x: float, y: float) -> np.ndarray:
    if (m1 * x).imag > 30:
        warnings.warn('Im(m1*x) > 30. This could result in bad solutions. See '
            'B&H p. 485. Compare the output of miecoated_ab1, miecoated_ab2, '
            'and miecoated_ab3 and be very careful if they give significantly '
            'different results.', InputDomainWarning)
    elif (m2 * x).imag > 30:
        warnings.warn('Im(m2*x) > 30. This could result in bad solutions. See '
            'B&H p. 485. Compare the output of miecoated_ab1, miecoated_ab2, '
            'and miecoated_ab3 and be very careful if they give significantly '
            'different results.', InputDomainWarning)
    elif (m2 * y).imag > 30:
        warnings.warn('Im(m2*y) > 30. This could result in bad solutions. See '
            'B&H p. 485. Compare the output of miecoated_ab1, miecoated_ab2, '
            'and miecoated_ab3 and be very careful if they give significantly '
            'different results.', InputDomainWarning)
    if x == y:
        raise ValueError('x == y, size parameters cannot be the same size.')
    elif x > y:
        raise ValueError('x > y, "inner" sphere larger than "outer"')
    m = m2 / m1
    u = m1 * x
    v = m2 * x
    w = m2 * y
    nmax = np.round(2 + y + 4 * y**(1/3)).astype(int)
    mx = np.amax([np.abs(m1 * y), np.abs(m2 * y)])
    nmx = np.round(np.amax([nmax, mx]) + 16).astype(int)
    nmax1 = nmax - 1
    n = np.arange(1, nmax + 1).astype(int)
    dnx = 1j * np.zeros(nmx)
    #I think this next line is unnecessary in python
    #dnx[nmx - 1] = complex(0,0)
    z = u
    for j in np.arange(nmx, 1, -1):
        dnx[j - 2] = j / z - 1 / (dnx[j - 1] + j / z)
    dnu = dnx[n - 1]
    z = v
    for j in np.arange(nmx, 1, -1):
        dnx[j - 2] = j / z - 1 / (dnx[j - 1] + j / z)
    dnv = dnx[n - 1]
    z = w
    for j in np.arange(nmx, 1, -1):
        dnx[j - 2] = j / z - 1 / (dnx[j - 1] + j / z)
    dnw = dnx[n - 1]

    nu = n + 0.5
    sv = np.sqrt(0.5 * np.pi * v)
    pv = sv * jv(nu, v)
    sw = np.sqrt(0.5 * np.pi * w)
    pw = sw * jv(nu, w)
    sy = np.sqrt(0.5 * np.pi * y)
    py = sy * jv(nu, y)
    p1y = np.append(np.sin(y), py[0:nmax1])
    chv = -sv * yv(nu, v)
    chw = -sw * yv(nu, w)
    chy = -sy * yv(nu, y)
    ch1y = np.append(np.cos(y), chy[0:nmax1])
    gsy = py - 1j * chy
    gs1y = p1y - 1j * ch1y

    uu = m * dnu - dnv
    vv = dnu / m - dnv
    fv = pv / chv
    fw = pw / chw
    ku1 = uu * fv / pw
    kv1 = vv * fv / pw
    ku2 = uu * (pw - chw * fv) + (pw / pv) / chv
    kv2 = vv * (pw - chw * fv) + (pw / pv) / chv
    dns1 = ku1 / ku2
    gns1 = kv1 / kv2

    dns = dns1 + dnw
    gns = gns1 + dnw
    a1 = dns / m2 + n / y
    b1 = m2 * gns + n / y

    an = (py * a1 - p1y) / (gsy * a1 - gs1y)
    bn = (py * b1 - p1y) / (gsy * b1 - gs1y)
    return np.array([an, bn])

def miecoated_ab2(m1: complex, m2: complex, x: float, y: float) -> np.ndarray:
    if (m1 * x).imag > 30:
        warnings.warn('Im(m1*x) > 30. This could result in bad solutions. See '
            'B&H p. 485. Compare the output of miecoated_ab1, miecoated_ab2, '
            'and miecoated_ab3 and be very careful if they give significantly '
            'different results.', InputDomainWarning)
    elif (m2 * x).imag > 30:
        warnings.warn('Im(m2*x) > 30. This could result in bad solutions. See '
            'B&H p. 485. Compare the output of miecoated_ab1, miecoated_ab2, '
            'and miecoated_ab3 and be very careful if they give significantly '
            'different results.', InputDomainWarning)
    elif (m2 * y).imag > 30:
        warnings.warn('Im(m2*y) > 30. This could result in bad solutions. See '
            'B&H p. 485. Compare the output of miecoated_ab1, miecoated_ab2, '
            'and miecoated_ab3 and be very careful if they give significantly '
            'different results.', InputDomainWarning)
    if x == y:
        raise ValueError('x == y, size parameters cannot be the same size.')
    elif x > y:
        raise ValueError('x > y, "inner" sphere larger than "outer"')

    m = m2 / m1
    nmax = np.int(np.round(2 + y + 4 * y**(1/3)))
    n = np.arange(1, nmax + 1).astype(int)
    nu = n + 0.5
    u = m1 * x
    v = m2 * x
    w = m2 * y
    su = np.sqrt(0.5 * np.pi * u)
    sv = np.sqrt(0.5 * np.pi * v)
    sw = np.sqrt(0.5 * np.pi * w)
    sy = np.sqrt(0.5 * np.pi * y)
    pu = su * jv(nu, u)
    py = sy * jv(nu, y)
    pv = sv * jv(nu, v)
    pw = sw * jv(nu, w)
    p1u = np.append(np.sin(u), pu[0:nmax - 1])
    p1y = np.append(np.sin(y), py[0:nmax - 1])
    p1v = np.append(np.sin(v), pv[0:nmax - 1])
    p1w = np.append(np.sin(w), pw[0:nmax - 1])
    ppv = p1v - n * pv / v
    ppw = p1w - n * pw / w
    ppy = p1y - n * py / y
    chv = -sv * yv(nu, v)
    chw = -sw * yv(nu, w)
    chy = -sy * yv(nu, y)
    ch1v = np.append(np.cos(v), chv[0:nmax - 1])
    ch1w = np.append(np.cos(w), chw[0:nmax - 1])
    ch1y = np.append(np.cos(y), chy[0:nmax - 1])
    gsy = py - 1j * chy
    gs1y = p1y - 1j * ch1y
    gspy = gs1y - n * gsy / y
    du = p1u / pu - n / u
    dv = p1v / pv - n / v
    dw = p1w / pw - n / w
    chpv = ch1v - n * chv / v
    chpw = ch1w - n * chw / w

    aan = pv * (m * du - dv) / (m * du * chv - chpv)
    bbn = pv * (m * dv - du) / (m * chpv - du * chv)

    a1 = ppw - aan * chpw
    a2 = pw - aan * chw
    b1 = ppw - bbn * chpw
    b2 = pw - bbn * chw
    an = (py * a1 - m2 * ppy * a2) / (gsy * a1 - m2 * gspy * a2)
    bn = (m2 * py * b1 - ppy * b2) / (m2 * gsy * b1 - gspy * b2)

    return np.array([an, bn])

def miecoated_ab3(m1: complex, m2: complex, x: float, y: float) -> np.ndarray:
    if (m1 * x).imag > 30:
        warnings.warn('Im(m1*x) > 30. This could result in bad solutions. See '
            'B&H p. 485. Compare the output of miecoated_ab1, miecoated_ab2, '
            'and miecoated_ab3 and be very careful if they give significantly '
            'different results.', InputDomainWarning)
    elif (m2 * x).imag > 30:
        warnings.warn('Im(m2*x) > 30. This could result in bad solutions. See '
            'B&H p. 485. Compare the output of miecoated_ab1, miecoated_ab2, '
            'and miecoated_ab3 and be very careful if they give significantly '
            'different results.', InputDomainWarning)
    elif (m2 * y).imag > 30:
        warnings.warn('Im(m2*y) > 30. This could result in bad solutions. See '
            'B&H p. 485. Compare the output of miecoated_ab1, miecoated_ab2, '
            'and miecoated_ab3 and be very careful if they give significantly '
            'different results.', InputDomainWarning)
    if x == y:
        raise ValueError('x == y, size parameters cannot be the same size.')
    elif x > y:
        raise ValueError('x > y, "inner" sphere larger than "outer"')

    m = m2 / m1
    nmax = np.int(np.round(2 + y + 4 * y**(1/3)))
    n = np.arange(1, nmax + 1).astype(int)
    nu = n + 0.5
    u = m1 * x
    v = m2 * x
    w = m2 * y
    su = np.sqrt(0.5 * np.pi * u)
    sv = np.sqrt(0.5 * np.pi * v)
    sw = np.sqrt(0.5 * np.pi * w)
    sy = np.sqrt(0.5 * np.pi * y)
    pu = su * jv(nu, u)
    py = sy * jv(nu, y)
    pv = sv * jv(nu, v)
    pw = sw * jv(nu, w)
    p1u = np.append(np.sin(u), pu[0:nmax - 1])
    p1y = np.append(np.sin(y), py[0:nmax - 1])
    p1v = np.append(np.sin(v), pv[0:nmax - 1])
    p1w = np.append(np.sin(w), pw[0:nmax - 1])
    ppv = p1v - n * pv / v
    ppw = p1w - n * pw / w
    ppy = p1y - n * py / y
    chv = -sv * yv(nu, v)
    chw = -sw * yv(nu, w)
    chy = -sy * yv(nu, y)
    ch1y = np.append(np.cos(y), chy[0:nmax - 1])
    gsy = py - 1j * chy
    gs1y = p1y - 1j * ch1y
    gspy = gs1y - n * gsy / y
    du = p1u / pu - n / u
    dv = p1v / pv - n / v
    dw = p1w / pw - n / w
    chpw = chw * dw - 1 / pw
    uu = m * du - dv
    vv = du / m - dv
    pvi = 1 / pv
    aaa = pv * uu / (chv * uu + pvi)
    bbb = pv * vv / (chv * vv + pvi)
    aa1 = ppw - aaa * chpw
    aa2 = pw - aaa * chw
    bb1 = ppw - bbb * chpw
    bb2 = pw - bbb * chw
    aa = (py * aa1 - m2 * ppy * aa2) / (gsy * aa1 - m2 * gspy * aa2)
    bb = (m2 * py * bb1 - ppy * bb2) / (m2 * gsy * bb1 - gspy * bb2)

    return np.array([aa, bb])

def miecoated(m1: complex, m2: complex, x: float, y: float, opt: int = 1) -> np.ndarray:
    if (m1 * x).imag > 30:
        warnings.warn('Im(m1*x) > 30. This could result in bad solutions. See '
            'B&H p. 485. Compare the output of miecoated_ab1, miecoated_ab2, '
            'and miecoated_ab3 and be very careful if they give significantly '
            'different results.', InputDomainWarning)
    elif (m2 * x).imag > 30:
        warnings.warn('Im(m2*x) > 30. This could result in bad solutions. See '
            'B&H p. 485. Compare the output of miecoated_ab1, miecoated_ab2, '
            'and miecoated_ab3 and be very careful if they give significantly '
            'different results.', InputDomainWarning)
    elif (m2 * y).imag > 30:
        warnings.warn('Im(m2*y) > 30. This could result in bad solutions. See '
            'B&H p. 485. Compare the output of miecoated_ab1, miecoated_ab2, '
            'and miecoated_ab3 and be very careful if they give significantly '
            'different results.', InputDomainWarning)
    if x == y:
        return mie(m1, y)
    elif x == 0:
        return mie(m2, y)
    elif m1 == m2:
        return mie(m1, y)
    elif x > 0:
        nmax = np.int(np.round(2 + y + 4 * y**(1/3)))
        n1 = nmax - 1
        n = np.arange(1, nmax + 1)
        cn = 2 * n + 1
        c1n = n * (n + 2) / (n + 1)
        c2n = cn / n / (n + 1)
        y2 = y**2
    elif x > y:
        raise ValueError('x > y, "inner" sphere larger than "outer"')
    if opt == 1:
        f = miecoated_ab1(m1, m2, x, y)
    elif opt == 2:
        f = miecoated_ab2(m1, m2, x, y)
    elif opt == 3:
        f = miecoated_ab3(m1, m2, x, y)
    else:
        raise ValueError('Invalid value passed as "opt" must be 1, 2, or 3')
    anp = f[0, :].real
    anpp = f[0, :].imag
    bnp = f[1, :].real
    bnpp = f[1, :].imag
    g1 = np.zeros((4, nmax))
    g1[0, 0:n1] = anp[1:nmax]
    g1[1, 0:n1] = anpp[1:nmax]
    g1[2, 0:n1] = bnp[1:nmax]
    g1[3, 0:n1] = bnpp[1:nmax]
    dn = cn * (anp + bnp)
    q = np.sum(dn)
    qext = 2 * q / y2
    en = cn * (anp**2 + anpp**2 + bnp**2 + bnpp**2)
    q = np.sum(en)
    qsca = 2 * q / y2
    qabs = qext - qsca
    fn = (f[0,:] - f[1,:]) * cn
    gn = np.power(-1, n)
    f = np.vstack((f, fn * gn))
    q = np.sum(f[2,:])
    qb = (q * q.conjugate() / y2).real
    asy1 = c1n * (anp * g1[0, :] + anpp * g1[1, :] + bnp * g1[2, :] + bnpp * g1[3,:])
    asy2 = c2n * (anp * bnp + anpp * bnpp)
    asy = 4 / y2 * np.sum(asy1 + asy2) / qsca
    qratio = qb / qsca
    return np.array([qext, qsca, qabs, qb, asy, qratio])

def miecoated_S12(m1, m2, x, y, u):
    nmax = np.int(np.round(2 + y + 4 * y**(1/3)))
    ab = miecoated_ab1(m1, m2, x, y)
    an = ab[0, :]
    bn = ab[1, :]
    pt = mie_pt(u, nmax)
    pin = pt[0,:]
    tin = pt[1,:]
    n = np.array(range(1, nmax + 1))
    n2 = (2*n+1) / (n*(n+1))
    pin = n2 * pin
    tin = n2 * tin
    S1 = (an @ pin.transpose()) + (bn @ tin.transpose())
    S2 = (an @ tin.transpose()) + (bn @ pin.transpose())
    return np.array([S1, S2])
