"""
This script constructs the diffraction spike information.
"""

import asdf
import numpy as np
from astropy.io import fits
from psfsim.polychrom import PolychromaticPSF
from scipy.interpolate import RegularGridInterpolator

# Filter information
romanfilters = "WFHJYZRK"
wl1 = [0.93, 1.68, 1.38, 1.13, 0.93, 0.76, 0.48, 1.95]
wl2 = [2.00, 2.00, 1.77, 1.45, 1.19, 0.98, 0.76, 2.30]

# Number of filters and number of SCAs to use
nf = len(romanfilters)
nsca = 18
nxpos = 3  # draw on nxpos x nxpos grid
dsp = 256  # number of pixels off the edge to draw a star
rmax = 200  # maximum radius in native pixels
ov = 4  # PSF oversampling
nt = 360  # number of angles
wedge = np.pi / 180 * 6.5  # half-opening angle to consider spike

# thresholds (highest to lowest, as a fraction of the total flux)
thresh = [1.0e-4, 1.0e-5, 1.0e-6]

# get coordinates
_s = np.linspace(-1, 1, nxpos) * (2044 + dsp) + 2043.5
x, y = np.meshgrid(_s, _s)
x = x.ravel()
y = y.ravel()
print(x)

angles = np.zeros((8, 18, 5, 12), dtype=np.float32)
lengths = np.zeros((8, 18, 5, 12, 3), dtype=np.float32)
widths = np.zeros((8, 18, 5, 12, 3), dtype=np.float32)

for ifilt in range(nf)[::-1]:
    for sca in range(1, nsca + 1):
        # positions on each filter
        for j in range(0, nxpos**2, 2):
            psf = PolychromaticPSF(
                sca, x[j], y[j], np.array([wl1[ifilt], wl2[ifilt]]), frame="science"
            ).compute_poly_psf(
                postage_stamp_size=2 * (rmax + 1),
                ovsamp=ov,
                use_filter=romanfilters[ifilt],
                reflect=False,
                req_in_bandpass=False,
                centerpix=False,
                cycle=10,
            )
            print(romanfilters[ifilt], sca, j, np.sum(psf))

            _r, _theta = np.meshgrid(np.array(range(ov * rmax)), np.array(range(nt)) * 2 * np.pi / nt)
            _x = (rmax + 1) * ov - 0.5 + _r * np.cos(_theta)
            _y = (rmax + 1) * ov - 0.5 + _r * np.sin(_theta)
            polarpsf = RegularGridInterpolator(
                (np.array(range(2 * (rmax + 1) * ov)), np.array(range(2 * (rmax + 1) * ov))),
                psf,
                method="linear",
            )(np.vstack((_y.ravel(), _x.ravel())).T).reshape(np.shape(_r))

            binpsf = np.sum(polarpsf[:, ov * 16 :], axis=1)
            peaks = np.where(
                np.logical_and(
                    np.logical_and(
                        np.logical_and(binpsf > np.roll(binpsf, 1), binpsf >= np.roll(binpsf, -1)),
                        np.cos(12 * _theta[:, 0]) < 0.1,
                    ),
                    binpsf >= 0.4 * np.amax(binpsf),
                )
            )[0]

            if len(peaks) < 12:
                peaks = np.where(
                    np.logical_and(
                        np.logical_and(
                            np.logical_and(binpsf > np.roll(binpsf, 1), binpsf >= np.roll(binpsf, -1)),
                            np.cos(12 * _theta[:, 0]) < 0.1,
                        ),
                        binpsf >= 0.2 * np.amax(binpsf),
                    )
                )[0]
            if len(peaks) > 12:
                peaks = np.where(
                    np.logical_and(
                        np.logical_and(
                            np.logical_and(
                                np.logical_and(binpsf > np.roll(binpsf, 1), binpsf >= np.roll(binpsf, -1)),
                                np.logical_and(binpsf > np.roll(binpsf, 2), binpsf >= np.roll(binpsf, -2)),
                            ),
                            np.cos(12 * _theta[:, 0]) < 0.1,
                        ),
                        binpsf >= 0.2 * np.amax(binpsf),
                    )
                )[0]
            if len(peaks) != 12:
                print(peaks)
                print(binpsf[peaks])
            assert len(peaks) == 12  # should happen

            theta = np.zeros(12)
            spikelen = np.zeros((12, len(thresh)))
            spikewid = np.zeros((12, len(thresh)))
            for k in range(12):
                p0 = binpsf[peaks[k]]
                pp = binpsf[(peaks[k] + 1) % nt]
                pm = binpsf[peaks[k] - 1]
                theta[k] = (peaks[k] + (pp - pm) / 2.0 / (2 * p0 - pm - pp)) * 2 * np.pi / nt

                # now make another grid
                along = np.array(range(rmax))
                cross = np.array(range(-10, 10))
                _xr, _yr = np.meshgrid(along, cross)
                _x = (rmax + 1) * ov - 0.5 + ov * (_xr * np.cos(theta[k]) - _yr * np.sin(theta[k]))
                _y = (rmax + 1) * ov - 0.5 + ov * (_xr * np.sin(theta[k]) + _yr * np.cos(theta[k]))
                tiltpsf = (
                    RegularGridInterpolator(
                        (np.array(range(2 * (rmax + 1) * ov)), np.array(range(2 * (rmax + 1) * ov))),
                        psf,
                        method="linear",
                    )(np.vstack((_y.ravel(), _x.ravel())).T).reshape(np.shape(_xr))
                    * ov**2
                )

                tiltpsf_trunc = np.where(np.abs(np.arctan2(_yr, _xr)) > wedge, 0.0, tiltpsf)
                tiltpsf_trunc = np.maximum(tiltpsf_trunc, tiltpsf_trunc[::-1, :])
                for b in range(len(thresh)):
                    spikelen[k, b] = np.amax(np.where(np.amax(tiltpsf, axis=0) > thresh[b])[0])
                    spikewid[k, b] = np.amax(
                        np.where(np.amax(tiltpsf_trunc, axis=1) > thresh[b])[0]
                    ) - np.amax(cross)

            print(np.round(theta * 180 / np.pi, 1))
            print(spikelen.T)
            print(spikewid.T)

            angles[ifilt, sca - 1, j // 2, :] = theta
            lengths[ifilt, sca - 1, j // 2, :, :] = spikelen
            widths[ifilt, sca - 1, j // 2, :, :] = spikewid


# Save information
asdf.AsdfFile(
    {
        "filters": romanfilters,
        "dsp": dsp,
        "wedge": wedge,
        "thresh": np.array(thresh),
        "angles": angles,
        "lengths": lengths.astype(np.float16),
        "widths": widths.astype(np.float16),
    }
).write_to("spikedata.asdf")

fits.HDUList(
    [
        fits.PrimaryHDU(angles.reshape(-1, 60)),
        fits.ImageHDU(lengths.reshape(-1, 180)),
        fits.ImageHDU(widths.reshape(-1, 180)),
    ]
).writeto("spikedata.fits", overwrite=True)
