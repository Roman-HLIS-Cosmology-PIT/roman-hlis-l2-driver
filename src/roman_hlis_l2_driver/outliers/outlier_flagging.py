import gc
import json
import os
import sys
import warnings
from concurrent.futures import ThreadPoolExecutor
import pyarrow.parquet as pq
import galsim
import asdf
import numpy as np
import pyimcom
from astropy.io import fits
from astropy.coordinates import SkyCoord
from pyimcom.config import Settings as Stn
from pyimcom.utils.compareutils import get_overlap_matrix, map_sca2sca
from pyimcom.wcsutil import PyIMCOM_WCS
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import maximum_filter
from scipy.signal import convolve2d

from ..name_util import stem_l2


def _load_img_wcs(arg):
    """Simple utility to load WCSs in parallel."""

    (file, im) = arg

    print(file)
    sys.stdout.flush()
    with asdf.open(file, memmap=True) as a:
        im[:, :] = a["roman"]["data"]
        return PyIMCOM_WCS(a["roman"]["meta"]["wcs"], use_float32=True, niter=1)


def _load_mask(arg):
    """Simple utility to load WCSs in parallel."""

    (file, hdu, bits, im) = arg

    print(file)
    sys.stdout.flush()
    with fits.open(file) as f:
        im[:, :] = f[hdu].data & bits != 0


def _blot(arr, a=0.25):
    """
    Two-axis in-place blotting of an array to suppress aliasing.

    The kernel is (a, 1-2a, a) on each axis. Here a=0.25 leads to
    zero response at 0.5 cy/pix.

    Parameters
    ----------
    arr : np.ndarray
        The array to blot.
    a : float, optional
        The parameter of the blotting (controls where the zero-response frequencies are).

    Returns
    -------
    None

    """

    k1d = np.array([a, 1 - 2 * a, a])
    kernel = np.outer(k1d, k1d)
    arr[:, :] = convolve2d(arr, kernel, mode="same", boundary="symm")


def _blot_mask(arr):
    """
    Two-axis in-place 3x3 blotting of the mask.

    Parameters
    ----------
    arr : np.ndarray
        The array to blot.

    Returns
    -------
    None

    See Also
    --------
    _blot
        The "science image" version of this function.

    """

    arr2 = np.copy(arr)
    arr2[1:, :] |= arr[:-1, :]
    arr2[:-1, :] |= arr[1:, :]
    arr[:, :] = arr2
    arr2[:, 1:] |= arr[:, :-1]
    arr2[:, :-1] |= arr[:, 1:]
    arr[:, :] = arr2


def update_files(info, verbose=False):
    """
    Updates the mask information in the targeted files.

    An "l2_mask" array is added in place to the ASDF files.

    Parameters
    ----------
    info : dict
        The information needed for the update. The keys need to include:

        - ``files`` : list of str; the files to be updated.

        - ``mask`` : np.ndarray of bool; the additional pixels to mask

    verbose : bool, optional
        Talk to the output.

    Returns
    -------
    None

    Notes
    -----
    The mask array is as follows:

    - bit 0 (0x01): masked based on Level 1->2 processing

    - bit 1 (0x02): masked based on outlier detection with overlap images

    """

    N = len(info["files"])
    for j in range(N):
        fj = info["files"][j]
        if verbose:
            print(j, fj)
            sys.stdout.flush()
        with asdf.open(fj, lazy_load=False) as a:
            mytree = a.tree
        if "mask" not in mytree:
            mytree["mask"] = np.zeros((Stn.sca_nside, Stn.sca_nside), dtype=np.uint8)
        with fits.open(fj[:-5] + "_mask.fits") as f:
            mytree["mask"] &= ~np.uint8(3)
            mytree["mask"] |= np.where(f["MASK"].data > 0, 1, 0).astype(np.uint8)
        mytree["mask"] |= info["mask"][j, :, :].astype(np.uint8) << 1
        mytree["processinfo"]["outlier_complete"] = True
        asdf.AsdfFile(mytree).write_to(fj)


class OutlierMap:
    """
    Class for making a map of outliers in the mosaic.

    Parameters
    ----------
    cfg_file : str
        The configuration file.
    max_files : int, optional
        If specified, limits reading to the specified number of files.
    max_workers : int, optional
        If specified, gives the maximum number of workers for the thread pool.
    run_and_save : str, optional
        If given, runs the outlier map and saves to the indicated ASDF file.

    Attributes
    ----------
    cfg : pyimcom.config.Config
        The configuration object.
    n_obs : int
        The number of SCA images.
    files : list of str
        The input files.
    wcs : list of pyimcom.wcsutil.PyIMCOM_WCS
        The WCSs of the observations.
    ovmat : np.ndarray of float
        The overlap matrix, shape: (`n_obs`, `n_obs`).
    image : np.ndarray of float
        The SCA images as a cube (number, y, x).
    mask : np.ndarray of bool
        The input mask as a cube (number, y, x), expanded 3x3.
    mask_noblot : np.ndarray of bool
        The input mask as a cube (number, y, x), not expanded.
    max_workers : int
        Number of workers for thread parallelization.

    Methods
    -------
    __init__
        Constructor.
    overlap
        Builds the data on overlapping images for this image.
    outlier_mask
        Constructs an outlier mask.
    build_masks
        Builds the masks for this set of observations.
    save_data
        Save data for these observations as an ASDF file.

    """

    def __init__(self, cfg_file, max_files=None, max_workers=8, run_and_save=None):
        # config file
        self.cfg = pyimcom.config.Config(cfg_file)

        # get the input files
        # first get the file prefix and information
        with open(cfg_file, "r") as file:
            cfg = json.load(file)
        file_prefix = cfg["INDATA"][0]
        file_format = cfg["INDATA"][1]
        filter = Stn.RomanFilters[cfg["FILTER"]]
        star_catalog = cfg["STAR_CATALOG"]  # KL 

        # add tails as needed
        file_prefix += stem_l2(file_format, filter)

        fdir, fileprefix = os.path.split(file_prefix)
        n = len(fileprefix)
        numus = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "_"}
        self.files = []
        for f in sorted(os.listdir(fdir)):
            if f[:n] == fileprefix and f[-5:] == ".asdf" and all(c in numus for c in f[n:-5]):
                self.files.append(os.path.join(fdir, f))
            if len(self.files) == max_files:
                break
        self.n_obs = len(self.files)
        self.star_catalog = star_catalog  # KL

        print(self.n_obs)

        # image arrays
        self.image = np.zeros((self.n_obs, Stn.sca_nside, Stn.sca_nside), dtype=np.float32)

        # WCSs and overlap matrix
        with ThreadPoolExecutor(max_workers=max_workers) as e:
            self.wcs = list(
                e.map(_load_img_wcs, [(self.files[i], self.image[i, :, :]) for i in range(self.n_obs)])
            )
        gc.collect()
        print(self.image[:4, :4, :4])
        sys.stdout.flush()

        # mask
        if file_format == "L2_2506":
            tail = "_mask.fits"
            hdu = "MASK"
            bits = np.uint8(1)
        else:
            print("Please add the new format to outlier_flagging.")
            raise ValueError(f"unsupported format in outlier_flagging: {file_format}")
        self.maskfiles = [f[:-5] + tail for f in self.files]
        self.mask = np.zeros((self.n_obs, Stn.sca_nside, Stn.sca_nside), dtype=bool)
        with ThreadPoolExecutor(max_workers=max_workers) as e:
            e.map(_load_mask, [(self.maskfiles[i], hdu, bits, self.mask[i, :, :]) for i in range(self.n_obs)])
        gc.collect()

        print("Loaded data")
        sys.stdout.flush()

        # blotting
        self.mask_noblot = np.copy(self.mask)
        with ThreadPoolExecutor(max_workers=max_workers) as e:
            e.map(_blot, [self.image[i, :, :] for i in range(self.n_obs)])
        with ThreadPoolExecutor(max_workers=8) as e:
            e.map(_blot_mask, [self.mask[i, :, :] for i in range(self.n_obs)])

        self.ovmat = get_overlap_matrix(self.wcs, subsamp=18, verbose=True)

        self.max_workers = max_workers  # save the maximum number of workers to use

        if run_and_save is not None:
            self.build_masks()
            self.save_data(run_and_save, update=True, verbose=True)

    def overlap(self, iobs):
        """
        Builds the data on overlapping images for this image.

        Parameters
        ----------
        iobs : int
            The observation number we are considering.

        Returns
        -------
        np.ndarray
            A (6, nside, nside) array; the first axis indicates:

            - 0: the number of unmasked observations overlapping
            - 1: the median of these observations
            - 2: the maximum of these observations
            - 3: the minimum of these observations
            - 4: the gradient of the median map
            - 5: the observation itself (nan if masked)

        """

        # get the list of other ("j") observations
        jlist = list(np.where(self.ovmat[iobs, :] > 1e-3)[0])
        stack = np.full((1, Stn.sca_nside, Stn.sca_nside), np.nan, dtype=np.float32)
        ct = np.zeros((Stn.sca_nside, Stn.sca_nside), dtype=np.uint8)
        ns = 1

        # KL
        # filter star catalog to include only stars within this SCA's ra, dec limits
        if self.star_catalog is not None:
            with pq.ParquetFile(self.star_catalog) as pf:
                star_catalog_df = pf.read().to_pandas()
            ny, nx = Stn.sca_nside, Stn.sca_nside
            corners = np.array([[0,  0 ], [nx, 0 ], [0,  ny], [nx, ny]])
            sky = self.wcs[iobs].pixel_to_world(corners[:, 0], corners[:, 1])

            ra  = sky.ra.deg
            dec = sky.dec.deg

            ra_min,  ra_max  = ra.min(),  ra.max()
            dec_min, dec_max = dec.min(), dec.max()

            star_catalog_df = star_catalog_df[
                (star_catalog_df['ra'] >= ra_min) & (star_catalog_df['ra'] <= ra_max) &
                (star_catalog_df['dec'] >= dec_min) & (star_catalog_df['dec'] <= dec_max)
            ]


        # now build the set of other observations
        for jobs in jlist:
            if jobs == iobs:
                continue  # don't compare to yourself!
            x_target, y_target, is_in_ref = map_sca2sca(self.wcs[iobs], self.wcs[jobs], pad=0)
            coords = np.column_stack((y_target.ravel(), x_target.ravel()))
            interpolator = RegularGridInterpolator(
                (np.arange(Stn.sca_nside), np.arange(Stn.sca_nside)),
                self.image[jobs, :, :],
                method="linear",
                bounds_error=False,
                fill_value=0.0,
            )
            interp_data = interpolator(coords).reshape((Stn.sca_nside, Stn.sca_nside))
            
            interpolator = RegularGridInterpolator(
                (np.arange(Stn.sca_nside), np.arange(Stn.sca_nside)),
                1.0 - self.mask[jobs, :, :],
                method="linear",
                bounds_error=False,
                fill_value=0.0,
            )
            interp_good = interpolator(coords).reshape((Stn.sca_nside, Stn.sca_nside)) > 1e-9
            del interpolator
            ct += interp_good
            if np.any(ct > ns):
                ns = ns + 1
                stack = np.vstack(
                    (stack, np.full((1, Stn.sca_nside, Stn.sca_nside), np.nan, dtype=np.float32))
                )
            for k in range(ns):
                stack[k, :, :] = np.where(
                    np.logical_and(ct == k + 1, interp_good), interp_data, stack[k, :, :]
                )

        # KL
        if self.star_catalog is not None:
            for index, row in star_catalog_df.iterrows():
                x_star, y_star = self.wcs[iobs].world_to_pixel(SkyCoord(row['ra']*u.degree, row['dec']*u.degree))
                x_star = int(x_star)
                y_star = int(y_star)
                if x_star < 50 or x_star > Stn.sca_nside - 50 or y_star < 50 or y_star > Stn.sca_nside - 50:
                    continue
                cutout_iobs = self.image[iobs, y_star-25:y_star+25, x_star-25:x_star+25]
                cutout_interp = interp_data[y_star-25:y_star+25, x_star-25:x_star+25]
                
                # measure the moments using galsim hsm
                try:
                    hsm_iobs = galsim.hsm.FindAdaptiveMom(cutout_iobs)
                    hsm_interp = galsim.hsm.FindAdaptiveMom(cutout_interp)
                    star_catalog_df.loc[index, f'hsm_sigma_iobs_{iobs}'] = hsm_iobs.moments_sigma
                    star_catalog_df.loc[index, f'hsm_centroid_iobs_{iobs}'] = hsm_iobs.moments_centroid
                    star_catalog_df.loc[index, f'hsm_sigma_interp_{jobs}'] = hsm_interp.moments_sigma
                    star_catalog_df.loc[index, f'hsm_centroid_interp_{jobs}'] = hsm_interp.moments_centroid
                except Exception as e:
                    print(f"Error measuring moments for star at index {index}: {e}")


        del interp_data

        # get the median and maximum
        out = np.zeros((6, Stn.sca_nside, Stn.sca_nside), dtype=np.float32)
        out[0, :, :] = ct
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            out[1, :, :] = np.nanmedian(stack, axis=0)
            out[2, :, :] = np.nanmax(stack, axis=0)
            out[3, :, :] = np.nanmin(stack, axis=0)
        del stack

        # gradient computation
        medim = np.pad(out[1, :, :], 1, mode="edge")
        dmdx = np.maximum(np.abs(out[1, :, :] - medim[1:-1, 2:]), np.abs(out[1, :, :] - medim[1:-1, :-2]))
        dmdy = np.maximum(np.abs(out[1, :, :] - medim[2:, 1:-1]), np.abs(out[1, :, :] - medim[:-2, 1:-1]))
        out[4, :, :] = np.maximum(dmdx, dmdy)
        del medim, dmdx, dmdy

        # fill in the observation
        out[-1, :, :] = np.where(self.mask[iobs, :, :], np.float32(np.nan), self.image[iobs, :, :])

        # write star catalog to file
        # KL need a better solution than this bc this is one per sca...
        if self.star_catalog is not None:
            star_catalog_df.to_parquet(f'star_catalog_with_moments_iobs_{iobs}.parquet', index=False)

        return out

    def outlier_mask(
        self,
        iobs,
        cut_c=4.0,
        cut_m=0.2,
        cut_g=0.2,
        dyn4=0.2,
        dyn8=0.1,
        star_flux=300.0,
        star_rad=20,
        mdscale=0.333333,
    ):
        """
        Constructs an outlier mask.

        Parameters
        ----------
        iobs : int
            The SCA number to mask.
        cut_c, cut_m, cut_g : float, optional
            The value from median at which to cut as an outlier.
        dyn4, dyn8 : float, optional
            "Save" pixels whose outlier rate is less than this times the brightest pixel in the median
            blot in a 4- or 8-pixel radius.
        star_flux : float, optional
            "Saves" pixels around stars that are this bright (so we don't mask the inner
            diffraction spikes).
        star_rad : int, optional
            Radius around the star (in native pixels) for the cut above.
        mdscale : float, optional
            Cut if the image multiplied by this factor is above the median image by more than the noise.

        Returns
        -------
        mask : np.ndarray of bool
            The *new* pixels to mask.
        out : np.ndarray of float
            The data on overlap of this image.

        See Also
        --------
        overlap
            This computes `out` (see documentation for conventions).

        Notes
        -----
        The mask is constructed according to the criteria:

        * Pixel is not already masked.

        * The difference from the median of the other images should be at least the sum of:

          * `cut_c` times the interquartile range of the image

          * `cut_m` times the median image (clipped at 0)

          * `cut_g` times the gradient image

          and it should be at least `dyn4` / `dyn8` times the brightest pixel in a 4 / 8-pixel radius.

        At the end, it is grown by one pixel in every direction (including diagonal).

        Pixels are then "saved" if close to a bright star as determined by `star_flux` and `star_rad`
        (this prevents the stars from themselves being masked, and then we don't know they are there).

        There is also a single pixel level masking stage where the pixel is cut if it is above
        the smoothed median mask by at least `cut_c` interquartile ranges.

        """

        iqr = np.nanpercentile(self.image[iobs], 75) - np.nanpercentile(self.image[iobs], 25)
        out = self.overlap(iobs)  # the overlap cube from this image.

        diff = np.clip(
            np.maximum(self.image[iobs] - out[2, :, :], out[3, :, :] - self.image[iobs]), 0, None
        )  # diff from either max or min
        cutim = cut_c * iqr + cut_m * np.clip(out[1, :, :], 0.0, None) + cut_g * out[4, :, :]
        mask = diff > cutim
        del cutim

        # clipping based on 4 and 8-pixel radius
        R = 4
        s = np.arange(-R, R + 1)
        dx_, dy_ = np.meshgrid(s, s)
        mask &= diff > dyn4 * maximum_filter(out[1, :, :], footprint=dx_**2 + dy_**2 < R**2)
        R = 8
        s = np.arange(-R, R + 1)
        dx_, dy_ = np.meshgrid(s, s)
        mask &= diff > dyn8 * maximum_filter(out[1, :, :], footprint=dx_**2 + dy_**2 < R**2)
        del diff

        mask &= out[0, :, :] >= 2  # only mask if there are at least 2 others in constructing the median
        mask &= np.logical_not(self.mask[iobs])  # don't mask if the pixel is already masked

        # grow mask into a group of 3x3 pixels
        # (ensures that if this pixel was masked because of a problem with one of its neighbors that was not
        # previously masked, and then "grew" into this pixel via blotting, that pixel will be masked)
        _blot_mask(mask)

        # see if this pixel is surrounded by masked pixels
        # mask |= convolve(self.mask_noblot[iobs].astype(np.int8), np.ones((3,3), dtype=np.int8),
        #   mode="same", method="direct") >= surround

        # if too bright relative to median image
        with asdf.open(self.files[iobs], memmap=True) as a:
            mask |= np.logical_and(
                a["roman"]["data"] * mdscale > out[1, :, :] + cut_c * iqr, out[0, :, :] >= 2
            )

        mask &= np.logical_not(self.mask_noblot[iobs])  # don't mask if the pixel is already masked

        # now clear inner regions from bright stars
        s = np.arange(-star_rad, star_rad + 1)
        dx_, dy_ = np.meshgrid(s, s)
        mask &= maximum_filter(out[1, :, :] > star_flux, footprint=dx_**2 + dy_**2 <= star_rad**2) == 0

        print(iobs, self.files[iobs], iqr, np.count_nonzero(mask))

        return mask, out

    def build_masks(self, **kwargs):
        """
        Builds the masks for this set of observations.

        Parameters
        ----------
        **kwargs : various, optional
            Arguments to pass to the outlier mask function.

        Returns
        -------
        None

        See Also
        --------
        outlier_mask
            See documentation for the description of the optional keywords.

        """

        # which images to run
        ilist = np.arange(self.n_obs)

        # initialize the output
        self.outmask = np.zeros((len(ilist), Stn.sca_nside, Stn.sca_nside), dtype=bool)

        def _outlier_mask(di):
            """Wrapper for outlier_mask"""
            print("Running:", di)
            sys.stdout.flush()
            r, _ = self.outlier_mask(ilist[0] + di, **kwargs)
            self.outmask[di, :, :] = r

        # 1/3 as many workers since this is memory intensive
        with ThreadPoolExecutor(max_workers=(self.max_workers + 2) // 3) as e:
            e.map(_outlier_mask, list(range(len(ilist))))

    def save_data(self, outfile, update=False, verbose=False):
        """
        Save data for these observations as an ASDF file.

        The file contains the branches:

        * ``N`` : int; number of observations

        * ``files`` : list of str; the list of files

        * ``overlap`` : np.ndarray of float; 2D fractional overlap matrix

        * ``mask`` : np.ndarray of bool; 3D mask image

        Parameters
        ----------
        outfile : str
            The output file (must end in '.asdf').
        update : bool, optional
            If set to True, updates the original files.
        verbose : bool, optional
            Whether to talk while updating files.

        Returns
        -------
        None

        """

        # check we're writing to ASDF.
        if outfile[-5:] != ".asdf":
            raise ValueError("Incorrect file ending: only ASDF supported.")

        tree = {"N": self.n_obs, "files": self.files, "overlap": self.ovmat, "mask": self.outmask}
        asdf.AsdfFile(tree).write_to(outfile)

        if update:
            update_files(tree, verbose=verbose)


# this stuff will eventually be moved to another script
if __name__ == "__main__":
    u = OutlierMap(sys.argv[1], max_files=32, max_workers=12, run_and_save="M1.asdf")
    # fits.PrimaryHDU(np.where(u.outmask, 1, 0)).writeto("maskpile.fits", overwrite=True)

    # for i in [164]:
    #    mask, out = u.outlier_mask(i)
    #    mask = mask.astype(np.int8)
    #    fits.PrimaryHDU(out).writeto(f"im{i:d}.fits", overwrite=True)
    #    fits.PrimaryHDU(mask).writeto(f"mask{i:d}.fits", overwrite=True)
