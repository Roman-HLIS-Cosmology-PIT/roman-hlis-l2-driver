"""
Star catalog reader.

Classes
-------
StarCat
    Class to build an external star catalog with objects that will be near this observation.

"""

import re

import healpy
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from astropy.wcs import WCS
from pyimcom.config import Settings

from .. import pars


class StarCat:
    """
    Class to build an external star catalog with objects that will be near this observation.

    Parameters
    ----------
    img : roman-hlis-l2-driver.fullfovexport.FullFoVImage
        A full field of view image.
    starfile : str
        A formatted string, so that ``infile.format(*pars)`` returns
        the path to a star file whose properties are in the tuple ``pars``.
    star_src : dict, optional
        Parameters for extracting the star catalog.
    r : float, optional
        Minimum radius beyond the SCAs to search (in units of SCA size).
    verbose : bool, optional
        Whether to print lots of information during the initialization.

    Attributes
    ----------
    starfile : str
        A formatted string, so that ``infile.format(*pars)`` returns
        the path to a star file whose properties are in the tuple ``pars``.
    star_src : dict
        Parameters for extracting the star catalog.
    pos : np.ndarray of float
        A grid of points to search when looking for catalog objects.
        Format is that `pos[i, j]` corresponds to the `i`th point on the
        `j`th axis (`j`=0 for x, `j`=1 for y, `j`=2 for z).
    searchrad : float
        The search radius around the above points in radians.
    tiles : np.ndarray of int
        The tile numbers that should be searched when looking for stars.
    catalog : numpy.rec.recarray
        The star catalog.

    Methods
    -------
    __init__
        Constructor.
    _build_catalog
        Extract the star catalog from the files.
    to_parquet
        Save catalog as a parquet file.

    Notes
    -----
    The `star_src` dictionary may contain the following items (with defaults):

    * ``"GRID"``: The grid type (e.g., ``"healpix32"`` for HEALPix nested nside=32)

    * ``"CATTYPE"``: The catalog type. Current options are:

      * ``"ou2024star"``: Stars in the OpenUniverse 2024 parquet format. In this case,
        there must be a ``_flux_*.parquet`` file as well in the same directory.

    """

    def __init__(self, img, starfile, star_src=None, r=0.125, verbose=False):
        # save input information
        self.starfile = starfile
        self.star_src = star_src if star_src is not None else {}

        # default grid if not otherwise requested
        if "GRID" not in self.star_src:
            self.star_src["GRID"] = "healpix32"
        if "CATTYPE" not in self.star_src:
            self.star_src["CATTYPE"] = ""

        search_radius = Settings.pixscale_native * Settings.sca_nside * max(r, 0.125)

        # get a grid of points on the SCAs to project onto the sky
        s_ = np.linspace(0, pars.n_side - 1, 9)
        x_, y_ = np.meshgrid(s_, s_)
        x_ = x_.ravel()
        y_ = y_.ravel()

        # First, get the range of RA and Dec that we need.
        is_first = True
        for sca in range(1, 1 + pars.n_sca):
            if img.hdulist[f"WFI{sca:02d}"].header["HASWCS"]:
                thiswcs = WCS(img.hdulist[f"WFI{sca:02d}"].header)
                ra_, dec_ = thiswcs.all_pix2world(x_, y_, 0)
                if is_first:
                    ra = ra_
                    dec = dec_
                    is_first = False
                else:
                    ra = np.concatenate((ra, ra_))
                    dec = np.concatenate((dec, dec_))

        deg = np.pi / 180.0
        self.pos = np.stack(
            (np.cos(dec * deg) * np.cos(ra * deg), np.cos(dec * deg) * np.sin(ra * deg), np.sin(dec * deg))
        ).T
        self.searchrad = search_radius

        # Which catalog tiles do we need to read?
        tiles = []
        still_looking = True

        # HEALPix option
        m = re.match(r"healpix(\d+)$", self.star_src["GRID"].lower())
        if m and still_looking:
            nside = int(m.group(1))
            for i in range(len(dec)):
                tiles += healpy.query_disc(
                    nside, self.pos[i, :], self.searchrad + 1.0 / (8 * nside), inclusive=True, fact=8
                ).tolist()
            still_looking = False

        if still_looking:
            # Complain if none of the search options were successful.
            raise ValueError("Unrecognized grid format: " + self.star_src["GRID"])

        # Turn into a numpy array of which tiles we need
        tiles = list(set(tiles))
        tiles.sort()
        self.tiles = tiles = np.array(tiles)

        # Pull catalog from file
        self._build_catalog(verbose=verbose)

    def _build_catalog(self, verbose=False):
        """
        Extract the star catalog from the files.

        Parameters
        ----------
        verbose : bool, optional
            Whether to print lots of information during the initialization.

        Returns
        -------
        None

        Notes
        -----
        Extract the star catalog from the files.
        You can add parameters to the catalog. But in any case, it must have "ra" and "dec"
        columns.

        """

        pdlist = []  # pandas data frame list

        for i in list(self.tiles):
            thisfile = self.starfile.format(i)
            if verbose:
                print("Reading <<", thisfile)

            # what happens next depends on the type of catalog
            if self.star_src["CATTYPE"].lower() == "ou2024star":
                m = re.match(r"^(.+)_(\d+)\.parquet", thisfile)
                if not m:
                    raise ValueError("Can't parse file name: " + thisfile)
                thisfile2 = m.group(1) + "_flux_" + m.group(2) + ".parquet"
                if verbose:
                    print("    secondary file <<", thisfile2)

                # merge tables
                p1 = pq.read_table(thisfile).to_pandas()
                p2 = pq.read_table(thisfile2).to_pandas()
                pdlist.append(p1.merge(p2, on="id"))

            else:
                raise ValueError("Unrecognized catalog type: " + self.star_src["CATTYPE"])

        # convert to record array
        self.catalog = pd.concat(pdlist).to_records()

        if verbose:
            print(f"Original catalog: rows = {self.catalog.shape[0]}, columns = {self.catalog.dtype.names}")

        # trim catalog basesd on the caps provided
        N = self.catalog.shape[0]
        catpos = np.zeros((N, 3))
        deg = np.pi / 180.0
        catpos[:, 0] = np.cos(self.catalog["dec"] * deg) * np.cos(self.catalog["ra"] * deg)
        catpos[:, 1] = np.cos(self.catalog["dec"] * deg) * np.sin(self.catalog["ra"] * deg)
        catpos[:, 2] = np.sin(self.catalog["dec"] * deg)
        keep = np.zeros((N,), dtype=bool)
        r2max = (2 * np.sin(self.searchrad / 2.0)) ** 2
        for i in range(np.shape(self.pos)[0]):
            r2 = np.sum((catpos - self.pos[i, :][None, :]) ** 2, axis=1)
            keep |= r2 < r2max
        self.catalog = self.catalog[np.where(keep)]

        if verbose:
            print(f"Trimmed catalog: rows = {self.catalog.shape[0]}, columns = {self.catalog.dtype.names}")

    def to_parquet(self, outfile, columns=None):
        """
        Saves the catalog as a parquet file.

        Parameters
        ----------
        outfile : str or str-like
            Where to write the catalog.
        columns : list, optional
            Which columns to save in the file?

        Returns
        -------
        None

        """

        cat = self.catalog
        if columns is not None:
            cat = cat[columns]
        pd.DataFrame(cat).to_parquet(outfile)
