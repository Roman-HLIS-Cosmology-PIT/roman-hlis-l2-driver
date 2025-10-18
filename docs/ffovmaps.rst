Making full field of view images
################################

Full field of view images are one of the early data products in the Level 2 workflow.

Quick start guide
=================

The shortest way to build a full field of view image from Level 2 files is with the ``FullFoVImage`` tool::

  from roman_hlis_l2_driver.fullfovexport.fullfov import FullFoVImage

  # replace "mydir" with the location of your L2 files
  FullFoVImage("mydir/L2_H158_13906_{:d}.asdf",
        maskfile = "mydir/L2_H158_13906_{:d}_mask.fits"
        ).to_file("ffov_13906.fits", overwrite=True)
  # writes the output to "ffov_13906.fits"

The path arguments are regular expressions, which are evaluated with the argument (SCA) being from 1 through 18.

Output file structure
=====================

The full field of view files have a defined structure because they will be handed off to other tools, particularly `Piff <https://github.com/rmjarvis/Piff>`_.

The file is a FITS file with 19 HDUs: a Primary HDU and 18 Image HDUs (one for each SCA). The format is designed to package the information that will be needed for PSF fitting in as compact a form as possible.

Primary HDU
-----------

The Primary HDU contains the following keywords:

+------------+-------------------------------------------------+
| Keyword    | Description                                     |
+============+=================================================+
|``SOFTBIAS``| Code in the image corresponding to sky          |
|            | background.                                     |
+------------+-------------------------------------------------+
|``SLOPEMIN``| Signal in DN_lin/s corresponding to code        |
|            |      1.                                         |
+------------+-------------------------------------------------+
|``SLOPEMAX``| Signal in DN_lin/s corresponding to code        |
|            | 65534.                                          |
+------------+-------------------------------------------------+
|``DSLOPE``  | Difference in signal in DN_lin/s corresponding  |
|            | to successive codes.                            |
+------------+-------------------------------------------------+
|``MJD``     | The modified Julian date                        |
|            | at the start of the exposure.                   |
+------------+-------------------------------------------------+
|``TSTART``  | The observation start time as a string.         |
+------------+-------------------------------------------------+

Image extension HDUs
--------------------

There are 18 image extension HDUs, with names ``WFI01`` ... ``WFI18``. Each one contains the image from the corresponding chip. The images are stored as 2D (4088x4088) uint16 arrays.

Image data
^^^^^^^^^^
For pixels with values between 1 and 65534 inclusive, the signal can be obtained as:

    slope (in DN_lin/s) = DSLOPE * (code - SOFTBIAS),

where DSLOPE and SOFTBIAS are from the above table. There are two special values:

* 0 : Indicates the pixel was masked (for reasons other than saturation).

* 65535 : Indicates the pixel was saturated before the end of the exposure (but not flagged for any other reason).

Note that the FITS convention stores images internally as int16 using the ``BSCALE=1`` and ``BZERO=32768`` keywords. The conventions above assume the data is "unpacked" by the reader. Since we allow a ``SOFTBIAS`` we didn't have to choose such a convention, but it is convenient for those of us accustomed to CCD surveys.

Keywords
^^^^^^^^

The following keywords *not* related to the image contents or the WCS model are included:

+------------+-------------------------------------------------+
| Keyword    | Description                                     |
+============+=================================================+
|``ISVALID`` | Is this image valid? (T/F) This could be F      |
|            | if the data from one of the SCAs is missing or  |
|            | affected by a data problem. Usually this will   |
|            | be T.                                           |
+------------+-------------------------------------------------+
|``HASMASK`` | Was a mask applied? (T/F) Usually this will be  |
|            | T.                                              |
+------------+-------------------------------------------------+
|``HASWCS``  | Is a WCS provided? (T/F) Usually this will be T.|
+------------+-------------------------------------------------+

WCS information
^^^^^^^^^^^^^^^

If ``HASWCS`` is True, then a WCS is provided.

The WCS is provided in as a composition of 3 transformations (with example physical effects that could be incorporated as described below)::

  # Science frame
  #   (x,y)_sci
  #     |             Pixel-level error map
  #     |             2D array of (delta x, delta y) for each pixel
  #     V             ex: pixel lithography errors, asymmetric cross-talk
  # Regular frame
  #   (x,y)_reg
  #     |             SIP coefficients
  #     |             Low-order polynomial on each chip
  #     V             ex: optical distortions, tilt of each chip
  # Linearized frame
  #   (x,y)_lin
  #     |             TAN projection
  #     |             Affine transformation -> projection -> Euler angles
  #     V             ex: observatory pointing, relativistic aberration 
  # Celestial frame
  #   (RA,Dec)

Of these, the TAN projection (part of the FITS WCS standard) and the SIP coefficients (in common usage but not part of the FITS WCS standard) have their full data stored in the header, as is done with the usual "TAN-SIP" approach to storing WCSs (see `here <https://aspbooks.org/custom/publications/paper/347-0491.html>`_ for details). The pixel-level error map is stored separately, and a hash code is used to indicate which one should be used.

The logic to this approach is that each exposure gets its own unique version of the TAN and SIP coefficients (the observatory pointing and relativistic aberration are different in every exposure, and we expect some amount of focus variation, etc.). The pixel-level error map requires much more information to store, but we expect it to be highly stable (e.g., lithography errors are literally built into the structure of the chip).

The TAN-SIP part of the projection is stored in the FITS header in the usual way for such a WCS. This includes the usual FITS TAN-SIP WCS keywords, including the projection types ``CTYPE1  = 'RA---TAN-SIP'`` and ``CTYPE2  = 'DEC--TAN-SIP'``.

Some additional keywords provided are:

+------------+-------------------------------------------------+
| Keyword    | Description                                     |
+============+=================================================+
|``MAXWCSER``| The maximum of the pixel-level error map (in    |
|            | units of pixels).                               |
+------------+-------------------------------------------------+
|``ERRMAP``  | A string of 1--16 alphanumeric characters       |
|            | indicating the pixel-level error map            |
|            | reference file to use. ``ERRMAP="NULL"`` is     |
|            | reserved and indicates that the pixel-level     |
|            | error map is zero (hence no file is needed).    |
+------------+-------------------------------------------------+

We can't have ``ERRMAP`` provide the absolute path to the reference file since these may be in different locations on different systems. However, on each platform, the driver *will* have to know where to find the files. The format of the error map files is described `below <#pixel-level-error-maps>`_.

**Warning**: A FITS reader (e.g., ``ds9``) expecting TAN-SIP convention will only read the TAN-SIP part of the data, and so the coordinates on the display do not include the pixel-level error map. This is not a limitation we can eliminate easily since a standard FITS reader won't know to go look up the calibration reference file containing the pixel-level error map. However, we expect the pixel-level error maps to be small (<<1 pixel, but possibly still important for weak lensing or precision astrometry), so a visual display with these errors is still useful even if it doesn't represent the full accuracy of the WCS.

Pixel-level error maps
======================

The pixel-level error maps are stored in a separate file (except if ``ERRMAP=="NULL"``). A pixel-level error map has a single Primary HDU, consisting of a 3D float32 array. The shape is ``NAXIS3>=2`` (number of error maps), ``NAXIS2=4088`` (number of rows), and ``NAXIS1=4088`` (number of columns).

Normally an error map would be stored with a file name such as ``f"RomanPixelLevelError_{errmapname:s}_WFI{sca:02d}.fits"``.

The following keywords are required in the PrimaryHDU:

+------------+-------------------------------------------------+
| Keyword    | Description                                     |
+============+=================================================+
|``ERRMAP``  | The string of 1--16 alphanumeric characters     |
|            | indicating which pixel-level error map          |
|            | reference file this is. (This is also part of   |
|            | the file naming convention, but it is useful to |
|            | have in the header to protect against a file    |
|            | name getting changed.) ``"NULL"`` is not        |
|            | allowed.                                        |
+------------+-------------------------------------------------+

The convention is that the offset between the regular frame position of pixel (i,j) (in Python) is related to the Primary HDU data via::

  # x_reg(i,j) - (origin+i) = data[0,j,i]
  # y_reg(i,j) - (origin+j) = data[1,j,i]

where ``origin`` is 0 or 1 depending on whether you are using C or Fortran convention for the regular frame.

Additional slices through the error map (``data[2,j,i]``, etc.) are not part of the WCS, but are permitted. We intend to use them in the future if we discover we need to represent higher moments of the pixel-to-pixel variation in the response function.
