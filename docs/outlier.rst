Outlier rejection
#################

This step flags pixels that are outliers in a mosaic. This might occur for several reasons --- a cosmic ray that got missed by the jump detection, a temporarily hot pixel, or a moving or transient astronomical object.

Overview
========

The outlier rejection step is carried out for each mosaic. Typically, there are ≈320 exposure images per square degree, so a 1.25 deg^2 mosaic has ≈400 exposure images. The current way the code is implemented, the science layer and boolean mask for the whole mosaic is read into memory (≈30 GB), and multithreading is used to mask multiple exposures in parallel. This minimizes disk access, and the memory footprint is well within the capabilities of modern HPC systems.

The main interface is the ``OutlierMap`` class, which can be built from a PyIMCOM configuration. The simplest option is:

.. code-block:: python

    OutlierMap("myconfig.json", max_workers=24, run_and_save="mask_cube.asdf")

which will build the combined map from the PyIMCOM configuration fuile ``myconfig.json``, using a maximum of 24 workers, and save the output mask to ``mask_cube.asdf``. Note that since the Level 2 files are updated in place (with the new ``tree["mask"]`` object and ``tree["processinfo"]["outlier_complete"]`` flag), you don't actually need to keep the cube. You can specify parameters via:

.. code-block:: python

    OutlierMap("myconfig.json", max_workers=24, run_and_save="mask_cube.asdf",
        mask_kwargs={"cut_c": 1.0, "star_rad": 25})

(See `Outlier masks`_ for the full list of options.)

Lower-level control over the individual sub-steps is possible; the following would have the same effect:

.. code-block:: python

    # build the map structure (including reading in data)
    omap = OutlierMap("myconfig.json")

    # set the number of workers
    omap.max_workers = 24

    # build the masks
    # you can provide kwargs here to customize the algorithm
    omap.build_masks()

    # save the results; update keyword says to update
    # the original files
    omap.save_data("mask_cube.asdf", update=True,
        verbose=True)

Detailed processes
==================

We now describe the algorithms that are used here. Everything is implemented in the methods of the ``roman_hlis_l2_driver.outliers.OutlierMap`` class.

(**Note:** This algorithm isn't perfect, and we will be adding options to customize it in the future.)

Initialization
--------------

The initialization function (``.__init__``) reads the configuration file, and extracts the list of input ASDF files matching the configuration. This list is stored as ``self.files`` and the number of files is ``self.n_obs``. The science data is loaded as a 3D numpy array ``self.image``, so that ``self.image[f,y,x]`` is pixel (x,y) in the f-th file.

Mask information is loaded from the mask file, which depends on the ``file_format`` indicated in the configuration file. The options are:

+-------------+------------------------------------------+
| Format      | Mask file                                |
+-------------+------------------------------------------+
| ``L2_2506`` | Replace ``.asdf`` extension with         |
|             | ``_mask.fits``.                          |
+-------------+------------------------------------------+

(*More options will be added in the future.*)

The images stored for later use are "blotted" with a 3x3 kernel to faciliate comparison of interpolated undersampled images: the science data are smoothed by::

     [  1/16   1/8    1/16  ]
     [  1/8    1/4    1/8   ]
     [  1/16   1/8    1/16  ]

which has a transfer function that has exact nulls at spatial frequencies of 0.5 cycles/pixel. The mask is expanded 3x3 as well.

Finally, the WCS is extracted from the ASDF GWCS, converted to a ``PyIMCOM_WCS`` object, and stored in a list. The WCS for the i-th image is stored in ``self.wcs[i]``.

Overlap matrix
--------------

The overlap matrix of the exposures is computed next. This is an n_obs x n_obs matrix: values of ``self.ovmat[i,j]`` range from 0 (no overlap between exposures in and j) to 1 (full overlap).

Outlier masks
-------------

The outlier masks themselves are computed using the ``outlier_mask`` method::

    outlier_mask(
        self,
        iobs,             # which exposure
        cut_c=4.0,        # here & below are settings
        cut_m=0.2,
        cut_g=0.2,
        dyn4=0.2,
        dyn8=0.1,
        star_flux=300.0,
        star_rad=20,
        mdscale=0.333333,
    ):

The algorithm first uses the ``overlap`` method to build a 6 x 4088 x 4088 image, where the layers refer to the interpolated overlapping images:

+----------+---------------------------------------+
| Layer    | Usage                                 |
+----------+---------------------------------------+
| 0        | Number of ummasked overlapping images |
+----------+---------------------------------------+
| 1        | Median of overlapping images          |
+----------+---------------------------------------+
| 2        | Maximum of overlapping images         |
+----------+---------------------------------------+
| 3        | Minimum of overlapping images         |
+----------+---------------------------------------+
| 4        | Norm of gradient of the medium map    |
+----------+---------------------------------------+
| 5        | The observation itself                |
+----------+---------------------------------------+

Many entries in layers 1--5 may be nan if there are no overlapping images or if the observation itself is masked. When considering a given observation, we exclude it from the statistics in layers 0--4.

We also estimate the noise level based on the interquartile range of the image (``iqr``).

The current algorithm works as follows:

* A pixel is flagged as an outlier if it is larger than the maximum of other images or smaller than the minimum of other images by an amount::

    cut_c * iqr + cut_m * {median, clipped at 0} + cut_g * {gradient}

* A pixel is "rescued" unless it is at least ``dyn4`` times the maximum of the median image in a 4-pixel radius and ``dyn8`` times the maximum of the median image in an 8-pixel radius.

* A pixel is "rescued" if there are only 0 or 1 overlapping images (not sufficient information to build an outlier mask --- shouldn't happen except near survey boundaries or unlucky combinations of missed observations).

* The mask bit is turned off if the pixel is already masked in the input image.

* The mask is grown 3x3 so that we don't miss, e.g., edges of a defect.

* Individual pixels are masked if the original *unsmoothed* image is at least ``mdscale`` times the median image and the pixel isn't already masked (this sometimes helps with unflagged cosmic rays, although the current implementation isn't perfect).

* A region of radius ``star_rad`` is rescued from around bright stars (with a pixel exceeding ``star_flux`` DN/s) --- this is so that those stars aren't going to be missing from the coadded images (even though we'll eventually want to mask them for object detection).

The resulting mask is stored as a Boolean 3D numpy array in ``self.outmask``, so that ``self.outmask[f,y,x]`` is True if pixel (x,y) in file f is flagged as an outlier.

Saving the mask cube
--------------------

The ``save_data`` method stores the data to an output ASDF file; the tree contains::

   N : number of observations (int)
   files : list of image files (list of str)
   overlap : the overlap matrix (np.ndarray of float)
   mask : the output mask (np.ndarray of bool)

If the ``update`` keyword is set to True, then the original files will be updated by calling the ``update_files`` method. This method will:

* add the ``mask`` field if it isn't already in the ASDF file (and if it is not present, it sets the 0x1 bit based on the pre-existing Level 1->2 mask);

* set the 0x2 bit based on the outlier rejection; and

* set ``tree["processinfo"]["outlier_complete"]`` to True.
