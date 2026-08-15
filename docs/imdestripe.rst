Image destriping
################

This is the ``roman_hlis_l2_driver`` wrapper for the image destriping routines (``pyimcom.imdestripe``). For detailed descriptions of the destriping implementation, see the pyimcom documentation. This page describes how to use the driver.

Overview
========

The destriping works at the mosaic level, and takes a PyIMCOM configuration file. The interface is simple::

    from roman_hlis_l2_driver.destripe_interface.destripe import destripe_all_layers

    destripe_all_layers(cfg_file, verbose=True)

where ``cfg_file`` is the location of the configuration file.

**Important notes on compatibility**:

- The file stems in ``DSOUT`` and ``DSOBSFILE`` should be unique among all versions that are running at the same time.

  - Right now, only hacks are available to place the temporary observations on on-node storage.

- You will have a collision of the log files if you run two copies of imdestripe in the same directory.

- The directories in ``DSOOUT`` and ``DSOBSFILE`` will get cleared of certain file types, so don't use directories where you are storing important information!

(We are working on these.)

Configuration file
==================

There are several additional fields in the configuration file that control the destriping:

- Fields for controlling the algorithm (see ``pyimcom.imdestripe`` documentation for details on what these do; note "None" is a legal argument in some cases):

  - ``DSMODEL``: str, int

    Destriping model and number of ds_rows

  - ``CGMODEL``: str, int, float

    Conjugate-gradient method, maximum number of iterations, and convergence tolerance.

  - ``DSCOST``: str, float, float

    Cost function type, prior coefficient, and transition (used for the Huber function).

- Fields for controlling where temporary data are stored:

  - ``DSOUT``: str, str

    Stem and path for temporary storage of de-striped images.

  - ``DSOBSFILE``: str

    Stem for temporary storage of original (not destriped) images.

Some fields (``GAINDIR`` and ``DSNOISEFILE``) were used for testing and may be used in the future, but should not be used in the current interface.

Workflow details
================

The workflow in ``roman_hlis_l2_driver.destripe_interface.destripe.destripe_all_layers`` is as follows:

- First, the science layer is destriped. This involves a call to ``roman_hlis_l2_driver.destripe_interface.destripe.destripe_one_layer``. This in turn performs the functions:

  - Convert the input ASDF files to FITS files starting with the ``DSOUT`` stem. (We found that imdestripe was much faster when operating on FITS files.) This uses the ``roman_hlis_l2_driver.destripe_interface.destripe_setup.setup_all_files`` function; the GWCS is converted to a TAN-SIP approximation, which has acceptable error for destriping.

  - Checks the noise cube, and extracts the number of noise layers, ``n_noise_layer`` (this is the return value of ``destripe_one_layer``).

  - Clears files in the ``cfg["DSOUT"][0]`` directory (where ``cfg`` is the PyIMCOM JSON configuration).

  - Calls the main de-striping routine, ``pyimcom.imdestripe.main``.

  - Inserts the destriped image in the original ASDF file as ``tree[destripe_orig]``. The destriping counter is also set: ``tree["processinfo"]["destripe"] = 0``, and a completion flag is inserted, ``tree["processinfo"]["destripe_complete"]`` (initialized to False).
  
- Now the noise layers (in ``*_noise.asdf``) are destriped. This is done by addditional calls to ``destripe_one_layer`` (one call per noise layer, with the ``noiseid`` keyword being used to indicate which layer is being destriped). Differences from the science layer call are:

  - This time, a FITS file is written with the *sum* of the science and noise layers; then ``pyimcom.imdestripe.main`` is run, and the de-striped science image is subtracted. (This ensures that bright object masking, etc. occurs when run on the noise images as well.)

  - We check that the destriping counter ``tree["processinfo"]["destripe"]`` matches the noise layer index, and increment it at the end (this will raise an exception if there is a file-writing error instead of blindly continuing).

  - The ``*_noise.asdf`` image is updated with the de-striped noise image overwriting the original.

- When the last noise layer is run, ``destripe_one_layer`` recognizes this based the counter (``tree["processinfo"]["destripe"]``). It then moves the de-striped science image into the data field (``tree["roman"]["data"]``), and sets the completion flag ``tree["processinfo"]["destripe_complete"]`` to True.