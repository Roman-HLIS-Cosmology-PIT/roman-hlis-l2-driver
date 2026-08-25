Persistence flagging
####################

This step flags pixels that are affected by persistence from previous exposures.

Overview
========
Since persistence is the retention of signal from previous exposures, this module must be able to find and retrieve previous exposures for each observation ASDF file in the L2.1 directory.

As the L2.1 directory is populated in batches that may not necessarily be in sequential or chronological order, the module must be agnostic towards the sorting, or lack thereof, of the directories it searches through.

Upon finding associated previous exposures that are within the expected persistence lifetime set by 'delta_t_prime_max' (default kwarg set to 1,200 seconds) and have pixels whose values meet or exceed the 'signal_threshold' criteria (default kwarg set to 20,000 DN), the module will set the 0x04 flag ON for that pixel so that it will be excluded from later steps of the pipeline.

You can run this via:

.. code-block:: python

    from roman_hlis_l2_driver.persistence.main import run

    # Here, cfg is a json dictionary with an 'INDATA' key storing the path to the L2.1 directory
    # and an 'OBSFILE' key storing the path to the corresponding observation table FITS file
    # containing the exposure date and time of the observations presently in the L2.1 directory.
    # We set the update kwarg to True so that the module automatically sets flags
    # for pixels showing signs of persistence.
    persistence_mask, prev_obsids, current_obsid = run(cfg, update = True)

    # persistence_mask = the cumulative np.array of boolean values
    # prev_obsids = list of the obsids for the previous exposures used to build the mask
    # current_obsid = observation ID of the exposure/file being processed


The files are updated with the new flag in place.
