Diffraction spike rejection
###########################

This step flags pixels that are affected by far diffraction spikes from bright stars.

Overview
========

The diffraction spike rejection step is carried out for each mosaic; this way, we can flag spikes from stars that fall into a chip gap, but would be observed in another exposure. There are 2 sub-steps: building the catalog of stars and flagging the spikes.

You can run this via:

.. code-block:: python

    from roman_hlis_l2_driver.starcatalogs.main import spike_driver

    # here 0.1 is the threshold level to use in setting the mask width/length
    # and the Level 2.1 files are in the indicated location, e.g.,:
    #   /path_to_observations/prefix_1442_7.asdf
    # for obsid=1442, sca=7
    starcat, tot_nmask = spike_driver("/path_to_observations/prefix_{0:d}_{1:d}.asdf", 0.1)

    # starcat = np.recarray of the stars whose spikes were masked
    # tot_nmask = total number of pixels masked

The files are updated with the new flag in place.

Detailed processes
==================

(*coming soon*)
