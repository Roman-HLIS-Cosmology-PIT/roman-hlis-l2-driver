|badge1| |badge2|

.. |badge1| image:: https://codecov.io/gh/Roman-HLIS-Cosmology-PIT/roman-hlis-l2-driver/graph/badge.svg

.. |badge2| image:: https://github.com/Roman-HLIS-Cosmology-PIT/roman-hlis-l2-driver/actions/workflows/smoke-test.yml/badge.svg

Level 2 Processing Tools
########################

These are the tools to convert Roman "Level 2" images (individual calibrated exposures, but with no information from other exposures) into "Level 2.1" images (with additional corrections/flags baed on other images). There are also tools to generate information needed by the PSF fitting tools; these additional steps that do not correct/modify the files themselves produce "Level 2.2" products.

Steps available:

* Level 2.0 --> 2.1:

  * `Image destriping <docs/imdestripe.rst>`_.

  * `Outlier rejection <docs/outlier.rst>`_.

  * `Diffraction spike masking <docs/spike.rst>`_.

* Level 2.1 --> 2.2:

  * `Full field of view map information <docs/ffovmaps.rst>`_.

Planned additions:

* Persistence flagging.

Level 2.0-->2.1 conventions
===========================

The Level 2.0-->2.1 processing is carried out for each mosaic prior to running through PyIMCOM. It takes in a list of ASDF files, and these are updated in place. Steps add to the ``processinfo`` fields in the ASDF files to indicate that they have run successfully, e.g., they set::

    tree["processinfo"]["destripe_complete"] = True
    tree["processinfo"]["outlier_complete"] = True
    tree["processinfo"]["spike_complete"] = True

Note that in addition to the main data files ``*_{obsid:d}_{sca:d}.asdf``, the destriping also acts on the noise files, ``*_{obsid:d}_{sca:d}_noise.asdf``. The ``outlier_complete`` flag in the main data file is updated at the end of the step.

The Level 2.0-->2.1 processing also adds a top-level "mask" image, ``tree["mask"]``, of type uint8. The bits are as follows:

+-------------+-----------------------------------------+
| Bit         | Meaning                                 |
+-------------+-----------------------------------------+
| ``0x01``    | Masked based on Level 1-->2 processing  |
|             | (original data quality array remains in |
|             | ``tree["roman"]["dq"]``).               |
+-------------+-----------------------------------------+
| ``0x02``    | Flagged by outlier rejection step.      |
+-------------+-----------------------------------------+
| ``0x04``    | Flagged for persistence (*not yet       |
|             | implemented; reserved*)                 |
+-------------+-----------------------------------------+
| ``0x08``    | Flagged for diffraction spikes or       |
|             | other optical phenomena.                |
+-------------+-----------------------------------------+

The convention is that a "good" pixel is denoted by a 0, and a "bad" pixel (or at least not to be included in the PyIMCOM coadds) is denoted by a 1.
    
