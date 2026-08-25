Persistence flagging
####################

This step flags pixels that are affected by persistence from previous exposures.

Overview
========
Since persistence is the retention of signal from previous exposures, this module must be able to find and retrieve previous exposures for each observation ASDF file in the L2.1 directory.

As the L2.1 directory is populated in batches that may not necessarily be in sequential or chronological order, the module must be agnostic towards the sorting, or lack thereof, of the directories it searches through.

Upon finding associated previous exposures that are within the expected persistence lifetime set by 'delta_t_prime_max' (default argument set to 1,200 seconds) and have pixels whose values meet or exceed the 'signal_threshold' criteria (default argument set to 20,000 DN), the module will set the 0x04 flag ON for that pixel so that it will be excluded from later steps of the pipeline.
