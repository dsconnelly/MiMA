#!/bin/bash

singularity exec --overlay $overlay:ro --bind /usr,/etc $image ./compilescript.csh
