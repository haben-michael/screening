R scripts supporting the manuscript "The Effect of Screening for Publication Bias on the Outcomes of Meta-Analyses" by Haben Michael.

``sim.R`` implements the simulation described in the simulation section of the manuscript. 

``process.R`` generates the figure and table in the manuscript from the output of ``sim.R``.

``utils.R`` contains helper functions called by ``sim.R``.

``submit.sh`` is a shell script submitting ``sim.R`` to a job scheduler (slurm).
