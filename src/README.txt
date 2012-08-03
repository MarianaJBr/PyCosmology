********PYCOSMOLOGY**********************

DESCRIPTION

PyCosmology is a collection of routines which perform cosmology-specific tasks, especially in relation to simulation data. There are currently three main sub packages - sims, structure and shape lets. 

The sub package sims deals with statistical analyses on groups in cosmological simulations (specifically Gadget simulations). It contains routines for reading in simulation data, along with groups and subgroups, and then performing statistical tests such as finding velocity dispersions, axis ratio distributions and much much more…

The sub package shapelets uses the simulation data from sims and evaluates shapelet coefficients for the various groups therein - then performs statistical tests thereon. Note it is not working yet.

The sub package structure is more theoretical. It computes power spectra, mass variances, and mass functions for a range of cosmologies. 


HISTORY
v4.0		-	sims contains routines for:
			- importing simulation data
			- re-centering and rotating data into optimal co-ordinates
			- plotting the projected density of each group
			- plotting density profiles of each group
			- plotting velocity dispersions of each group
			- plotting axis ratio distributions
			- plotting centre of mass offsets
			- performing structure in shells analysis via ratio of objects in subgroups to total AND
			- by calculating density variances for varying angles in shells.

			- documentation improved heavily and added comments.

<  v3.9	-	sims and structure work but probably not optimally. Shapelets doesn't work at all. 