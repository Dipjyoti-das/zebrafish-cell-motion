# zebrafish-cell-motion

The fortran90 code produces a list of cell positions and velocities over time, which obey 'Vicsek' dynamics.
The 'Cells' are soft particles in 2D, which interact via a short-range adhesive-repulsive potential. The cells also
align their velocities in presence of angular noise.
The cells are inside a horse-shoe geometry resembling a zebrafish tailbud, with rigid and reflective boundaries.


The matlab code is for visualization of the cell flow. The output file of the fortran90 code is the input for the matlab code.
