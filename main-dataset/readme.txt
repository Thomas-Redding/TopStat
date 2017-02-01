


Mathematical summary

Let SO(3) be the Lie group of special orthogonal 3x3 real matrices (i.e., R such that R R^T = I and det R > 0). The bi-invariant metric between two points R, Q in SO(3) is

	d(R, Q) = arccos((tr Q R^T - 1) / 2).

Intuitively, d(R, Q) is the angle (in radians) of the rotation Q R^T. Consider the subgroup G = {I, J}, where

			-1		0		0
	J  =	0		1		0	.
			0		0		-1

Let SO(3) / G be the quotient space consisting of right cosets G R. This space inherits the metric as

	d(G R, G Q) = min(d(R, Q), d(J R, Q)).

The attached CSV file is a 151 x 151 symmetric matrix, describing the pairwise distances among 151 points in SO(3) / G. This information should be enough to compute the persistent homology of the points.



Geological summary

In geology, a 'fault' is a plane separating two blocks of rock that have moved relative to each other by sliding along that plane. The positions, orientations, and interconnections of faults throughout the Earth give insight into how the Earth has deformed, for example in accommodating tectonic plate motions. By 'give insight' I mean that such data can be incorporated into detailed physical-mathematical models of deformation.

The orientation of a fault involves two pieces of information. The first piece of information is a pole P, meaning a unit vector perpendicular to the plane. Notice that P is defined only up to sign. That is, -P is an equally good pole. The second piece of information is a 'vorticity vector' W, which is a unit vector, in the plane, perpendicular to the relative movement direction. This description specifies W only up to sign (+W or -W), but the sign is unambiguously determined from the observed fault according to a right-hand rule, the details of which I'll omit right now.

We don't want to analyze Ps and Ws separately, because they are perpendicular and hence correlated. Working within a fixed Cartesian coordinate system, we form a rotation matrix R with rows P, W, and P x W that completely describes the orientation of the fault. The matrix with rows -P, W, and -P x W is an equally good descriptor, and this matrix is J R. So it turns out that fault orientations are in one-to-one correspondence with SO(3) / G. By adopting this view, we can take advantage of a large body of literature called 'orientation statistics'.

In the attached CSV file, the 151 points in SO(3) / G are fault orientations measured by Sarah Titus in Cyprus.


