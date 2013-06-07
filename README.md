BLAND
=====

BLAND provides recursive generic rings in C on top
of FLINT (http://flintlib.org).

Starting from various FLINT ground types, it allows constructing
complex numbers, fractions, residues, polynomials, and matrices.
In particular, it allows working with multivariate polynomials
constructed in a recursive fashion.

Some (outdated) documentation can currently be found at
http://sage.math.washington.edu/home/fredrik/flint/generics/

"BLAND" can be read as a recursive backronym for
"BLAND is a Library for Algebra over Nested Domains"
(the word "bland" was chosen because it is likely to describe
something that is "generic").

To install it in FLINT, clone the repository into a directory
on your filesystem, e.g. /home/username/bland

Now build FLINT in the usual way,
passing --extensions=/home/username/bland to FLINT's configure.

Bug reports, feature requests and other comments are welcome on the
FLINT mailing list: flint-devel@googlegroups.com

Author: Fredrik Johansson <fredrik.johansson@gmail.com>

