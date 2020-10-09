Given the lower and upper limits of integration `x1` and `x2`, this routine returns arrays `x` and `w` of length `n`,
containing the abscissas and weights of the Gauss-Legendre `n`-point
[quadrature formula](https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature).

Usage
```javascript
import gauleg from 'gauss-legendre';

const { x, w } = gauleg(n, x1, x2);
```
