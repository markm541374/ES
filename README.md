<h1>Gaussian process global optimisation using entropy techniques</h1>

<h2>ISSUES</h2>

<p>boundary value problem. Currently only using functionswith truemin in the interior 99% of each dimension. If the min lies
on the boundary then the search will only get as close as DIRECT can get to the boundary with the number of points allocated.
Consider alternative search method (which), setting DIRECT to search to within eps (will be slow)</p>

<p>I need to switch to slice sample for draws from teh min as rejection is too slow</p>

<p>Currenly only constraining the diagonal elements of the hessian to be +ve, rather than a +ve definite constraint on H</p>