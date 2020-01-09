*[Reply to this [Mmm.SE question](https://mathematica.stackexchange.com/questions/212167/how-much-can-we-trust-nintegrate-when-nintegrateeincr-is-shown)]*

## Introduction

I'll give a generic answer to a generic question (and mostly ignore Cuba).

What is supposed to happen in the global-adaptive strategy is that each recursive refinement, which subdivides the interval with the largest integration error, should reduce the error estimate.  But it doesn't always happen. Why not? Let's consider a few common, single-variable cases:

- The integrand is analytic or at least sufficiently smooth (the derivatives exist and are continuous to a higher order that the integration rule).
- The integrand has a singularity in the interval of integration.
- The integral is zero or near zero.
- The integrand is highly oscillatory.

Before we get started, let me mention that the integral and error estimate are based on the evaluation of the integrand $f(x)$ at a discrete set of nodes $x_1,\dots,x_n$.  Note that
$$
f(x) + (x-x_1)\cdots(x-x_n)\,g(x)
$$
has the same function values at the nodes as $f$ and therefore has the same integral and error estimate as $f$.  Since $g(x)$ is arbitrary, the actual error (not the computed error estimate) can be arbitrary.  In practice it is rare that you get this sort of aliasing problem except for a few of the subdivisions, since the nodes change when a subinterval is divided. (The bigger danger here is that the error is falsely estimated as close to zero and the subinterval is never divided.  But that has nothing to do with the `NIntegrate::eincr` warning.)



## Three cases


### Smooth integrand

In interpolatory function approximation, which is used in most integration rules,
there is a pre-convergent phase followed by a convergent phase.  

<sup> Example.
A three-point interpolatory rule (using a quadratic polynomial interpolant) cannot make a good approximation to a function that has five local extrema (or inflection points).  Either you have to use a higher-order rule (use more points) or subdivide the interval so that in each subinterval the function has at most one extremum; and if an inflection point cannot be avoided, make the interval very small so that the bad approximation will have little effect. (Those of you who know Simpson's rule know that because of its symmetry when integrated, we get a superconvergent rule with degree of precision 3 that can integrate cubic inflection points exactly.  But an asymmetric 3-point rule won't have this superconvergence, so bear with me.)  The main point is that the global-adaptive method has to break up the interval into small enough pieces that the integration rule starts to make good approximations of the integrand or its integral. After that point, the method will enter the convergent phase.
</sup>

During the preconvergent phase, one can expect the error estimate to increase sometimes.  For a highly oscillatory integrand, it might increase often.  If it increases too many times (for `NIntegrate`), in such a case, raising `"MaxErrorIncreases"` would theoretically fix the problem (and practically would, too, if the number of increases were not so great as to take too much time).  Another approach is to increase the order of the interpolatory rule using the `"Points"` suboption (see the docs; beware that Newton-Cotes becomes less stable as the order increases); this can shorten the preconvergent phase. 


### Singularity

A singularity often results in a `NIntegrate::slwcon` warning, and sometimes in a `NIntegrate::eincr` one, too.  Singularities tend to yield larger error estimates. As an interval is subdivided, the effect of the singularity might be come more noticeable.  Another way of saying it is that singularities can extend the preconvergent phase and error increases can accumulate.


### Zero integrals

Convergence to the desired `PrecisionGoal` can be particularly troublesome for integrals that are equal to zero or very close to zero (compared to the amplitude of the function, the working precision, and the `PrecisionGoal`). Recall that a precision goal $p$ is met by the error being many times smaller $10^{-p}$ than the result.  Potentially round-off error will dominate the truncation error of the integration method, and the error will jitter about causing many error increases.  This last problem is handled by setting `AccuracyGoal` to a finite value $a$ that you are comfortable that an estimate less $10^{-a}$ should be accepted as close enough to zero to be considered zero and further increase of precision should not be sought.

### Oscillatory integrand

This is actually a special subcase of a smooth integrand. Oscillatory integrals usually require a special method, since they have infinitely many extrema, and you can't just break up the interval into pieces in which the number of extrema is somewhat less than the order of the integration rule. The preconvergent phase could be arbitrarily long, depending on how fast the amplitude decays (which could depend on a transformation), during which error increases may be numerous.


## Examples

I've picked some contrived, simple-ish examples that can be investigated analytically.  In some cases, I force a bad method to illustrate how `eincr` works.

### Smooth integrand, long preconvergent phase

    nn = 2406;
    ff[x_?NumericQ] := ChebyshevT[nn, x]^2;
    NIntegrate[ff[x], {x, -1, 1}, MaxRecursion -> 16, 
     Method -> {"GlobalAdaptive", Method -> "GaussKronrodRule"}]
    % - Integrate[Cos[nn t]^2 Sin[t], {t, 0, Pi}] (* exact error *)
    
> NIntegrate::slwcon: Numerical integration converging too slowly; suspect one of the following: ...highly oscillatory integrand....
    
> NIntegrate::eincr: ...NIntegrate obtained 0.9294557079749661 and 0.0013688527738452687 for the integral and error estimates.
    
    (*
      0.929456
      -0.0705442
    *)

Note the exact error is quite a bit larger than the error estimate.
This can be fixed with `"MaxErrorIncreases"`:
    
    NIntegrate[ff[x], {x, -1, 1}, MaxRecursion -> 16, 
     Method -> {"GlobalAdaptive", Method -> "GaussKronrodRule", 
       "MaxErrorIncreases" -> 1000}]
    % - Integrate[Cos[nn t]^2 Sin[t], {t, 0, Pi}]  (* exact error *)
    
> NIntegrate::slwcon: ....
    
    (*
      0.9999999565975315` 
      3.32068*10^-13
    *)

The problem can be handled by pre-subdividing the interval.  The error estimate is not computed for this, and so it does not count in the number of error increases.  It also does not count in the heuristics that lead to the `NIntegrate:slwcon` warning, and sometimes it prevents that message.
One way is to use `MinRecursion`:

    NIntegrate[ff[x], {x, -1, 1}, MinRecursion -> 10, MaxRecursion -> 16, 
     Method -> {"GlobalAdaptive", Method -> "GaussKronrodRule"}]
    % - Integrate[Cos[nn t]^2 Sin[t], {t, 0, Pi}]
    (*
    1.
    1.62093*10^-14
    *)

Another way is to subdivide the interval manually at, say, the zeros of `ff[x]`:

    NIntegrate[ff[x], 
     Evaluate@Flatten@{x, Join[{-1}, -Cos[Pi ( Range[nn] - 1/2)/nn], {1}]}, 
     Method -> {"GlobalAdaptive", Method -> "GaussKronrodRule"}]
    % - Integrate[Cos[nn t]^2 Sin[t], {t, 0, Pi}]
    (*
      0.9999999568134225` 
      -6.66134*10^-16
    *)


### Singularities

    integrand[x_?NumericQ] := Abs[Sin[x]]; (* prevents symbolic preprocessing *)
    NIntegrate[integrand[x], {x, 0, 800 Pi}, MaxRecursion -> 20]
    
    > NIntegrate::eincr: ...NIntegrate obtained 1600.083569906413 and 0.18856172376838135 for the integral and error estimates.
    
    (*  1600.08  *)

This can be handled by pre-subdividing the interval.  The error estimate is not computed for this, and so it does not count in the number of error increases.  (Note this approach can be used when there is a long preconvergent phase, too.)
    
    NIntegrate[integrand[x], {x, 0, 800 Pi},
     MinRecursion -> 10, MaxRecursion -> 20]
    (*  1600.  *)

    NIntegrate[integrand[x], Evaluate@Flatten@{x, Range[0, 800] Pi}]
    (*  1600.  *)


### Integral near zero


    NIntegrate[Sin[x], {x, 0, 2 Pi}, MaxRecursion -> 16]

> NIntegrate::slwcon: Numerical integration converging too slowly; suspect one of the following: ...value of the integration is 0....

> NIntegrate::eincr: ...NIntegrate obtained 4.631711680858075`*^-16 and 5.169593219712382`*^-16 for the integral and error estimates.

    (*  4.63171*10^-16  *)

### Smooth oscillatory integrand, nonzero integral

This example has an unexpected "coincidence," 
which took me off guard.  But it shows that when numerical methods go bad,
you should be open and separate what you actually *know* from what usually happens.
The interval comprises many, many periods of the integrand.
Predictably, during the first several subdivisions, the error bounces around.
Since the interval is so long, it takes many subdivisions during which too many
error increases occur.

    NIntegrate[Sin[x]^2, {x, 0, 2^12 Pi}, MaxRecursion -> 11, 
     Method -> "GaussKronrodRule"]
    % - 2^12 Pi/2. (* actual error *)
    
> NIntegrate::eincr: ...NIntegrate obtained 7328.140790022457` and 295.19159332164276` for the integral and error estimates.
    
    (*
      7328.14
      894.159
    *)

A slight change in the interval and we get a very good estimate of the integral automatically. This turns out to be a surprising accident. It computes the correct result with just one subdivision, which you can check by adding `MaxRecursion -> 1`.
    
    NIntegrate[Sin[x]^2, {x, 0, (2^12 - 1) Pi}, 
     Method -> "GaussKronrodRule"]
    % - (2^12 - 1) Pi/2. (* actual error *)
    (*
      6432.41
      -5.63887*10^-11
    *)

Change it slightly again, just by adding some preliminary subdivision, and we're in trouble again.

    NIntegrate[Sin[x]^2, {x, 0, (2^12 - 1) Pi}, MinRecursion -> 2, 
     Method -> "GaussKronrodRule"]
    % - (2^12 - 1) Pi/2. (* actual error *)

> NIntegrate::eincr: ...NIntegrate obtained 6377.785070697375 and 143.07090946442491 for the integral and error estimates.

    (*
      6377.79
      -54.6259
    *)

Just to illustrate that the Gauss-Kronrod rule I chose above is contrived to produce the `NIntegrate::eincr` problem: The Levin rule is the automatic rule chosen by `NIntegrate`, and it produces a good result, with much more than the eight digits of precision sought with the default `PrecisionGoal`.
    
    NIntegrate[Sin[x]^2, {x, 0, 2^12 Pi}, Method -> "LevinRule"]
    % - 2^12 Pi/2. (* actual error *)
    *&
     6433.98
     6.89033*10^-9
    *)

No doubt you're curious what is going on.  Here is the coincidence:  It happens that a (pretty useless) property of the Gauss-Kronrod rule is that it calculates the integral of $\sin^2 x$ over an interval of the form $[m \pi/2, n \pi/2]$, where $m$ and $n$ are positive integers such that $m+n$ is odd, exactly and estimates the error to be zero (we don't get zero exactly because of round-off error; but run it with `WorkingPrecision -> 16` and you do get zero).  In the second integral, the first subdivision creates two subintervals of that form, $[0, (2^12-1)\pi/2]$ and $[(2^12-1)\pi/2, 2^13 \pi/2]$,
and the integral is computed exactly over each of them.  In the third integral, the interval is subdivided twice before any integration is computed.  This bypasses the special form intervals, and the Gauss-Kronrod rule is no longer exact.  `NIntegrate` proceeds as in the first integral, and we get an `eincr` message.


## Visualization

`NIntegrate` comes with an undocumented tool, the option [`IntegrationMonitor`](https://mathematica.stackexchange.com/questions/26401/determining-which-rule-nintegrate-selects-automatically/96663#96663).  With it, we can see the error increases (and could even show which subdivisions cause each one, but I'll omit that).

The following is a function that uses `IntegrationMonitor` to plot the error after each subdivision.  You can use it to illustrate any of the examples above.  You can also use it on any `NIntegrate` command whose method or rule utilizes `IntegrationMonitor` (not all do).

    ClearAll[errorPlot];
    SetAttributes[errorPlot, HoldAll];
    errorPlot[nint_NIntegrate, plotopts___?OptionQ] := Block[{integral, errors},
      {integral, errors} = Reap[
        Hold[nint] /. 
         Hold[Verbatim[NIntegrate][args___]] :> NIntegrate[args,
           IntegrationMonitor :> (Sow[Total@Through[#@"Error"]] &)]
        ];
      errors = Flatten@errors/integral;
      ListLinePlot[errors // RealExponent, PlotRange -> All, 
       Frame -> True, FrameLabel -> {"subdivisions", "log error"},
       PlotLabel -> 
        Row[{UnitStep@Differences@errors // Total, " error increases"}],
       DataRange -> {0, Length@errors - 1},
       plotopts]
      ]

We'll apply it to the last example. First a shorter interval, so we can see and count the error increases easily. After around 16 subdivisions, the preconvergent phase ends. One can check that in this (very symmetric) integral, we have 16 equal subintervals.  To see the convergence phase, each of these has to be subdivided. After 32 subdivisions the error decreases dramatically.  It will do it again after 64 subdivisions.  (There is an obvious blip at subdivision 31, but I did not investigate it. You can use `IntegrationMonitor` to do so, if curious.)

    errorPlot[NIntegrate[Sin[x]^2, {x, 0, 2^5 Pi},
      Method -> "GaussKronrodRule"]]

[![enter image description here][1]][1]

In the main example, we see that we're still in the preconvergent phase after 800+ subdivisions.
    
    errorPlot[NIntegrate[Sin[x]^2, {x, 0, 2^12 Pi},
      Method -> "GaussKronrodRule", MaxRecursion -> 11]]

[![enter image description here][3]][3]

Going further, we see that the preconvergent phase ends after about 2000 subdivisions (or `2^11`), and a big leap in convergence happens after another 2000 steps. A second leap occurs after another 4000 steps.  (Please keep in mind this doubling of the number of subdivisions comes from the symmetry of the integral and is not at all typical.)
    
    errorPlot[NIntegrate[Sin[x]^2, {x, 0, 2^12 Pi},
      Method -> {"GlobalAdaptive", Method -> "GaussKronrodRule", 
        "MaxErrorIncreases" -> 1000}, MaxRecursion -> 16]]

[![enter image description here][4]][4]


  [1]: https://i.stack.imgur.com/umXv7.png
  [3]: https://i.stack.imgur.com/NuOng.png
  [4]: https://i.stack.imgur.com/rWTN2.png
