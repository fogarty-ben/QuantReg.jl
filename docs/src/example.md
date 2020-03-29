# Example

To further explain the use of this package, consider the following short example. We use a
9,800 observation subset of the [1996 US natality dataset used in Abrevaya (2006).]
(http://qed.econ.queensu.ca/jae/2006-v21.4/abrevaya/)

To load this dataset and run quantile regressions, we first load in the Julia CSV.jl package
and the QuantReg.jl package then load the dataset from file.

```
julia> using CSV, QuantReg
[ Info: Precompiling QuantReg [a0becc08-653f-40d2-91e7-721373d1053f]

julia> df = CSV.read("./bwght.csv")
9800×14 DataFrames.DataFrame
│ Row  │ birthweight │ boy   │ married │ black │ age   │ highschool │ somecollege │ college │
│      │ Int64       │ Int64 │ Int64   │ Int64 │ Int64 │ Int64      │ Int64       │ Int64   │ 
├──────┼─────────────┼───────┼─────────┼───────┼───────┼────────────┼─────────────┼─────────┼─
│ 1    │ 2926        │ 0     │ 1       │ 0     │ 25    │ 1          │ 0           │ 0       │ 
│ 2    │ 3595        │ 0     │ 0       │ 1     │ 17    │ 0          │ 0           │ 0       │ 
│ 9798 │ 3325        │ 1     │ 1       │ 1     │ 26    │ 1          │ 0           │ 0       │ 
│ 9799 │ 3232        │ 1     │ 0       │ 1     │ 21    │ 0          │ 1           │ 0       │ 
│ 9800 │ 2495        │ 0     │ 0       │ 0     │ 18    │ 0          │ 0           │ 0       │ 

│ prenone │ presecond │ prethird │ smoker │ cigsdaily │ weightgain │
│ Int64   │ Int64     │ Int64    │ Int64  │ Int64     │ Int64      │
│ ────────┼───────────┼──────────┼────────┼───────────┼────────────┤
│ 0       │ 0         │ 0        │ 0      │ 3         │ 22         │
│ 0       │ 0         │ 0        │ 1      │ 0         │ 44         │
│ 0       │ 0         │ 0        │ 1      │ 0         │ 99         │
│ 0       │ 0         │ 0        │ 1      │ 0         │ 40         │
│ 0       │ 0         │ 0        │ 0      │ 2         │ 25         │
```

Using the `rq` function, we the fit three quantile regression models at the 0.25th, 0.50th,
and 0.75 quantiles. We specify that the models be fit using the Barrodale-Roberts simplex
algorithm, and otherwise allow the program to choose sensible defaults (see
[QuantReg.QuantRegModel](@ref)).

```
julia> models = rq(@formula(birthweight ~ boy + married + black + age + highschool +
                            somecollege + college + prenone + presecond + prethird + smoker
                            + cigsdaily + weightgain), df, τ=[0.25:0.25:0.75;],
                   fitmethod="br")

birthweight ~ 1 + boy + married + black + age + highschool + somecollege + college + prenone
+ presecond + prethird + smoker + cigsdaily + weightgain, τ=0.25
───────────────────────────────────────────────────────────────────────────────────────────
             Coefficient  Std. Error           t      P(>|t|)   95% CI Lower   95% CI Upper
───────────────────────────────────────────────────────────────────────────────────────────
(Intercept)  2649.39       51.4538     51.4905    0.0             2548.53        2750.25
boy            95.6971     14.2654      6.70834   2.07738e-11       67.734        123.66
married        29.5964     19.5786      1.51167   0.13065           -8.78172       67.9746
black        -237.98       18.8469    -12.627     0.0             -274.923       -201.036
age            -0.927394    1.54661    -0.599631  0.548766          -3.95906        2.10428
highschool     92.5733     20.0475      4.6177    3.92964e-6        53.2761       131.87
somecollege   128.694      23.1254      5.56505   2.68988e-8        83.3634       174.025
college       102.463      28.3025      3.62028   0.000295771       46.9842       157.941
prenone      -157.604     171.45       -0.919242  0.357992        -493.681        178.473
presecond      -9.87082    21.7078     -0.454712  0.649326         -52.4227        32.681
prethird       15.253      32.98        0.462493  0.643738         -49.3946        79.9006
smoker        230.066      37.497       6.13558   8.8137e-10       156.564        303.568
cigsdaily      -2.19568     2.19156    -1.00188   0.316426          -6.49159        2.10023
weightgain      2.63786     0.322655    8.17549   4.44089e-16        2.00539        3.27033
───────────────────────────────────────────────────────────────────────────────────────────

Degrees of freedom: 9800 total; 9786 residual

birthweight ~ 1 + boy + married + black + age + highschool + somecollege + college + prenone
+ presecond + prethird + smoker + cigsdaily + weightgain, τ=0.5
───────────────────────────────────────────────────────────────────────────────────────────
             Coefficient  Std. Error           t      P(>|t|)   95% CI Lower   95% CI Upper
───────────────────────────────────────────────────────────────────────────────────────────
(Intercept)   2845.8       47.0243     60.5175    0.0             2753.62        2937.97
boy            121.298     12.2698      9.88595   0.0               97.247        145.35
married         14.1307    17.4178      0.81128   0.417225         -20.0117        48.2731
black         -227.1       15.9263    -14.2594    0.0             -258.319       -195.882
age              4.79392    1.34459     3.56535   0.000365106        2.15825        7.42959
highschool      62.8769    17.9189      3.50897   0.000451868       27.7521        98.0016
somecollege     95.2196    20.0219      4.75577   2.00515e-6        55.9725       134.467
college         76.7584    24.0546      3.191     0.00142229        29.6063       123.91
prenone        -87.0081    75.4319     -1.15346   0.248748        -234.87          60.8541
presecond       31.7744    19.3375      1.64315   0.100384          -6.13109       69.6799
prethird        17.4735    42.536       0.410792  0.681234         -65.9059       100.853
smoker         230.295     34.8011      6.61748   3.84524e-11      162.078        298.513
cigsdaily       -2.19744    2.04107    -1.07661   0.28168           -6.19837        1.80348
weightgain       2.83568    0.374805    7.56575   4.21885e-14        2.10099        3.57038
───────────────────────────────────────────────────────────────────────────────────────────

Degrees of freedom: 9800 total; 9786 residual

birthweight ~ 1 + boy + married + black + age + highschool + somecollege + college + prenone
+ presecond + prethird + smoker + cigsdaily + weightgain, τ=0.75
───────────────────────────────────────────────────────────────────────────────────────────
             Coefficient  Std. Error           t      P(>|t|)   95% CI Lower   95% CI Upper
───────────────────────────────────────────────────────────────────────────────────────────
(Intercept)   3170.49      47.8338     66.2814    0.0            3076.73         3264.26
boy            126.695     13.553       9.34808   0.0             100.128         153.262
married          7.32203   18.7681      0.390132  0.696448        -29.4673         44.1114
black         -208.254     18.1268    -11.4887    0.0            -243.787        -172.722
age              6.28814    1.44357     4.35595   1.3384e-5         3.45843         9.11784
highschool      36.4407    19.0143      1.91649   0.0553318        -0.831229       73.7126
somecollege     76.7797    22.2705      3.4476    0.000567972      33.1249        120.434
college         53.1864    25.0569      2.12263   0.0338101         4.06978       102.303
prenone       -161.627     44.68       -3.61744   0.00029903     -249.209         -74.0451
presecond       46.6102    20.2305      2.30396   0.0212458         6.95427        86.2661
prethird        -6.30508   35.184      -0.179203  0.857782        -75.2731         62.6629
smoker         154.847     34.4751      4.49158   7.1507e-6        87.2692        222.426
cigsdaily       -5.60339    1.78481    -3.13949   0.00169747       -9.10199        -2.10479
weightgain       3.98305    0.355627   11.2001    0.0               3.28595         4.68015
───────────────────────────────────────────────────────────────────────────────────────────

Degrees of freedom: 9800 total; 9786 residual
```

Above, we can see the results of the three models printed to console. These results plus
deeper information about each of the models is also stored in the `models` object. We can
directly access each of the models by indexing the `models` object by τ values. For example,
to access the median regression object, we would do the following:

```
julia> models[0.5]

birthweight ~ 1 + boy + married + black + age + highschool + somecollege + college + prenone +
presecond + prethird + smoker + cigsdaily + weightgain, τ=0.5
───────────────────────────────────────────────────────────────────────────────────────────
             Coefficient  Std. Error           t      P(>|t|)   95% CI Lower   95% CI Upper
───────────────────────────────────────────────────────────────────────────────────────────
(Intercept)   2845.8       47.0243     60.5175    0.0             2753.62        2937.97
boy            121.298     12.2698      9.88595   0.0               97.247        145.35
married         14.1307    17.4178      0.81128   0.417225         -20.0117        48.2731
black         -227.1       15.9263    -14.2594    0.0             -258.319       -195.882
age              4.79392    1.34459     3.56535   0.000365106        2.15825        7.42959
highschool      62.8769    17.9189      3.50897   0.000451868       27.7521        98.0016
somecollege     95.2196    20.0219      4.75577   2.00515e-6        55.9725       134.467
college         76.7584    24.0546      3.191     0.00142229        29.6063       123.91
prenone        -87.0081    75.4319     -1.15346   0.248748        -234.87          60.8541
presecond       31.7744    19.3375      1.64315   0.100384          -6.13109       69.6799
prethird        17.4735    42.536       0.410792  0.681234         -65.9059       100.853
smoker         230.295     34.8011      6.61748   3.84524e-11      162.078        298.513
cigsdaily       -2.19744    2.04107    -1.07661   0.28168           -6.19837        1.80348
weightgain       2.83568    0.374805    7.56575   4.21885e-14        2.10099        3.57038
───────────────────────────────────────────────────────────────────────────────────────────

Degrees of freedom: 9800 total; 9786 residual
```

Say we wanted to make some change to the model from above. We could do so directly from this
model by using the provided `QuantRegModel` convenience constructor. In the following code
block, we generate a new model where regression errors are assumed to be iid from the
existing model at the 0.50th quantile.

```
julia> altmodel = QuantReg.QuantRegModel(models[0.5], iid=true)

birthweight ~ 1 + boy + married + black + age + highschool + somecollege + college + prenone
+ presecond + prethird + smoker + cigsdaily + weightgain, τ=0.5
────────────────────────
             Coefficient
────────────────────────
(Intercept)   2845.8
boy            121.298
married         14.1307
black         -227.1
age              4.79392
highschool      62.8769
somecollege     95.2196
college         76.7584
prenone        -87.0081
presecond       31.7744
prethird        17.4735
smoker         230.295
cigsdaily       -2.19744
weightgain       2.83568
────────────────────────
```

Notice that this model, stored in `altmodel` still has its fit but does not have inference
computed; this is because changing the assumption about the distribution of the regression
errors only affects computing inference. To compute inference for `altmodel`, we run the
following command:

```
julia> QuantReg.compute_inf!(altmodel)

birthweight ~ 1 + boy + married + black + age + highschool + somecollege + college + prenone
+ presecond + prethird + smoker + cigsdaily + weightgain, τ=0.5
───────────────────────────────────────────────────────────────────────────────────────────
             Coefficient  Std. Error           t      P(>|t|)   95% CI Lower   95% CI Upper
───────────────────────────────────────────────────────────────────────────────────────────
(Intercept)   2845.8        45.6468    62.3438    0.0             2756.32        2935.27
boy            121.298      12.5634     9.65487   0.0               96.6714       145.925
married         14.1307     17.352      0.814353  0.415463         -19.8829        48.1443
black         -227.1        16.2947   -13.9371    0.0             -259.041       -195.159
age              4.79392     1.29373    3.70551   0.000212125        2.25795        7.3299
highschool      62.8769     17.5895     3.57469   0.000352329       28.3979        97.3558
somecollege     95.2196     20.4008     4.66744   3.09033e-6        55.2297       135.209
college         76.7584     24.1841     3.17392   0.00150859        29.3526       124.164
prenone        -87.0081     58.4498    -1.4886    0.136626        -201.582         27.5656
presecond       31.7744     18.4814     1.71926   0.0855979         -4.45295       68.0018
prethird        17.4735     39.3888     0.443615  0.657331         -59.7368        94.6837
smoker         230.295      33.8347     6.80649   1.05878e-11      163.972        296.618
cigsdaily       -2.19744     2.19289   -1.00208   0.316332          -6.49596        2.10108
weightgain       2.83568     0.30841    9.19452   0.0                2.23114        3.44023
───────────────────────────────────────────────────────────────────────────────────────────

Degrees of freedom: 9800 total; 9786 residual
```

Returning to the original model fit at the 0.50th quantile, we can also access information
about its fit and inference directly. Since we fit via the Barrodale-Roberts simplex
algorithm, the program computes and stores dual solutions to the quantile regression
linear program. We can access those as follows:

```
julia> models[0.5].fit.dual
9800-element Array{Float64,1}:
 0.0
 1.0
 1.0
 0.0
 0.0
 ⋮
 1.0
 0.0
 0.0
 0.0
 0.0
```

There are a number of other values that can be accessed via the fields of the `QuantRegFit`
object stored at `models[0.5].fit`. For a full list, see the [QuantRegFit](@ref) section.

Similarly, we could access the column vector of standard errors direction as follows:

```
julia> models[0.5].inf.σ
14-element Array{Float64,1}:
 47.02433161327634
 12.26976221196334
 17.417768956074895
 15.92630617837017
  1.3445878800214834
  ⋮
 19.3375090921255
 42.53601653109182
 34.801095529285575
  2.0410717521211215
  0.3748053760950519
```

There are also a number of other values that can be accessed via the fields of the
`QuantRegInf` object stored at `models[0.5].inf`. For a full list, see the
[QuantRegInf](@ref) section.
