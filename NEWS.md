# merlin 0.0.2

`merlin` 0.0.2 is more robust and it includes several bug fixes that have been reported/discovered since its initial release.
Among others:

* `coef`, `vcov`, `logLik` methods are now implemented for models of class `merlin` and `mlsurv`;
* The `summary` method has been revamped, with new options to customise printouts;
* Model with a single parameter are fitted using Brent's method;
* Starting values for Gompertz models are improved to improve convergence; 
* Warnings are now thrown when models do not converge;
* Parameters that control the fitting process can now be passed through the `control` argument;
* A separate function named `mlsurv` has been implemented to fit a variety of survival models with `merlin` using a simplified interface.

# merlin 0.0.1

Initial release.
