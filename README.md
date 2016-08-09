# README #

### What is this repo? ###

This code base currently falls somewhere between a probabilistic program and an SMC/PMCMC-based inference toolbox.  Models are written using a common format that is both expressive and maintains a format in which inference is efficient, exploiting vectorization, compression and minimal memory overhead.  The core advantages of the package are

* The ability to run huge numbers of particles and iterations without suffering memory issues.  For example, on a 3 dimensional Kalman filter with 50 time-steps in the state sequence, with 32GB of availible RAM, one can run up to ~10e6 particles and store ~200e6 total samples.
* Vectorization allows for very fast performance (~20 times faster than Anglican on the above model).
* Matlab's debugger provides easy debugging for the models.
* Can incorporate arbitrary deterministic external Matlab code (this may sometimes require a small amount of wrapper code and may be less efficient than writing directly in the template).
* Provides output in a common format with automatic output processing functions provided.
* Allows state-of-the-art general purpose inference in the form of interacting particle Markov chain Monte Carlo (iPMCMC).

### How do I get set up? ###

* Paths are added automatically when calling infer.m or example_run.m

### Usage ###

* You need to write a model using the model_template.m format (or model_template_commented.m which is indentical but provides more information).
* Inference on the model is performed using the infer.m function.
* Postprocessing of the outputs can be done using method functions of @stack_object, for example, result summarizing and histogram plotting.  See also generic_toolbox.
* example_models gives numerous examples of how to define models, whilst example_inference_calls gives examples of how to call infer and on results post-processing.
* If in doubt about which inference algorithm / options to use, try first to run SMC with as many particles as fit in memory.  If this is sufficient great, otherwise switch to using iPMCMC, again using as many particles as possible.

### Supported inference algorithms ###

* Sequential Monte Carlo (SMC)
* Particle Gibbs
* Particle independent Metropolis Hastings (PIMH)
* Alternate move Particle Gibbs
* Arbitrary independent combinations of the above algorithms 
* Interacting particle Markov Chain Monte Carlo (iPMCMC)

### Further information & Contribution guidelines ###

This work was developed as part of the paper 
Rainforth, T., Naesseth, C. A., Lindsten, F., Paige, B., van de Meent, J.-W., Doucet, A., & Wood, F. (2016). Interacting Particle Markov Chain Monte Carlo. In Proceedings of the 33rd International Conference on Machine Learning, JMLR: W&CP (Vol. 48).
which should be used as reference for the inference algorithms.  I am currently in the process of providing academic documentation of the innovations behind the package (rather than the inference algorithms).  Until then, please use the following citation if you use this package:
@inproceedings{rainforth-icml-2016,
  title = {Interacting Particle {M}arkov Chain {M}onte {C}arlo},
  author = {Rainforth, Tom and Naesseth, Christian A and Lindsten, Fredrik and Paige, Brooks and van de Meent, Jan-Willem and Doucet, Arnaud and Wood, Frank},
  booktitle = {Proceedings of the 33rd International Conference on Machine Learning},
  series = {JMLR: W\&CP},
  volume = {48},
  year = {2016}
}
 
### Who do I talk to? ###

* Tom Rainforth: twgr@robots.ox.ac.uk