# README #

### What is this repo? ###

This is a toolbox for automated inference in Matlab.  Models are written using a template format that is expressive and general purpose, but which forces an inference-friendly structuring.  Inference is based on SMC/PMCMC methods, with the enforced structuring producing fast inference and efficient memory usage; exploiting vectorization, compression, and efficient sample storage.  

The main advantages of the package are:

* Allows state-of-the-art general purpose inference in the form of interacting particle Markov chain Monte Carlo (iPMCMC) without requiring expertise in inference on the part of the user.
* It can run large numbers of particles and iterations without suffering memory issues for many models.  For example, on a 3 dimensional Kalman filter with 50 time-steps in the state sequence, with 32GB of available RAM, one can run up to ~10e6 particles and store ~200e6 total samples.
* Vectorization gives fast performance for many models.
* Can incorporate arbitrary deterministic external Matlab code.
* Provides output in a common format with automatic output processing functions provided.

### How do I get set up? ###

* Paths are added automatically when calling infer.m or example_run.m

### Usage ###

* You need to write a model using the model_template.m format.
* Inference on the model is performed using the infer.m function.
* example_models gives numerous examples of how to define models, whilst example_inference_calls gives examples of how to call infer and on results post-processing.
* Post-processing of the outputs can be done using method functions of @stack_object, for example, result summaries and histogram plotting.  There are also a few functions in generic_toolbox that are useful but do not operate directly on the @stack_object.
* If in doubt about which inference algorithm / options to use, try first to run SMC with as many particles as fit in memory.  If this is not sufficient, switch to using iPMCMC, again using as many particles as possible.
* See the model_template_commented.m and infer.m for more information.

When I get the time, I intend to write a compiler to convert Matlab code into required the template format, such that the package can be used as a probabilistic programming system.  Nonetheless, the current template-based model defintion format is simple to use it is own right and allows specification of a very general class of models. The main things that are not currently supported are:

* Memoization (though I hope to implement this in the future).
* Infinite recursion.
* The resampling points must remain fixed, though in some cases the resampling itself can be negated.  This should not in practise restrict the models that can be coded more than fork-based SMC/PMCMC methods as a) SMC/PMCMC based inferences require the number state-sequence steps to be constant across all particles and b) one can always define likelihood functions that branch depending on the variables present if necessary.

### Supported inference algorithms ###

* Sequential Monte Carlo ('smc')
* Particle Gibbs ('pgibbs')
* Particle independent Metropolis Hastings ('pimh')
* Alternate move Particle Gibbs ('a_pgibbs')
* Arbitrary independent combinations of the above algorithms  ('independent_nodes')
* Interacting particle Markov Chain Monte Carlo ('ipmcmc')

### Further information & Contribution guidelines ###

This work was developed as part of / in complement to the paper 

Rainforth, T., Naesseth, C. A., Lindsten, F., Paige, B., van de Meent, J.-W., Doucet, A., & Wood, F. (2016). Interacting Particle Markov Chain Monte Carlo. In Proceedings of the 33rd International Conference on Machine Learning, JMLR: W&CP (Vol. 48).

which should be used as reference for the inference algorithms.  Please use the following citation if you use this package:

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