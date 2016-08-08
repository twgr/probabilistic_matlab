clear all
close all

this_loc = fileparts(mfilename('fullpath'));
addpath([this_loc filesep() 'example_inference_calls']);

run_branching;
run_gaussian;
run_hmm;
run_nlss;
run_lgss;