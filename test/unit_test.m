clear variables

inference_algorithms = {'pgibbs','a_pgibbs','pimh','smc','ipmcmc','independent_nodes'};

% Fixed options
n_particles = 1000;
n_iter = 100;
b_parallel = false;
M = 32;
P = 16;

common_options_fixed = {};
% Fixed local options
local_options_fixed.ipmcmc = {'M',M;
                              'P',P;
                              'b_parallel',b_parallel;
                              'n_particles',n_particles;
                              'n_iter',n_iter};
local_options_fixed.pgibbs = {};
local_options_fixed.pimh =   {'n_particles',n_particles;
                              'n_iter',n_iter};
local_options_fixed.a_pgibbs = {'n_particles',n_particles;
                                'n_iter',n_iter};
local_options_fixed.smc = {'n_particles',n_particles;
                           'n_iter',n_iter*M;
                           'b_parallel',b_parallel};
local_options_fixed.independent_nodes = {'n_particles',n_particles;
                                         'n_iter',n_iter;
                                         'b_parallel',b_parallel};
                                         
% Common options to enumerate over
common_options_enumerate = {'resample_method',{'residual','stratified','multinomial','systematic'};
                            'b_compress',{true, false}};
% Local options to enumerate over
local_options_enumerate.ipmcmc = {'b_Rao_Black',{true,false}};
local_options_enumerate.pgibbs = {'b_Rao_Black',{true,false};
                                  {'n_particles','n_iter'},{{n_particles,n_iter},{n_particles*M,n_iter},{n_particles,n_iter*M}}};
local_options_enumerate.pimh = {'b_Rao_Black',{true,false}};
local_options_enumerate.a_pgibbs = {'b_Rao_Black',{true,false}};
local_options_enumerate.smc = {'n_islands',{1,M};
                               'prop_sub_sample',{1,1/n_particles}};
local_options_enumerate.independent_nodes = {'b_Rao_Black',{true,false};
                                             'Ms',{[32,0,0],[0,32,0],[0,0,32],[16,0,16]}};
                                         
%% Enumerate all the options sets

for n=1:numel(inference_algorithms)
    alg = inference_algorithms{n};
    fixed_opt = [common_options_fixed;local_options_fixed.(alg)]';
    base = {alg,fixed_opt{:}};
    enum_options = [common_options_enumerate;local_options_enumerate.(alg)];
    names_enum = enum_options(:,1);
    b_multi = cellfun(@iscell,names_enum);    
    values_unenum = enum_options(:,2);
    values_enum = enumerate_cells(values_unenum{:});
    options.(alg) = cell(numel(values_enum),1);
    for m=1:numel(values_enum)
        opts = base;
        for o = 1:numel(names_enum)
            if b_multi(o)
                for so = 1:numel(names_enum{o})
                    opts = [opts,names_enum{o}{so},values_enum{m}{o}(so)];
                end
            else
                opts = [opts,names_enum{o},values_enum{m}(o)];
            end
        end
        options.(alg){m} = opts;
    end
end
    
%% Enumerate tests

for n=1:numel(inference_algorithms)
    alg = inference_algorithms{n};
    [branching_KL.(alg),gaussian_KL.(alg),hmm_distance.(alg)] = deal(NaN(numel(options.(alg)),1));
    for m=1:numel(options.(alg))
        disp([options.(alg){m}]);
        tTest = tic;
        [branching_KL.(alg)(m),gaussian_KL.(alg)(m),hmm_distance.(alg)(m)] = run_tests_for_given_options(options.(alg){m});
        timeTest = toc(tTest);        
        disp(['Results ' num2str([branching_KL.(alg)(m),gaussian_KL.(alg)(m),hmm_distance.(alg)(m)])]);
        disp(['Time taken ' num2str(timeTest) ' seconds']);
    end
end