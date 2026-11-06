% Step 1: Compute the mediation tensor for one subject.
% Set roi (time points x parcels), epsilonT, epsilonI and epsilonD before Run.
% Optional: icTimeCourses (time points x 20) and seed (default 1).
% Set icTimeCourses = [] to calculate ICs using the included bmn_fastica.m.
assert(exist('roi','var') && exist('epsilonT','var') && ...
    exist('epsilonI','var') && exist('epsilonD','var'), ...
    'Set roi, epsilonT, epsilonI and epsilonD in the workspace before Run.');
if ~exist('icTimeCourses','var')
    icTimeCourses = [];
end
if ~exist('seed','var')
    seed = 1;
end

[tensor, info] = diaomedia(roi, epsilonT, epsilonI, epsilonD, icTimeCourses, seed);
