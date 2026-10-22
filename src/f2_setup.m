function Th = f2_setup(P)
%F2_SETUP  Initialise Field II and build the linear array.
%   Safe to call repeatedly in one MATLAB session: it tears down a previous
%   Field II state first, then re-initialises.  Never call field_end from a
%   script that later needs the mex files (it crashed the headless session).

if exist('field_init', 'file') ~= 2
    addpath(P.fieldIIPath);
end

% NB: field_end is deliberately never called.  On this build a field_end
% followed by a new field_init aborted the MATLAB process, and field_init(-1)
% on its own is sufficient to (re)initialise the library.
field_init(-1);
set_field('c', P.c);
set_field('fs', P.fs);
set_field('show_times', 0);

Th = xdc_linear_array(P.N, P.width, P.height, P.kerf, 1, 1, [0 0 P.elevFocus]);

fc  = P.pulseCarrier;
t   = (0:round(P.pulseHalfLen/fc*P.fs))' / P.fs;
sig = P.pulseSigma / fc;
imp = exp(-((t - 3/fc).^2) / (2*sig^2)) .* sin(2*pi*fc*(t - 3/fc));

xdc_impulse(Th, imp');          % Field II wants a row vector
xdc_excitation(Th, 1);          % delta excitation: one-way pulse = imp
xdc_apodization(Th, 1, ones(1, P.N));
xdc_apodization(Th, 2, ones(1, P.N));
xdc_center_focus(Th, [0 0 0]);
end
