function P = params(root)
%PARAMS  Single source of truth for all simulation parameters.
%   P = PARAMS() builds the parameter struct; ROOT is the project directory.
%
%   Reference: G. Montaldo et al., "Coherent Plane-Wave Compounding for Very
%   High Frame Rate Ultrasonography and Transient Elastography", IEEE TUFFC
%   56(3):489-506, 2009.

if nargin < 1 || isempty(root)
    root = pwd;
end

P.root    = root;
P.dataDir = fullfile(root, 'data');
P.figDir  = fullfile(root, 'figures');
P.docDir  = fullfile(root, 'docs');

% ---------------------------------------------------------------- acoustics
P.c      = 1540;          % m/s
P.fs     = 30e6;          % Hz, simulation sampling rate
P.f0     = 4.5e6;         % Hz, nominal centre frequency
P.bwGoal = 0.60;          % -6 dB relative bandwidth (paper)

% ------------------------------------------------------------------- array
P.N        = 128;         % elements (paper: 128, 124 active)
P.pitch    = 0.33e-3;     % m
P.width    = 0.30e-3;     % m (kerf 0.03 mm)
P.height   = 5.0e-3;      % m
P.kerf     = 0.03e-3;     % m
P.elevFocus= 30e-3;       % m, elevation lens focus
P.Fnum     = 1.75;        % receive F-number
P.xe       = ((0:P.N-1) - (P.N-1)/2) * P.pitch;   % element centres [m]

% ------------------------------------------------ transmit pulse (calibrated)
% The two-way pulse is the transducer impulse response applied twice, so the
% single-way Gaussian must be centred above f0.  E0 sweep on a single scatterer
% at 30 mm: carrier 4.90 MHz, sigma 0.32/f0 -> fc 4.45 MHz, -6 dB BW 62.8 %
% (target 4.5 +- 0.1 MHz, 60 +- 5 %).
P.pulseSigma   = 0.32;      % sigma = pulseSigma / f0
P.pulseCarrier = 4.90e6;    % Hz
P.pulseHalfLen = 5;         % pulse length = pulseHalfLen / carrier

% ---------------------------------------------------------- imaging grid
P.lambda = P.c / P.f0;
P.L      = P.N * P.pitch;     % 42.24 mm
P.lines  = 256;               % compounding / display lines
P.dx     = P.lambda / 2;      % 0.1711 mm
P.zMin   = 5e-3;
P.zMax   = 60e-3;
P.dz     = 0.05e-3;
P.x      = ((0:P.lines-1) - (P.lines-1)/2) * P.dx;
P.z      = (P.zMin:P.dz:P.zMax)';

% ------------------------------------------------------- conventional B-mode
P.txFocus  = 30e-3;                         % m
P.txLines  = 128;
P.txPitch  = 0.33e-3;
P.xf       = ((0:P.txLines-1) - (P.txLines-1)/2) * P.txPitch;
P.txApodW  = P.txFocus / P.Fnum;            % 17.14 mm
P.psfLines = 256;                           % point-target phantom uses 256
P.xpsf     = ((0:P.psfLines-1) - (P.psfLines-1)/2) * P.dx;

% --------------------------------------------------------- plane-wave angles
P.iAng    = -35:35;                          % i in sin(alpha) = i*lambda/L
P.angles  = asin(P.iAng * P.lambda / P.L);
P.nAngles = numel(P.angles);
P.subsetN = [1 5 9 17 33 45 71];
P.nOpt    = P.L / (P.lambda * P.Fnum);       % ~70.5

% ---------------------------------------------------------------- phantoms
P.seed   = 2009;
P.scatN  = 20000;
P.phX    = [-22 22] * 1e-3;    % slightly wider than the display (+-21.8 mm)
P.phZ    = [5 60] * 1e-3;
P.cystR  = 3e-3;
P.cystX  = [-12 0 12] * 1e-3;
P.cystZ  = [20 35 50] * 1e-3;
P.ptX    = [0 0 0 0 0 0 -12 12] * 1e-3;
P.ptZ    = [15 25 30 35 45 55 30 30] * 1e-3;
P.ptAmp  = 30;

% ------------------------------------------------------------------- noise
P.snrPW30 = 20;      % dB target for single plane wave at z = 30 mm
P.nMC     = 20;      % Monte-Carlo repeats for the analytic noise check

% --------------------------------------------------------------- displays
P.drDB   = 60;       % display dynamic range

% ------------------------------------------------------------ field maps
P.fmX = (-20:1:20) * 1e-3;
P.fmZ = (5:1:60) * 1e-3;

% --------------------------------------------------------------- E7 (opt.)
P.e7.warmup = 3;
P.e7.repeat = 10;
P.e7.nList  = [1 5 9 17 33 45 71];

% ------------------------------------------------------- external toolbox
P.fieldIIPath = 'C:\Users\HP\Desktop\~\Field_II_ver_3_30_windows';
end
