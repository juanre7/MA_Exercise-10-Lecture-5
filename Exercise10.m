% ============================================================
% Exercise 10 – Moving Average FIRs (5-pt and 9-pt)
% tags: DSP, FIR, Moving Average, Frequency/Phase, Group Delay
% Author: Juan Rodriguez Esteban
% Toolbox-free version (no Signal Processing Toolbox required)
% ============================================================
clear; clc; close all;

% -----------------------------
% Coefficients (normalized MAFs)
% -----------------------------
M1 = 5;                 % 5-point moving average
M2 = 9;                 % 9-point moving average
b1 = (1/M1)*ones(1,M1); % DC gain = 1
b2 = (1/M2)*ones(1,M2);
a  = 1;                 % FIR => a = 1 (non-recursive)

% (Optional) Unnormalized versions (same zeros; scaled magnitude)
b1_raw = ones(1,M1);
b2_raw = ones(1,M2);

% -----------------------------
% Test signal (noisy sinusoid)
% -----------------------------
N = 700; n = 0:N-1;
x = sin(2*pi*0.04*n) + 0.25*sin(2*pi*0.28*n) + 0.35*randn(size(n));

% --------------------------------------------
% Time-domain filtering and convolution check
% --------------------------------------------
y1 = filter(b1, a, x);               % causal FIR (MA-5)
y2 = filter(b2, a, x);               % causal FIR (MA-9)
y1_conv = conv(x, b1, 'same');       % should match y1 (FIR == convolution)

% --------------------------------------------
% Optional zero-phase (offline, no group delay)
% --------------------------------------------
if exist('filtfilt','file') == 2
    y1_zero = filtfilt(b1, a, x);
    y2_zero = filtfilt(b2, a, x);
else
    % Toolbox-free forward-backward filtering with edge reflection
    y1_zero = filtfilt_fb(b1, a, x);
    y2_zero = filtfilt_fb(b2, a, x);
end

% -----------------------------
% Frequency responses (mag/dB)
% -----------------------------
nfft = 2048;
[H1, w] = safe_freqz(b1, a, nfft);  H1dB = 20*log10(abs(H1)+eps);
[H2, ~] = safe_freqz(b2, a, nfft);  H2dB = 20*log10(abs(H2)+eps);

% Also show unnormalized to compare amplitude scaling (optional)
[H1raw, ~] = safe_freqz(b1_raw, a, nfft);
[H2raw, ~] = safe_freqz(b2_raw, a, nfft);

% -----------------------------
% Zeros / poles
% -----------------------------
[z1, p1, ~] = safe_tf2zpk(b1_raw, a);   % scaling does not change zeros
[z2, p2, ~] = safe_tf2zpk(b2_raw, a);

% -----------------------------
% Group delay (should be (M-1)/2)
% -----------------------------
[gd1, wg] = safe_grpdelay(b1, a, 1024);
gd2 = safe_grpdelay(b2, a, 1024);
gd1_theory = (M1-1)/2;   % 2.0 samples
gd2_theory = (M2-1)/2;   % 4.0 samples

% ============================================================
%                       P L O T S
% ============================================================

% --- Time domain (causal vs zero-phase) ---
figure('Name','Time-Domain Responses');
subplot(3,1,1);
plot(n, x, 'Color', [0.7 0.7 0.7]); grid on;
title('Input x[n]'); xlabel('n'); ylabel('Amplitude');

subplot(3,1,2);
plot(n, y1, 'LineWidth', 1.25); hold on;
plot(n, y2, 'LineWidth', 1.25);
plot(n, y1_conv, '--', 'LineWidth', 1); % convolution check
grid on; legend('MA-5 (filter)', 'MA-9 (filter)', 'MA-5 (conv, same)', 'Location','best');
title('Causal FIR outputs & convolution equivalence'); xlabel('n');

subplot(3,1,3);
plot(n, y1_zero, 'LineWidth', 1.25); hold on;
plot(n, y2_zero, 'LineWidth', 1.25);
grid on; legend('MA-5 zero-phase', 'MA-9 zero-phase', 'Location','best');
title('Zero-phase (forward-backward): smoothing without delay'); xlabel('n');

% --- Frequency response (linear magnitude) ---
figure('Name','Frequency Responses');
subplot(2,2,1);
plot(w/pi, abs(H1), 'LineWidth', 1.4); hold on;
plot(w/pi, abs(H2), 'LineWidth', 1.4);
grid on; xlabel('\omega/\pi'); ylabel('|H(e^{j\omega})|');
title('Magnitude - normalized MAFs (DC gain = 1)');
legend('MA-5','MA-9','Location','best');

subplot(2,2,2);
plot(w/pi, H1dB, 'LineWidth', 1.25); hold on;
plot(w/pi, H2dB, 'LineWidth', 1.25);
grid on; xlabel('\omega/\pi'); ylabel('Magnitude (dB)');
title('Magnitude in dB - normalized');
legend('MA-5','MA-9','Location','best');

subplot(2,2,3);
plot(w/pi, abs(H1raw), 'LineWidth', 1.0); hold on;
plot(w/pi, abs(H2raw), 'LineWidth', 1.0);
grid on; xlabel('\omega/\pi'); ylabel('|H(e^{j\omega})|');
title('Optional: unnormalized (sum of taps = M)');
legend('MA-5 raw','MA-9 raw','Location','best');

% --- Pole-zero on unit circle ---
figure('Name','Pole-Zero Plots');
subplot(1,2,1); safe_zplane(z1, p1); title('Pole-zero: 5-point MAF');
subplot(1,2,2); safe_zplane(z2, p2); title('Pole-zero: 9-point MAF');

% --- Group delay ---
figure('Name','Group Delay');
plot(wg/pi, gd1, 'LineWidth', 1.25); hold on;
plot(wg/pi, gd2, 'LineWidth', 1.25);
% draw theoretical horizontal lines without yline (for older MATLAB)
plot(wg/pi, gd1_theory*ones(size(wg)), 'k:');
plot(wg/pi, gd2_theory*ones(size(wg)), 'k:');
grid on; xlabel('Normalized Frequency (\times\pi)'); ylabel('Samples');
title('Group Delay of Moving-Average FIRs');
legend('MA-5 GD','MA-9 GD','MA-5 theory','MA-9 theory','Location','best');

% ============================================================
% Helper functions (toolbox-free fallbacks)
% Place them at the end of the script.
% ============================================================

function [H, w] = safe_freqz(b, a, n)
% Use built-in freqz if available; else evaluate H(e^{jw}) directly.
    if exist('freqz','file') == 2
        if nargin < 3
            [H, w] = freqz(b, a);
        else
            [H, w] = freqz(b, a, n);
        end
        return;
    end
    if nargin < 3, n = 512; end
    w = linspace(0, pi, n).';
    kb = 0:numel(b)-1;    Ba = exp(-1j*w*kb) * b(:);
    ka = 0:numel(a)-1;    Aa = exp(-1j*w*ka) * a(:);
    H = Ba ./ Aa;
end

function [gd, w] = safe_grpdelay(b, a, n)
% Use grpdelay if available; else derive from unwrapped phase numerically.
    if exist('grpdelay','file') == 2
        if nargin < 3
            [gd, w] = grpdelay(b, a);
        else
            [gd, w] = grpdelay(b, a, n);
        end
        return;
    end
    if nargin < 3, n = 512; end
    [H, w] = safe_freqz(b, a, n);
    phi = unwrap(angle(H));
    dphi = diff(phi);
    dw = diff(w);
    gd = -[dphi ./ dw; dphi(end) ./ dw(end)];  % pad last value to match length
end

function [z, p, k] = safe_tf2zpk(b, a)
% Use tf2zpk if available; else compute with roots.
    if exist('tf2zpk','file') == 2
        [z, p, k] = tf2zpk(b, a);
        return;
    end
    z = roots(b(:).');   % ensure row
    p = roots(a(:).');
    k = b(1)/a(1);
end

function safe_zplane(z, p)
% Use zplane if available; else draw a simple pole-zero plot.
    if exist('zplane','file') == 2
        zplane(z, p); grid on;
        xlabel('Real Part'); ylabel('Imaginary Part');
        axis equal; xlim([-1.2 1.2]); ylim([-1.2 1.2]);
        return;
    end
    th = linspace(0, 2*pi, 512);
    plot(cos(th), sin(th), 'k:'); hold on; axis equal;
    plot(real(z), imag(z), 'o', 'MarkerSize', 6, 'LineWidth', 1.1);
    plot(real(p), imag(p), 'x', 'MarkerSize', 7, 'LineWidth', 1.1);
    grid on; xlabel('Real Part'); ylabel('Imaginary Part');
    xlim([-1.2 1.2]); ylim([-1.2 1.2]); hold off;
end

function y = filtfilt_fb(b, a, x)
% Zero-phase filtering via forward-backward with edge reflection.
    x = x(:);
    n = length(x);
    ord = max(length(b), length(a));
    npad = max(1, min(3*ord, n-1));
    xpad = [2*x(1)-x((npad+1):-1:2);
            x;
            2*x(end)-x((end-1):-1:(end-npad))];
    y = filter(b, a, xpad);
    y = flipud(y);
    y = filter(b, a, y);
    y = flipud(y);
    y = y(npad+1:npad+n);
end
