# EXP 1 A : COMPUTATION OF DFT USING DIRECT AND FFT

# AIM: 

# To Obtain DFT and FFT of a given sequence in SCILAB. 

# APPARATUS REQUIRED: 
PC installed with SCILAB. 

# PROGRAM: 
// DISCRETE FOURIER TRANSFORM 


x = [1, 2, 5, 6];
N = length(x);

// DFT
X_dft = zeros(1, N);
for k = 0:N-1
    sum_val = 0;
    for n = 0:N-1
        W = exp(-%i * 2 * %pi * k * n / N);
        sum_val = sum_val + x(n+1) * W;
    end
    X_dft(k+1) = sum_val;
end

// FFT
pow2N = 2^ceil(log2(N));
if pow2N ~= N then
    x($+1:pow2N) = 0;
    N = pow2N;
end

// Bit-reversal without bitshift
function y = bitrevorder(x)
    N = length(x);
    nbits = log2(N);
    y = zeros(1, N);
    for k = 0:N-1
        bin_str = dec2bin(k, nbits);          // binary string of index
        rev_str = part(bin_str, length(bin_str):-1:1); // reverse string
        rev = bin2dec(rev_str);               // convert back to decimal
        y(rev+1) = x(k+1);
    end
endfunction

X_fft = bitrevorder(x);

for s = 1:log2(N)
    m = 2^s;
    wm = exp(-%i * 2 * %pi / m);
    for k = 0:m:N-1
        for j = 0:(m/2 - 1)
            t = wm^j * X_fft(k + j + m/2 + 1);
            u = X_fft(k + j + 1);
            X_fft(k + j + 1) = u + t;
            X_fft(k + j + m/2 + 1) = u - t;
        end
    end
end

disp("DFT:");
disp(X_dft);
disp("FFT:");
disp(X_fft);


# OUTPUT: 
<img width="1600" height="900" alt="dtsp exp 2 output" src="https://github.com/user-attachments/assets/8cc7c416-2823-4d0a-9c98-d348cbf09bf8" />


# RESULT:
 Thus, the Discrete Fourier Transform and Fast Fourier Transform of the given sequence were
 obtained and its magnitude andphasespectrum were plotted.
