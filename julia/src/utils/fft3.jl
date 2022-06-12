
using FFTW
# using AbstractFFTs  

function fft3(x)
#FFT3 Summary of this function goes here
#   Detailed explanation goes here

    X = fft(fft(fft(x, 3), 2), 1);

    return X
end




# TEST CODE for fft
# A3 = zeros(Float64, 3,3,3)
# A3[:, :, 1] = [1 2 3; 4 5 6; 7 8 9 ]
# A3[:, :, 2] = [1 2 3; 4 5 6; 7 8 9 ]
# A3[:, :, 3] = [1 2 3; 4 5 6; 7 8 9 ]

# F1 = fft(A3,  3)
# println(F1)

# println("---------------------------------")
# F2 = fft(F1, 2)
# println(F2)

# println("---------------------------------")
# F3 = fft(F2, 1)
# println(F3)

# fft11 = fft3(A3)
# println(fft11)