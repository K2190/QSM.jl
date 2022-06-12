function ifft3(X)
#IFFT3 Summary of this function goes here
#   Detailed explanation goes here

    X = ifft(ifft(ifft(X, 3),  2),  1)
    return X
end
