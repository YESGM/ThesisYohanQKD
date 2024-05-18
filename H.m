function [Entropy] = H(x)

    Entropy = -x.*log2(x)-(1-x).*log2(1-x);
end
