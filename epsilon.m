function [result] = epsilon(x,fs,GrEXP)

    % In x we store all the unknowns (Ai, dampi, omegai, Rh and Rl)
    A = x(1:4);
    damp = x(5:8);
    omega = x(9:12);
    Rhr = x(13);
    Rhi = x(14);
    Rlr = x(15);
    Rli = x(16);
    % Result will store the sum of all number
    result = 0;
    omegafs = 2*pi*fs;
    for n = 1:size(A,2)
        for m = 1:size(fs,2)
            GrNUMnm = A(n)/(-omegafs(m)^2 +1i*2*damp(n)*omega(n)*omegafs(m) + omega(n)^2) + Rhr/(omegafs(m)^2) + 1i*Rhi/(omegafs(m)^2) + Rlr + 1i*Rli;
            valuenm = real(GrEXP(n,m)- GrNUMnm)^2 + imag(GrEXP(n,m)- GrNUMnm)^2;
            result = result + valuenm;
        end
    end
