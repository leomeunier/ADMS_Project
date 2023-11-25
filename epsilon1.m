function [result] = epsilon(x,fs,GrEXP)

    % In x we store all the unknowns (Ai, dampi, omegai, Rh and Rl)
    A = x(1:4);
    damp = x(5);
    omega = x(6);
    Rhr = x(7);
    Rhi=x(8);
    Rlr=x(9);
    Rli = x(10);
    % Result will store the sum of all number
    result = 0;
    omegafs = 2*pi*fs;
    for n = 1:4
        for m = 1:size(fs,2)
            GrNUMnm = A(n)/(-omegafs(m)^2 +1i*2*damp*omega*omegafs(m) + omega^2) + (Rhr+1i*Rhi)/omegafs(m).^2 + Rlr+1i*Rli;
            valuenm = real(GrEXP(n,m)- GrNUMnm)^2 + imag(GrEXP(n,m)- GrNUMnm)^2;
            result = result + valuenm;
        end
    end
