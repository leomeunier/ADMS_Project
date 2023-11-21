function [valuenm] = epsilon(x,fs,GrEXP)

    % In x we store all the unknowns (Ai, dampi, omegai, Rh and Rl)
    A = x(1);
    damp = x(2);
    omega = x(3);
    Rhr = x(4);
    Rhi=x(5);
    Rlr=x(6);
    Rli = x(7);
    % Result will store the sum of all number
    result = 0;
    omegafs = 2*pi*fs;
    for n=1:4
        for m = 1:size(fs,2);
            GrNUMnm = A/(-omegafs(m)^2 +1i*2*damp*omega*omegafs(m) + omega^2) + (Rhr+1i*Rhi)/omegafs(m).^2 + Rlr+1i*Rli;
            valuenm = real(GrEXP(n,m)- GrNUMnm)^2 + imag(GrEXP(n,m)- GrNUMnm)^2;
            result = result + valuenm;
        end
    end
