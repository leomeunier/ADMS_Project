function [valuenm] = epsilon1(x,omegafs,GrEXP)

    % In x we store all the unknowns (Ai, dampi, omegai, Rh and Rl)
    A = x(1);
    damp = x(2);
    omega = x(3);
    Rhr = x(4);
    Rhi=x(5);
    Rlr=x(6);
    Rli = x(7);
    % Result will store the sum of all number
    % result = 0;
    % omegafs = 2*pi*fs;
    for n = 1
        for m = 1:size(omegafs,2)
            GrNUMnm = A(n)/(-omegafs(m)^2 +1i*2*damp*omega*omegafs(m) + omega^2) + (Rhr(n)+1i*Rhi(n))/omegafs(m).^2 + Rlr(n)+1i*Rli(n);
            valuenm = sum(sum(real(GrEXP(n,m)- GrNUMnm)^2 + imag(GrEXP(n,m)- GrNUMnm)^2));
            % result = result + valuenm;
        end
    end
