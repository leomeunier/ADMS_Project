function [valuenm] = epsilon1(x,omegafs,GrEXP)

    % In x we store all the unknowns (Ai, dampi, omegai, Rh and Rl)
    A = x(1:4);
    damp = x(5);
    omega = x(6);
    Rhr = x(7:10);
    Rhi=x(11:14);
    Rlr=x(15:18);
    Rli = x(19:22);
    % Result will store the sum of all number
    % result = 0;
    % omegafs = 2*pi*fs;
    for n = 1:4
        for m = 1:size(omegafs,2)
            GrNUMnm = A(n)/(-omegafs(m)^2 +1i*2*damp*omega*omegafs(m) + omega^2) + (Rhr(n)+1i*Rhi(n))/omegafs(m).^2 + Rlr(n)+1i*Rli(n);
            valuenm = sum(sum(real(GrEXP(n,m)- GrNUMnm)^2 + imag(GrEXP(n,m)- GrNUMnm)^2));
            % result = result + valuenm;
        end
    end
